#!/usr/bin/env python
"""The one learning formulation never tried: weights conditioned on the pair.

Three things get called pairwise learning and two are already settled here.
E16 put a pairwise ranking surrogate under gradient descent to learn the metric
density lives in, and it failed structurally -- the surrogate is *negatively*
correlated with AUC across the sigma grid and only aligns as its temperature
goes to zero, where it becomes the 1-AUC step function and the gradients vanish.
E10 and E13 fitted pointwise learners on the features and lost to the
zero-parameter rank average. What was never built is a model whose training
example is a *pair*, using quantities that exist only for a pair.

The motivation is E22's decomposition. Error rate varies enormously by the
cluster state of a pair -- 0.207 when neither member has relatives against
0.673 when the positive is isolated and the negative is clustered -- and 62% of
all loss sits in the first class. A pointwise score cannot express "when both
are isolated, trust density less"; a pair-conditioned one can.

**Only four of E22's five classes are observable.** Separating
`pos_clu_neg_sing` from `pos_clu_neg_sing` requires knowing which member is the
positive, which is the label. Their union -- exactly one member clustered -- is
label-free, and so are the other three. Conditioning on the five would be
leakage, and the whole point is a rule that works on the private split, where
nothing is labelled.

The model is deliberately the smallest thing that can express the hypothesis:

    decision(i, j) = sum_k w[k, class(i,j)] * (f_k(i) - f_k(j))

over the components the frozen entry already uses, plus typicality. Every term
is an antisymmetric difference scaled by a symmetric pair property, so
P(i beats j) = 1 - P(j beats i) holds by construction, and w_struct = 0.3,
w_typ = 0 reproduces the frozen entry exactly -- which is the sanity anchor.

Because class partitions the pairs, each class's weights affect only its own
pairs, so the oracle separates into four independent searches and Part 1 can
price the ceiling of the entire idea cheaply. If an oracle that is told the
best weight per class cannot clear +0.02, no fitted version can, and nothing
needs to be fitted.

The trap to avoid is the effective sample size. GILGFVFTL offers 73,914 pairs,
but they come from 670 TCRs in roughly 400 families, and the independent unit is
the family, not the pair. This project has already been bitten by exactly that
once: E41's permutation null was computed on family-collapsed labels while its
observed statistic was computed at site level, and correcting it in E46 moved
+3.629 to -0.149. So Part 2's training objective weights each pair by the
inverse size of its (family_i, family_j) cell, making the effective count the
number of family pairs.

Run: .venv/bin/python scripts/47_pairwise_learning.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality  # noqa: E402

REPORT = ROOT / "reports/23_pairwise_learning.md"
FAMILY_RADIUS = 0.12
STRUCTURE_WEIGHT = 0.3     # the frozen entry's weight, the anchor at w_typ = 0
THRESHOLD = 0.02
STEP = 0.05
SEED = 0

CLASSES = ("both_singleton", "one_clustered", "same_cluster", "diff_cluster")


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def simplex(step: float = STEP):
    """Weights on (structure, typicality); density takes the remainder."""
    values = np.round(np.arange(0.0, 1.0 + 1e-9, step), 4)
    return [(s, t) for s, t in product(values, values) if s + t <= 1.0 + 1e-9]


def cohort_families(test: pd.DataFrame, vectors, position) -> dict:
    """Label-free, over the whole cohort, so the same partition is available on
    the private epitopes."""
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[position[t] for t in ids]]
        distances = 1.0 - block @ block.T
        _, family = connected_components(
            sp.csr_matrix(distances <= FAMILY_RADIUS), directed=False)
        sizes = pd.Series(family).value_counts()
        out[epitope] = {
            "family": dict(zip(ids, family)),
            "clustered": dict(zip(ids, np.isin(family, sizes[sizes >= 2].index))),
            "distance": pd.DataFrame(distances, index=ids, columns=ids),
        }
    return out


def build_pairs(labels, components, structure_info) -> dict:
    """Every positive-negative pair inside an epitope, with its antisymmetric
    feature differences, its label-free class, and its family cell."""
    out = {}
    for epitope, block in labels.groupby("epitope_aa"):
        info = structure_info[epitope]
        ids = block["tcr_id"].to_numpy()
        y = block["pos"].to_numpy()
        pos_ids, neg_ids = ids[y == 1], ids[y == 0]

        feats = np.stack([block["tcr_id"].map(c).to_numpy(dtype=float)
                          for c in components])      # (n_features, n_rows)
        index = {t: i for i, t in enumerate(ids)}
        pi = np.array([index[t] for t in pos_ids])
        ni = np.array([index[t] for t in neg_ids])

        # (n_features, n_pos, n_neg) differences; positive minus negative.
        delta = feats[:, pi][:, :, None] - feats[:, ni][:, None, :]

        clustered = np.array([info["clustered"][t] for t in ids])
        pc, nc = clustered[pi][:, None], clustered[ni][None, :]
        family = np.array([info["family"][t] for t in ids])
        same = (family[pi][:, None] == family[ni][None, :]) & pc & nc

        klass = np.empty((len(pi), len(ni)), dtype=object)
        klass[~pc & ~nc] = "both_singleton"
        klass[pc ^ nc] = "one_clustered"
        klass[pc & nc & same] = "same_cluster"
        klass[pc & nc & ~same] = "diff_cluster"

        dist = info["distance"].loc[pos_ids, neg_ids].to_numpy()
        cell = np.char.add(np.char.add(family[pi][:, None].astype(str), "|"),
                           family[ni][None, :].astype(str))
        out[epitope] = {"delta": delta, "klass": klass, "distance": dist,
                        "cell": cell, "n_pairs": len(pi) * len(ni)}
    return out


def class_contribution(pairs: dict, klass: str, weights, epitopes,
                       cell_weight: bool = False) -> float:
    """Mean over epitopes of (this class's concordant pairs) / (all its pairs).

    Summing this over the four classes gives the mean per-epitope AUC, because
    the classes partition the pairs and every epitope's denominator is its own
    total pair count. That separability is what makes the oracle cheap.
    """
    w_struct, w_typ = weights
    w = np.array([1.0 - w_struct - w_typ, w_struct, w_typ])
    total = 0.0
    for epitope in epitopes:
        data = pairs[epitope]
        mask = data["klass"] == klass
        if not mask.any():
            continue
        decision = np.tensordot(w, data["delta"], axes=(0, 0))[mask]
        credit = (decision > 0).astype(float) + 0.5 * (decision == 0)
        if cell_weight:
            # One effective observation per (family_i, family_j) cell.
            cells = data["cell"][mask]
            _, inverse, counts = np.unique(cells, return_inverse=True,
                                           return_counts=True)
            weightvec = 1.0 / counts[inverse]
            total += float((credit * weightvec).sum() / weightvec.sum()
                           * mask.sum() / data["n_pairs"])
        else:
            total += float(credit.sum() / data["n_pairs"])
    return total / len(epitopes)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}

    # The frozen entry's own two components, plus typicality as a third.
    frozen = decomp.submitted_score(test)
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure_rank = hyp.rank_within(boltz["pae_tcr_pep_min"], test)
    # density = (frozen - 0.3 * structure) / 0.7, exactly as it was built.
    density_rank = (frozen - STRUCTURE_WEIGHT * structure_rank) / (1 - STRUCTURE_WEIGHT)
    typ_rank = hyp.rank_within(typicality(test), test)
    components = [density_rank, structure_rank, typ_rank]

    info = cohort_families(test, vectors, position)
    pairs = build_pairs(labels, components, info)
    epitopes = sorted(labels["epitope_aa"].unique())

    anchor = (STRUCTURE_WEIGHT, 0.0)
    frozen_auc = sum(class_contribution(pairs, k, anchor, epitopes) for k in CLASSES)
    direct = hyp.mean_auc(labels, frozen)
    print("=== anchor: does the decomposition reproduce the frozen entry? ===")
    print(f"  summed over the four classes {frozen_auc:.4f}")
    print(f"  scored directly              {direct:.4f}")
    print(f"  difference                   {abs(frozen_auc - direct):.6f}")
    assert abs(frozen_auc - direct) < 2e-3, "the pair decomposition is not the entry"

    print("\n=== Part 1: the share of each label-free class ===")
    rows = []
    for klass in CLASSES:
        share = float(np.mean([
            (pairs[e]["klass"] == klass).sum() / pairs[e]["n_pairs"]
            for e in epitopes]))
        at_anchor = class_contribution(pairs, klass, anchor, epitopes)
        within = at_anchor / share if share else np.nan
        rows.append({"class": klass, "pair_share": share,
                     "contribution_at_frozen_weights": at_anchor,
                     "auc_within_class": within,
                     "loss": share * (1.0 - within)})
        print(f"  {klass:<15} {share:6.1%} of pairs   within-class AUC "
              f"{within:.3f}")
    total_loss = sum(r["loss"] for r in rows)
    for row in rows:
        row["loss_share"] = row["loss"] / total_loss
    print("\n  where the loss actually sits:")
    for row in rows:
        print(f"  {row['class']:<15} {row['loss_share']:6.1%} of all loss")

    print("\n=== Part 2: the oracle, one weight pair per class ===")
    print("  Each class is searched independently, which is exact: the classes")
    print("  partition the pairs, so a class's weights move only its own.\n")
    oracle_total, oracle_rows = 0.0, []
    for klass in CLASSES:
        best, best_w = -np.inf, None
        for weights in simplex():
            value = class_contribution(pairs, klass, weights, epitopes)
            if value > best:
                best, best_w = value, weights
        oracle_total += best
        at_anchor = class_contribution(pairs, klass, anchor, epitopes)
        share = rows[CLASSES.index(klass)]["pair_share"]
        oracle_rows.append({
            "class": klass, "oracle_w_structure": best_w[0],
            "oracle_w_typicality": best_w[1],
            "oracle_w_density": round(1 - sum(best_w), 4),
            "gain": best - at_anchor,
            "within_class_auc": best / share if share else np.nan})
        print(f"  {klass:<15} density {1 - sum(best_w):.2f}  structure "
              f"{best_w[0]:.2f}  typicality {best_w[1]:.2f}   "
              f"gain {best - at_anchor:+.4f}")
    print(f"\n  oracle class-conditioned weights {oracle_total:.4f}  "
          f"({oracle_total - frozen_auc:+.4f} vs frozen)")

    best_global, best_gw = -np.inf, None
    for weights in simplex():
        value = sum(class_contribution(pairs, k, weights, epitopes) for k in CLASSES)
        if value > best_global:
            best_global, best_gw = value, weights
    print(f"  oracle single global weight      {best_global:.4f}  "
          f"(structure {best_gw[0]:.2f}, typicality {best_gw[1]:.2f}, "
          f"{best_global - frozen_auc:+.4f})")
    print(f"  what pair-conditioning buys an oracle: "
          f"{oracle_total - best_global:+.4f}")

    nested = np.nan
    if oracle_total - frozen_auc > THRESHOLD:
        print("\n=== Part 3: nested LOEO, family-cell weighted ===")
        held = []
        for epitope in epitopes:
            inner = [e for e in epitopes if e != epitope]
            chosen = {}
            for klass in CLASSES:
                best, best_w = -np.inf, anchor
                for weights in simplex():
                    value = class_contribution(pairs, klass, weights, inner,
                                               cell_weight=True)
                    if value > best:
                        best, best_w = value, weights
                chosen[klass] = best_w
            held.append(sum(class_contribution(pairs, k, chosen[k], [epitope])
                            for k in CLASSES))
            print(f"  {epitope:<11} held-out {held[-1]:.4f}   "
                  + " ".join(f"{k[:4]}:s{chosen[k][0]:.2f}/t{chosen[k][1]:.2f}"
                             for k in CLASSES))
        nested = float(np.mean(held))
        print(f"\n  nested LOEO {nested:.4f} ({nested - frozen_auc:+.4f}) "
              f"against +{THRESHOLD}")
    else:
        print(f"\n  Oracle gain {oracle_total - frozen_auc:+.4f} is below "
              f"+{THRESHOLD}.")
        print("  An oracle told the best weights per class is the ceiling of")
        print("  every fitted version, so nothing was fitted.")
        print(f"\n  And the sharper comparison: E15 measured the cost of")
        print(f"  cross-validated selection at -0.0048 over a search of similar")
        print(f"  size ({len(CLASSES)} classes x {len(simplex())} weight pairs")
        print(f"  here against 323 configurations there). The entire oracle")
        print(f"  ceiling, {oracle_total - frozen_auc:+.4f}, is smaller than the")
        print("  price of choosing within it, so a fitted version is expected to")
        print("  come out net negative, not merely short of the bar.")
        biggest = max(rows, key=lambda r: r["loss_share"])
        gain_there = next(r["gain"] for r in oracle_rows
                          if r["class"] == biggest["class"])
        print(f"\n  The motivation was E22's largest loss bucket. Here "
              f"{biggest['class']}")
        print(f"  holds {biggest['loss_share']:.1%} of all loss, and the oracle "
              f"moves it by {gain_there:+.4f}:")
        print("  reweighting cannot reach the place the loss lives.")

    print("\n=== the finer conditioning, for the record ===")
    print("  Mutual distance is the one genuinely pairwise continuous quantity")
    print("  with no pointwise counterpart. Within both_singleton, does the")
    print("  oracle weight move with it?\n")
    tercile_rows = []
    for which, name in enumerate(("near", "middle", "far")):
        best, best_w = -np.inf, None
        for weights in simplex(0.1):
            total = 0.0
            for epitope in epitopes:
                data = pairs[epitope]
                mask = data["klass"] == "both_singleton"
                if not mask.any():
                    continue
                d = data["distance"]
                cuts = np.quantile(d[mask], [1 / 3, 2 / 3])
                band = (mask & (d <= cuts[0]) if which == 0 else
                        mask & (d > cuts[1]) if which == 2 else
                        mask & (d > cuts[0]) & (d <= cuts[1]))
                if not band.any():
                    continue
                w_s, w_t = weights
                w = np.array([1 - w_s - w_t, w_s, w_t])
                decision = np.tensordot(w, data["delta"], axes=(0, 0))[band]
                credit = (decision > 0).astype(float) + 0.5 * (decision == 0)
                total += float(credit.mean())
            if total > best:
                best, best_w = total, weights
        tercile_rows.append({"distance_tercile": name,
                             "oracle_w_structure": best_w[0],
                             "oracle_w_typicality": best_w[1],
                             "mean_within_band_auc": best / len(epitopes)})
        print(f"  {name:<7} oracle structure {best_w[0]:.2f}  typicality "
              f"{best_w[1]:.2f}   within-band AUC {best / len(epitopes):.3f}")

    gate_ok, gate_detail = private_variance_gate(test, frozen)
    print(f"\n  private variance gate on the entry: {gate_ok} ({gate_detail})")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Pair-conditioned weights: the pairwise formulation never tried\n\n"
        f"Anchor: the pair decomposition reproduces the frozen entry to "
        f"{abs(frozen_auc - direct):.6f}.\n\n"
        "## Label-free classes\n\n"
        + pd.DataFrame(rows).round(4).to_markdown(index=False)
        + "\n\n## Oracle, one weight pair per class\n\n"
        + pd.DataFrame(oracle_rows).round(4).to_markdown(index=False)
        + f"\n\n- frozen entry **{frozen_auc:.4f}**\n"
        f"- oracle class-conditioned **{oracle_total:.4f}** "
        f"({oracle_total - frozen_auc:+.4f})\n"
        f"- oracle single global weight {best_global:.4f} "
        f"(structure {best_gw[0]:.2f}, typicality {best_gw[1]:.2f})\n"
        f"- pair-conditioning buys an oracle {oracle_total - best_global:+.4f}\n"
        + (f"- nested LOEO **{nested:.4f}**\n" if not np.isnan(nested) else
           f"- oracle below +{THRESHOLD}, so nothing was fitted\n")
        + "\n## Mutual distance within both_singleton\n\n"
        + pd.DataFrame(tercile_rows).round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
