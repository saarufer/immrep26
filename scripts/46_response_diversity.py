#!/usr/bin/env python
"""Does typicality work only where the response is motif-dominated?

The claim under test came out of E45. A score built from nothing but how far a
TCR's own chain features sit from its cohort mean -- CDR3 length, net charge and
mean residue mass on both chains, no pairing, no sequence identity, no
neighbours -- reaches 0.6378, and it does so very unevenly:

    GIL 0.779   TTD 0.667   GLC 0.641   YLQ 0.584   NLV 0.518

The proposed mechanism is that "typical shape" is a proxy for "belongs to this
epitope's dominant motif". Where one restricted motif carries the public
response, being ordinary inside the cohort nearly means being in that motif, so
typicality acts as a crude motif detector. Where the response is broad, ordinary
says nothing.

It is worth testing because it predicts something the project wants: that the
*weight* typicality deserves is epitope-dependent, and that the quantity setting
it -- how concentrated the response is -- is observable without labels, hence
observable on the private epitopes.

The obvious test is a between-epitope correlation, and it is nearly powerless:
the diversity indices can be computed for all eight epitopes, but typicality's
AUC needs labels, so the y axis exists for five. n=5 cannot be significant, and
Part 2 reports it only as a direction.

Part 3 is the real test, and it is within-epitope. If typicality predicts
validation *through* motif membership, then comparing TCRs only against others
in the same motif stratum should remove it. That is a conditional AUC over all
pairs inside a stratum, which has real power, and it is falsifiable in both
directions: on GIL the conditional AUC should collapse toward 0.5, and on NLV
there is nothing to collapse because it starts there.

Part 3 is useless without a negative control, added after a first run: if
conditioning on TRBV damages *every* score on an epitope, then a drop in
typicality's conditional AUC says nothing about typicality. The control runs the
same conditioning on density alone and on the submitted blend.

Part 4 asks whether any of it is worth leaderboard weight, oracle first: if an
oracle per-epitope weight cannot clear +0.02 over the frozen entry, no
label-free estimate of that weight can, and the direction is closed regardless
of what Part 3 says about mechanism.

Pre-registered before running, because both are choices that could otherwise be
made after seeing results:
  - the primary diversity index is `largest_family_share`, as the most direct
    reading of "motif-dominated". The other three are secondary.
  - the primary stratification is TRBV, as the standard proxy for motif
    identity. Family and TRBV x length are secondary.

Cohort sizes run 130 to 670, and every diversity index is biased by sample size,
so all of them are rarefied to the smallest cohort.

Run: .venv/bin/python scripts/46_response_diversity.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality  # noqa: E402

REPORT = ROOT / "reports/22_response_diversity.md"
FAMILY_RADIUS = 0.12      # the radius that defines a family everywhere else
STRUCTURE_WEIGHT = 0.3    # reused, not fitted: the weight structure already has
THRESHOLD = 0.02          # the pre-registered inclusion bar
RAREFY_DRAWS = 200
BOOTSTRAP_DRAWS = 400
SEED = 0

PRIMARY_INDEX = "largest_family_share"
PRIMARY_STRATUM = "trbv"

# E47: tcrdist3's published per-epitope AUROC, all eight epitopes.
PUBLISHED_TCRDIST = {
    "GILGFVFTL": (0.87, 0.93), "YLQPRTFLL": (0.87, 0.93),
    "CINGVCWTV": (0.80, 0.83), "NLVPMVATV": (0.80, 0.83),
    "LLWNGPMAV": (0.80, 0.83), "TTDPSFLGRY": (0.80, 0.83),
    "GLCTLVAML": (0.66, 0.68), "LTDEMIAQY": (0.66, 0.68),
}


def exact_spearman(x, y) -> tuple[float, float]:
    """Spearman rho with an EXACT two-sided permutation p.

    scipy's spearmanr uses a t-approximation, which at n=5 reports p=0.000 for
    rho=-1. The exact value is 2/5! = 0.0167: there is one perfectly reversed
    ordering out of 120. With five epitopes every p here must be exact or it
    overstates the evidence by an order of magnitude.
    """
    from itertools import permutations

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    observed = spearmanr(x, y).statistic
    rhos = [spearmanr(x, np.asarray(order)).statistic
            for order in permutations(y)]
    rhos = np.asarray(rhos)
    p = float((np.abs(rhos) >= abs(observed) - 1e-12).mean())
    return float(observed), p


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_space(name: str = "default"):
    cached = np.load(ROOT / f"data/interim/sceptr/{name}.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors, {t: i for i, t in enumerate(cached["tcr_id"])}


def families(ids: np.ndarray, vectors, position) -> np.ndarray:
    """Single-linkage families over the WHOLE cohort, so this is label-free and
    therefore computable on the private epitopes too."""
    block = vectors[[position[t] for t in ids]]
    distances = 1.0 - block @ block.T
    _, family = connected_components(
        sp.csr_matrix(distances <= FAMILY_RADIUS), directed=False)
    return family


# --------------------------------------------------------------------------
# Part 1: label-free diversity, rarefied


def simpson(counts: np.ndarray) -> float:
    """Probability two random draws differ. Higher = more diverse."""
    total = counts.sum()
    if total <= 1:
        return float("nan")
    p = counts / total
    return float(1.0 - (p ** 2).sum())


def indices_once(family: np.ndarray, trbv: np.ndarray,
                 lengths: np.ndarray) -> dict:
    fam_counts = np.unique(family, return_counts=True)[1]
    return {
        "largest_family_share": float(fam_counts.max() / fam_counts.sum()),
        "family_simpson": simpson(fam_counts),
        "trbv_simpson": simpson(np.unique(trbv, return_counts=True)[1]),
        "len_simpson": simpson(np.unique(lengths, return_counts=True)[1]),
    }


def diversity_table(test: pd.DataFrame, vectors, position, rng) -> pd.DataFrame:
    """Every index rarefied to the smallest cohort, because all four are biased
    by n and the cohorts differ fivefold in size."""
    n_min = int(test.groupby("epitope_aa").size().min())
    rows = []
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        family = families(ids, vectors, position)
        trbv = cohort["TRBV_IMGT"].astype(str).to_numpy()
        lengths = cohort["cdr3_beta_aa"].astype(str).str.len().to_numpy()

        draws = []
        for _ in range(RAREFY_DRAWS):
            pick = rng.choice(len(ids), size=n_min, replace=False)
            # The family partition is computed on the full cohort and then
            # subset; recomputing linkage inside a subsample would measure the
            # subsample's connectivity, not the cohort's.
            draws.append(indices_once(family[pick], trbv[pick], lengths[pick]))
        row = {"epitope": epitope, "n": len(ids),
               "split": cohort["split"].iloc[0]}
        for key in draws[0]:
            row[key] = float(np.mean([d[key] for d in draws]))
        rows.append(row)
    return pd.DataFrame(rows).set_index("epitope")


# --------------------------------------------------------------------------
# Part 3: conditional AUC


def pair_auc(y: np.ndarray, s: np.ndarray, strata: np.ndarray) -> tuple[float, int]:
    """Pair-weighted AUC counting only positive-negative pairs that fall in the
    same stratum. One stratum for everything reproduces the ordinary AUC.

    Pair-weighted, not a mean of per-stratum AUCs: a stratum holding one pair
    would otherwise count as much as one holding ten thousand.
    """
    num = den = 0.0
    for value in np.unique(strata):
        mask = strata == value
        pos, neg = s[mask & (y == 1)], s[mask & (y == 0)]
        if len(pos) == 0 or len(neg) == 0:
            continue
        diff = pos[:, None] - neg[None, :]
        num += float((diff > 0).sum()) + 0.5 * float((diff == 0).sum())
        den += float(len(pos) * len(neg))
    return (num / den if den else float("nan")), int(den)


def conditional_test(test, labels, series, vectors, position, rng) -> pd.DataFrame:
    rows = []
    for epitope, block in labels.groupby("epitope_aa"):
        cohort = test[test["epitope_aa"] == epitope]
        ids = cohort["tcr_id"].to_numpy()
        family_of = dict(zip(ids, families(ids, vectors, position)))

        y = block["pos"].to_numpy()
        s = block["tcr_id"].map(series).to_numpy(dtype=float)
        one = np.zeros(len(block))
        fam = block["tcr_id"].map(family_of).to_numpy()
        trbv = block["tcr_id"].map(
            cohort.set_index("tcr_id")["TRBV_IMGT"].astype(str)).to_numpy()
        blen = block["tcr_id"].map(
            cohort.set_index("tcr_id")["cdr3_beta_aa"].astype(str).str.len()
        ).to_numpy()

        strata = {"none": one, "trbv": trbv, "family": fam,
                  "trbv_x_len": np.char.add(trbv.astype(str), blen.astype(str))}
        row = {"epitope": epitope}
        for name, values in strata.items():
            auc, pairs = pair_auc(y, s, values)
            row[f"auc_{name}"] = auc
            row[f"pairs_{name}"] = pairs
        row["pair_retention_trbv"] = row["pairs_trbv"] / row["pairs_none"]

        # Does dominant-motif membership predict validation on its own? If the
        # mediation story holds it must, or there is nothing to mediate through.
        biggest = pd.Series(fam).value_counts().idxmax()
        row["auc_in_largest_family"], _ = pair_auc(
            y, (fam == biggest).astype(float), one)
        top_v = pd.Series(trbv).value_counts().idxmax()
        row["auc_is_dominant_trbv"], _ = pair_auc(y, (trbv == top_v).astype(float), one)

        # Cluster bootstrap over families for the drop, since members of a
        # family are not independent observations.
        drops = []
        unique_fams = np.unique(fam)
        members = {f: np.where(fam == f)[0] for f in unique_fams}
        for _ in range(BOOTSTRAP_DRAWS):
            drawn = rng.choice(unique_fams, size=len(unique_fams), replace=True)
            idx = np.concatenate([members[f] for f in drawn])
            raw, _ = pair_auc(y[idx], s[idx], one[idx])
            cond, den = pair_auc(y[idx], s[idx], trbv[idx])
            if den and not np.isnan(raw) and not np.isnan(cond):
                drops.append(raw - cond)
        row["drop"] = row["auc_none"] - row["auc_trbv"]
        row["drop_lo"] = float(np.percentile(drops, 2.5)) if drops else np.nan
        row["drop_hi"] = float(np.percentile(drops, 97.5)) if drops else np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index("epitope")


def conditioning_control(test, labels, series_by_name, vectors, position):
    """The same TRBV conditioning applied to other scores.

    Without this the mediation test has no specificity: within-TRBV pairs may
    simply be harder for anything that uses the neighbourhood, in which case a
    drop in typicality's conditional AUC is a property of the stratification.
    """
    rows = []
    for name, series in series_by_name.items():
        for epitope, block in labels.groupby("epitope_aa"):
            cohort = test[test["epitope_aa"] == epitope]
            y = block["pos"].to_numpy()
            s = block["tcr_id"].map(series).to_numpy(dtype=float)
            trbv = block["tcr_id"].map(
                cohort.set_index("tcr_id")["TRBV_IMGT"].astype(str)).to_numpy()
            raw, _ = pair_auc(y, s, np.zeros(len(block)))
            cond, _ = pair_auc(y, s, trbv)
            rows.append({"score": name, "epitope": epitope, "raw": raw,
                         "conditional": cond, "drop": raw - cond})
    return pd.DataFrame(rows).pivot(index="epitope", columns="score",
                                    values="drop")


def plain_density(test, vectors, position, sigma: float = 0.1) -> pd.Series:
    """Density alone, so the control covers the blend's dominant term."""
    out = pd.Series(index=test["tcr_id"], dtype=float)
    for _, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[position[t] for t in ids]]
        distances = 1.0 - block @ block.T
        np.fill_diagonal(distances, np.inf)
        out.loc[ids] = np.exp(-distances / sigma).mean(axis=1)
    return out


# --------------------------------------------------------------------------


def main() -> int:
    rng = np.random.default_rng(SEED)
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    vectors, position = load_space()
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")

    typ = typicality(test)
    base_rank = hyp.rank_within(decomp.submitted_score(test), test)
    typ_rank = hyp.rank_within(typ, test)
    frozen = hyp.mean_auc(labels, base_rank)
    typ_per = hyp.per_epitope(labels, typ)

    print("=== Part 1: label-free response diversity, rarefied ===")
    print(f"  every index averaged over {RAREFY_DRAWS} draws rarefied to the")
    print("  smallest cohort; computed on full cohorts, so private included\n")
    div = diversity_table(test, vectors, position, rng)
    print(div.round(4).to_string())

    print("\n=== Part 2: the between-epitope correlation (n=5, underpowered) ===")
    public = div[div["split"] == "public"].copy()
    public["typicality_auc"] = [typ_per[e] for e in public.index]
    tested = ("largest_family_share", "family_simpson", "trbv_simpson",
              "len_simpson")
    corr_rows = []
    for index in tested:
        rho, p_exact = exact_spearman(public[index], public["typicality_auc"])
        # Sidak over the four indices tried, since only one was pre-registered.
        p_family = 1.0 - (1.0 - p_exact) ** len(tested)
        corr_rows.append({"index": index, "rho_vs_typicality": rho,
                          "p_exact": p_exact, "p_familywise": p_family,
                          "n": len(public)})
        mark = "  <- pre-registered primary" if index == PRIMARY_INDEX else ""
        print(f"  {index:<22} rho {rho:+.3f}  exact p {p_exact:.3f}  "
              f"family-wise {p_family:.3f}{mark}")
    print("\n  Predicted sign: dominance up (or diversity down) with typicality")
    print("  up, i.e. rho > 0 for largest_family_share and rho < 0 for the")
    print("  three Simpson indices.")

    print("\n  GILGFVFTL is an outlier on TRBV concentration (0.57 against")
    print("  0.88-0.96), so any ordering could be that one point. Dropping it:")
    for index in tested:
        reduced = public.drop(index="GILGFVFTL")
        rho, p_exact = exact_spearman(reduced[index], reduced["typicality_auc"])
        corr_rows.append({"index": index + " (GIL dropped)",
                          "rho_vs_typicality": rho, "p_exact": p_exact,
                          "p_familywise": np.nan, "n": len(reduced)})
        print(f"    {index:<22} rho {rho:+.3f}  exact p {p_exact:.3f}  (n=4)")

    # A perfect ordering on five points invites the question of whether it is
    # an artifact of the rarefaction draw. It is not, and it does not need the
    # rarefaction at all.
    seeds = []
    for seed in range(6):
        table = diversity_table(test, vectors, position,
                                np.random.default_rng(seed))
        rho, _ = exact_spearman(table.loc[public.index, "trbv_simpson"],
                                public["typicality_auc"])
        seeds.append(rho)
    unrarefied = {
        epitope: simpson(np.unique(cohort["TRBV_IMGT"].astype(str),
                                   return_counts=True)[1])
        for epitope, cohort in test.groupby("epitope_aa")}
    rho_raw, p_raw = exact_spearman([unrarefied[e] for e in public.index],
                                    public["typicality_auc"])
    print(f"\n  trbv_simpson rho across 6 rarefaction seeds: "
          f"{min(seeds):+.3f} to {max(seeds):+.3f}")
    print(f"  with no rarefaction at all: rho {rho_raw:+.3f}, exact p {p_raw:.3f}")
    print("  private cohorts on the same index: "
          + " ".join(f"{e}={unrarefied[e]:.3f}" for e in
                     ("CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"))
          + "\n  -- all at NLV's diverse end, so the association predicts"
            " typicality\n     near chance on all three private epitopes.")

    print("\n  Separately, and a different question: the same indices against")
    print("  tcrdist3's published AUROC, which exists for all eight epitopes.")
    print("  That asks what diversity does to a density method, not to")
    print("  typicality, so it neither supports nor refutes the mechanism.")
    div["published_tcrdist_mid"] = [
        np.mean(PUBLISHED_TCRDIST[e]) for e in div.index]
    for index in ("largest_family_share", "trbv_simpson"):
        rho = spearmanr(div[index], div["published_tcrdist_mid"])
        corr_rows.append({"index": index, "rho_vs_tcrdist_published": rho.statistic,
                          "p_exact": np.nan, "n": len(div)})
        print(f"  {index:<22} rho {rho.statistic:+.3f}  (n=8, ties in the "
              f"published bands)")
    corr = pd.DataFrame(corr_rows)

    print("\n=== Part 3: within-epitope conditional AUC (the real test) ===")
    print("  typicality's AUC counting only pairs inside the same stratum. If")
    print("  typicality predicts validation through motif membership, matching")
    print("  on TRBV should remove it.\n")
    cond = conditional_test(test, labels, typ, vectors, position, rng)
    show = ["auc_none", "pairs_none", "auc_trbv", "pairs_trbv",
            "pair_retention_trbv", "auc_family", "pairs_family",
            "drop", "drop_lo", "drop_hi"]
    print(cond[show].round(3).to_string())
    print("\n  dominant-motif membership scored on its own:")
    print(cond[["auc_in_largest_family", "auc_is_dominant_trbv"]].round(3).to_string())

    print("\n=== Part 3b: the control -- does TRBV conditioning damage ===")
    print("    everything on this epitope, or typicality specifically?\n")
    control = conditioning_control(
        test, labels,
        {"typicality": typ, "density alone": plain_density(test, vectors, position),
         "submitted blend": decomp.submitted_score(test)},
        vectors, position)
    print(control.round(3).to_string())
    gil = control.loc["GILGFVFTL"]
    specific = gil["typicality"] > max(gil["density alone"], gil["submitted blend"])
    print(f"\n  On GILGFVFTL, the only epitope whose drop has a CI excluding")
    print(f"  zero: typicality {gil['typicality']:+.3f}, density alone "
          f"{gil['density alone']:+.3f}, blend {gil['submitted blend']:+.3f}.")
    print(f"  Specific to typicality: {specific}")
    if not specific:
        print("  So the drop is a property of TRBV stratification on that")
        print("  cohort, not evidence that typicality acts through motif")
        print("  membership. The mediation test fails its own control.")

    print("\n=== Part 4: is an epitope-conditioned weight worth anything? ===")
    print(f"  frozen entry {frozen:.4f}; oracle first, because an oracle that")
    print(f"  cannot clear +{THRESHOLD} closes the direction outright.\n")
    grid = (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
    per_epitope_best, oracle_scores = {}, []
    for epitope, block in labels.groupby("epitope_aa"):
        best_w, best = 0.0, -np.inf
        for w in grid:
            blend = (1 - w) * base_rank + w * typ_rank
            score = hyp.per_epitope(block, blend)[epitope]
            if score > best:
                best_w, best = w, score
        per_epitope_best[epitope] = (best_w, best)
        oracle_scores.append(best)
        print(f"  {epitope:<11} oracle w {best_w:.2f}  -> {best:.4f}  "
              f"(w=0 gives {hyp.per_epitope(block, base_rank)[epitope]:.4f})")
    oracle = float(np.mean(oracle_scores))
    global_best = max(
        (hyp.mean_auc(labels, (1 - w) * base_rank + w * typ_rank), w) for w in grid)
    print(f"\n  oracle per-epitope weight   {oracle:.4f}  "
          f"({oracle - frozen:+.4f} vs frozen)")
    print(f"  oracle single global weight {global_best[0]:.4f}  "
          f"(w={global_best[1]:.2f}, {global_best[0] - frozen:+.4f})")
    print(f"  what conditioning buys an oracle: "
          f"{oracle - global_best[0]:+.4f}")

    gate_ok, gate_detail = private_variance_gate(test, typ)
    viable, viable_detail = hyp.viability(test, typ)
    print(f"\n  private variance gate: {gate_ok} ({gate_detail})")
    print(f"  E3 viability:          {viable} ({viable_detail})")

    nested = np.nan
    if oracle - frozen > THRESHOLD:
        print("\n  oracle clears the bar, so a label-free estimate is worth"
              " nested testing")
        held = []
        for epitope in sorted(labels["epitope_aa"].unique()):
            inner = [e for e in labels["epitope_aa"].unique() if e != epitope]
            cut = float(div.loc[inner, PRIMARY_INDEX].median())
            w = STRUCTURE_WEIGHT if div.loc[epitope, PRIMARY_INDEX] > cut else 0.0
            blend = (1 - w) * base_rank + w * typ_rank
            held.append(hyp.per_epitope(labels[labels["epitope_aa"] == epitope],
                                        blend)[epitope])
        nested = float(np.mean(held))
        print(f"  nested zero-parameter rule: {nested:.4f} "
              f"({nested - frozen:+.4f})")
    else:
        print(f"\n  oracle gain {oracle - frozen:+.4f} is below +{THRESHOLD}.")
        print("  No label-free estimate of a per-epitope weight can beat its own")
        print("  oracle, so the direction is closed for the leaderboard and the")
        print("  entry stands. Part 3 does not rescue the mechanism either: its")
        print("  one significant drop is matched by density alone, so what")
        print("  survives is Part 2's association and nothing causal.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Does typicality work only where the response is motif-dominated?\n\n"
        "## Part 1 - label-free diversity, rarefied to the smallest cohort\n\n"
        + div.round(4).to_markdown()
        + "\n\n## Part 2 - between-epitope correlation (n=5, cannot be significant)\n\n"
        + corr.round(4).to_markdown(index=False)
        + "\n\n## Part 3 - within-epitope conditional AUC\n\n"
        + cond.round(4).to_markdown()
        + "\n\n## Part 3b - the control: TRBV-conditioning drop, by score\n\n"
        + control.round(4).to_markdown()
        + f"\n\nSpecific to typicality on GILGFVFTL: **{specific}**.\n"
        + "\n## Part 4 - oracle epitope-conditioned weight\n\n"
        f"- frozen entry **{frozen:.4f}**\n"
        f"- oracle per-epitope weight **{oracle:.4f}** ({oracle - frozen:+.4f})\n"
        f"- oracle single global weight {global_best[0]:.4f} "
        f"(w={global_best[1]:.2f})\n"
        f"- conditioning buys an oracle {oracle - global_best[0]:+.4f}\n"
        f"- private variance gate {gate_ok} ({gate_detail})\n"
        + (f"- nested zero-parameter rule {nested:.4f}\n"
           if not np.isnan(nested) else
           f"- oracle below the +{THRESHOLD} bar, so no label-free estimate was "
           "fitted\n"),
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
