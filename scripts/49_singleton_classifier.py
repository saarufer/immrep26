#!/usr/bin/env python
"""Can a classifier trained only on singletons rank the TCRs that have no relatives?

The question is the user's, and it is different from every supervised attempt
already settled here in a way that matters. E10 and E13 fitted models on the
whole cohort, where the easy pairs dominate: `one_clustered` pairs have a
within-class AUC of 0.903 and 43.6% of the mass (E50), so almost all of a
fitted model's pressure goes to cases density already handles. Restricting the
training set to singletons replaces the objective with the hard class alone --
and that class is 59.8% of all loss and 61-93% of the private split's pairs.

One trap, stated before the run. "Singleton" at cosine 0.12 does not mean "no
neighbours": E20 measured density at 0.758 *within* singletons and `nn1_dist`
alone at 0.757, because the nearest-neighbour distance still varies and still
ranks. So a singleton-only model is not automatically neighbour-free. Both
variants are run and reported separately: one with every feature, one with
every neighbour-derived feature deleted.

The second problem is that the metric is a per-epitope AUC over *all* rows, so a
separate singleton ranking has to be put back without disturbing the rest. The
splice does that: within an epitope, the singletons keep the rank slots they
already occupy in the base score, and only their ordering among themselves is
replaced. Splicing a score into itself must reproduce it exactly, which is the
anchor checked first.

Features are within-epitope percentile ranks *among singletons*, so every
epitope's feature distribution is uniform on [0,1] and nothing epitope-specific
can leak into a learned coefficient. Selection of features and of the model
happens inside the training folds only.

Run: .venv/bin/python scripts/49_singleton_classifier.py
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
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/25_singleton_classifier.md"
FAMILY_RADIUS = 0.12
SIGMA = 0.1
THRESHOLD = 0.02
SEED = 0

CHARGE = {**{a: 1.0 for a in "KRH"}, **{a: -1.0 for a in "DE"}}
MASS = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103,
        "L": 113, "I": 113, "N": 114, "D": 115, "Q": 128, "K": 128, "E": 129,
        "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186}
# Anything computed from other cohort members. Deleting these is what makes the
# neighbour-free variant neighbour-free.
NEIGHBOUR_FEATURES = ("density", "nn1_dist", "nn3_dist", "frozen")


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def singleton_mask_and_neighbours(test):
    """Label-free, over the whole cohort, so it is available on private rows."""
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    singleton, density, nn1, nn3 = {}, {}, {}, {}
    for _, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, fam = connected_components(
            sp.csr_matrix(d <= FAMILY_RADIUS), directed=False)
        sizes = pd.Series(fam).value_counts()
        clustered = np.isin(fam, sizes[sizes >= 2].index)
        np.fill_diagonal(d, np.inf)
        ordered = np.sort(d, axis=1)
        for i, t in enumerate(ids):
            singleton[t] = not clustered[i]
            density[t] = float(np.exp(-d[i][np.isfinite(d[i])] / SIGMA).mean())
            nn1[t] = -float(ordered[i, 0])
            nn3[t] = -float(ordered[i, :3].mean())
    return singleton, density, nn1, nn3


def build_features(test, singleton, density, nn1, nn3, frozen) -> pd.DataFrame:
    frame = pd.DataFrame(index=test["tcr_id"])
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    for column in boltz.columns:
        # PAE is a distance: lower is better, so negate to make every feature
        # oriented the same way before ranking.
        values = boltz[column].reindex(frame.index)
        frame[column] = -values if column.startswith("pae") else values
    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    for column in pgen.select_dtypes("number").columns:
        frame[f"pgen_{column}"] = pgen[column].reindex(frame.index)
    consistency = pd.read_csv(ROOT / "data/interim/boltz_consistency.csv")
    if "tcr_id" in consistency:
        consistency = consistency.set_index("tcr_id")
        for column in consistency.select_dtypes("number").columns:
            frame[f"cons_{column}"] = consistency[column].reindex(frame.index)
    ordered = test.set_index("tcr_id")
    for side, cdr3 in (("a", "cdr3_alpha_aa"), ("b", "cdr3_beta_aa")):
        seqs = ordered[cdr3].astype(str).reindex(frame.index)
        frame[f"len_{side}"] = seqs.str.len()
        frame[f"charge_{side}"] = [sum(CHARGE.get(c, 0.0) for c in s) for s in seqs]
        frame[f"mass_{side}"] = [np.mean([MASS.get(c, 110) for c in s]) for s in seqs]
    frame["density"] = pd.Series(density).reindex(frame.index)
    frame["nn1_dist"] = pd.Series(nn1).reindex(frame.index)
    frame["nn3_dist"] = pd.Series(nn3).reindex(frame.index)
    frame["frozen"] = frozen.reindex(frame.index)
    return frame.astype(float)


def rank_among_singletons(frame, test, singleton) -> pd.DataFrame:
    """Percentile rank of each feature within its epitope, among singletons only.

    This is what makes a learned coefficient transferable: after it, every
    epitope's singleton feature distribution is uniform on [0, 1], so nothing
    epitope-specific survives into the model.
    """
    out = pd.DataFrame(index=frame.index, columns=frame.columns, dtype=float)
    ids_by_epitope = test.groupby("epitope_aa")["tcr_id"].apply(list)
    for _, ids in ids_by_epitope.items():
        singles = [t for t in ids if singleton[t]]
        if len(singles) < 2:
            continue
        out.loc[singles] = frame.loc[singles].rank(pct=True).to_numpy()
    return out


def splice(base: pd.Series, new: pd.Series, test, singleton) -> pd.Series:
    """Replace the singletons' ordering among themselves, keeping their slots.

    Within an epitope the singletons occupy a fixed set of positions in `base`.
    Those positions are handed out again in the order `new` implies. Clustered
    rows are untouched, so the only thing that can change is the
    both-singleton part of the pairwise AUC.
    """
    out = base.copy()
    for _, cohort in test.groupby("epitope_aa"):
        singles = [t for t in cohort["tcr_id"] if singleton[t]]
        if len(singles) < 2:
            continue
        slots = np.sort(base.loc[singles].to_numpy())
        order = np.argsort(np.argsort(new.loc[singles].to_numpy()))
        out.loc[singles] = slots[order]
    return out


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    frozen = decomp.submitted_score(test)
    epitopes = sorted(labels["epitope_aa"].unique())

    singleton, density, nn1, nn3 = singleton_mask_and_neighbours(test)
    raw = build_features(test, singleton, density, nn1, nn3, frozen)
    features = rank_among_singletons(raw, test, singleton)

    labelled = labels.copy()
    labelled["single"] = labelled["tcr_id"].map(singleton)
    singles = labelled[labelled["single"]].reset_index(drop=True)
    print("=== the training set ===")
    for epitope, block in singles.groupby("epitope_aa"):
        print(f"  {epitope:<11} {len(block):4d} labelled singletons, "
              f"{int(block['pos'].sum()):3d} positive "
              f"({block['pos'].mean():.1%})")
    print(f"  total      {len(singles):4d}")

    print("\n=== anchor: splicing a score into itself must change nothing ===")
    identity = splice(frozen, frozen, test, singleton)
    base_macro = hyp.mean_auc(labels, frozen)
    print(f"  frozen {base_macro:.6f}   spliced into itself "
          f"{hyp.mean_auc(labels, identity):.6f}")
    assert abs(hyp.mean_auc(labels, identity) - base_macro) < 1e-9, "splice is lossy"

    def singleton_auc(series) -> tuple[float, list]:
        per = []
        for epitope, block in singles.groupby("epitope_aa"):
            values = block["tcr_id"].map(series).to_numpy(dtype=float)
            per.append(roc_auc_score(block["pos"], values)
                       if len(np.unique(values)) > 1 else np.nan)
        return float(np.nanmean(per)), per

    print("\n=== Part 1: single features, within labelled singletons ===")
    print("  E20 put density at 0.758 here. Anything a classifier learns has to")
    print("  beat the best single feature, or it is not earning its parameters.\n")
    singles_rows = []
    for column in features.columns:
        value, per = singleton_auc(features[column])
        singles_rows.append({"feature": column, "mean_auc": value,
                             "worst_epitope": float(np.nanmin(per)),
                             "neighbour_free": column not in NEIGHBOUR_FEATURES})
    ladder = pd.DataFrame(singles_rows).sort_values("mean_auc", ascending=False)
    print(ladder.head(12).round(4).to_string(index=False))
    free = ladder[ladder["neighbour_free"]]
    print(f"\n  best overall         {ladder.iloc[0]['feature']} "
          f"{ladder.iloc[0]['mean_auc']:.4f}")
    print(f"  best neighbour-free  {free.iloc[0]['feature']} "
          f"{free.iloc[0]['mean_auc']:.4f}")

    print("\n=== Part 2: classifiers trained on singletons only, nested LOEO ===")
    print("  Fit on four epitopes' singletons, predict the fifth's. Feature")
    print("  selection happens inside the training folds.\n")
    variants = {
        "all features": list(features.columns),
        "neighbour-free only": [c for c in features.columns
                                if c not in NEIGHBOUR_FEATURES],
    }
    models = {
        "logistic L2": lambda: LogisticRegression(max_iter=2000, C=0.1),
        "logistic L2 (C=1)": lambda: LogisticRegression(max_iter=2000, C=1.0),
        "boosting depth 2": lambda: HistGradientBoostingClassifier(
            max_depth=2, max_iter=120, learning_rate=0.06,
            l2_regularization=1.0, random_state=SEED),
    }
    results = []
    for variant, columns in variants.items():
        for model_name, factory in models.items():
            for top_k in (5, 15, None):
                held, spliced_scores = [], frozen.copy()
                for epitope in epitopes:
                    train = singles[singles["epitope_aa"] != epitope]
                    held_block = singles[singles["epitope_aa"] == epitope]
                    use = columns
                    if top_k is not None:
                        inner = []
                        for column in columns:
                            values = train["tcr_id"].map(features[column]).to_numpy()
                            inner.append((roc_auc_score(train["pos"], values)
                                          if len(np.unique(values)) > 1 else 0.5,
                                          column))
                        use = [c for _, c in sorted(inner, reverse=True)[:top_k]]
                    x = features.loc[train["tcr_id"], use].to_numpy()
                    y = train["pos"].to_numpy()
                    model = factory().fit(np.nan_to_num(x, nan=0.5), y)
                    xo = features.loc[held_block["tcr_id"], use].to_numpy()
                    p = model.predict_proba(np.nan_to_num(xo, nan=0.5))[:, 1]
                    held.append(roc_auc_score(held_block["pos"], p)
                                if len(np.unique(p)) > 1 else np.nan)
                    spliced_scores.loc[held_block["tcr_id"]] = p
                # Splice the held-out predictions back, epitope by epitope.
                combined = splice(frozen, spliced_scores, test, singleton)
                results.append({
                    "features": variant, "model": model_name,
                    "top_k": "all" if top_k is None else top_k,
                    "singleton_auc": float(np.nanmean(held)),
                    "worst_epitope": float(np.nanmin(held)),
                    "spliced_macro": hyp.mean_auc(labels, combined),
                })
                print(f"  {variant:<20} {model_name:<18} top{str(top_k or 'all'):>4}  "
                      f"singleton {np.nanmean(held):.4f}  "
                      f"spliced macro {results[-1]['spliced_macro']:.4f}")

    table = pd.DataFrame(results).sort_values("spliced_macro", ascending=False)
    frozen_singleton, _ = singleton_auc(frozen)
    print(f"\n  the frozen entry within singletons: {frozen_singleton:.4f}")
    print(f"  the frozen entry's macro:           {base_macro:.4f}")
    best = table.iloc[0]
    print(f"\n  best spliced macro {best['spliced_macro']:.4f} "
          f"({best['spliced_macro'] - base_macro:+.4f}) from "
          f"{best['features']} / {best['model']} / top {best['top_k']}")
    print(f"  its within-singleton AUC {best['singleton_auc']:.4f} against the "
          f"frozen {frozen_singleton:.4f} ({best['singleton_auc'] - frozen_singleton:+.4f})")

    print("\n=== verdict against the pre-registered bar ===")
    gain = best["spliced_macro"] - base_macro
    print(f"  gain {gain:+.4f} against +{THRESHOLD}: "
          f"{'CLEARS' if gain > THRESHOLD else 'does not clear'}")
    if best["singleton_auc"] <= frozen_singleton:
        print("  And the classifier does not even beat the frozen score inside")
        print("  the class it was built for, so the restriction of the training")
        print("  set did not buy what it was supposed to buy.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Classifying singletons only\n\n"
        f"Labelled singletons: {len(singles)} across five epitopes. The splice "
        "reproduces the frozen score exactly when fed itself.\n\n"
        "## Single features within labelled singletons\n\n"
        + ladder.head(20).round(4).to_markdown(index=False)
        + "\n\n## Singleton-only classifiers, nested LOEO\n\n"
        + table.round(4).to_markdown(index=False)
        + f"\n\n- frozen entry within singletons **{frozen_singleton:.4f}**\n"
        f"- frozen entry macro **{base_macro:.4f}**\n"
        f"- best spliced macro **{best['spliced_macro']:.4f}** ({gain:+.4f})\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
