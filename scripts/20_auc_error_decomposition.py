#!/usr/bin/env python
"""Decompose 1 - AUC into mechanism classes, so effort goes where the error is.

Motivation (a reviewer's P0 proposal). Mean per-epitope AUC is exactly the
probability that a random validated TCR outranks a random non-validating one
from the same cohort, ties counting half. So the remaining 0.15 of error is a
set of specific positive-negative pairs that come out in the wrong order, and
those pairs can be grouped by what distinguishes the two TCRs.

The inference this guards against: 54% of positives are singletons, but that
does **not** imply 54% of the ranking error involves them. A singleton positive
can still beat a singleton negative, while a handful of negatives sitting inside
dense positive clusters can generate a large number of inverted comparisons on
their own. Only the decomposition settles which.

Partition of every positive-negative pair within an epitope, by cluster state
(single-linkage at SCEPTR cosine <= 0.12). It is mutually exclusive and
complete, so the class losses sum to 1 - AUC, which the script asserts:

  same_cluster       both clustered, in the same cluster  -> fine specificity
  diff_cluster       both clustered, different clusters
  pos_clu_neg_sing   positive clustered, negative isolated -> the easy case
  pos_sing_neg_clu   positive isolated, negative clustered -> adversarial
  both_singleton     neither clustered                     -> no neighbourhood

Run: python scripts/20_auc_error_decomposition.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

CLUSTER_RADIUS = 0.12
SIGMA = 0.1
STRUCTURE_WEIGHT = 0.3
OUT = ROOT / "reports/06_auc_error_decomposition.md"

CLASSES = ("pos_clu_neg_sing", "diff_cluster", "same_cluster", "both_singleton", "pos_sing_neg_clu")


def submitted_score(test: pd.DataFrame) -> pd.Series:
    """Reproduce the submitted scoring function, so the decomposition describes
    the model actually entered rather than a simplified stand-in."""
    def rank(values: pd.Series) -> pd.Series:
        frame = test[["tcr_id", "epitope_aa"]].copy()
        frame["v"] = frame["tcr_id"].map(values)
        return frame.assign(r=frame.groupby("epitope_aa")["v"].rank(pct=True)).set_index("tcr_id")["r"]

    components = []
    for space in ("default", "small"):
        cached = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
        pos = {t: i for i, t in enumerate(cached["tcr_id"])}
        vectors = cached["vectors"]
        values = pd.Series(index=test["tcr_id"], dtype=float)
        for _, cohort in test.groupby("epitope_aa"):
            block = vectors[[pos[t] for t in cohort["tcr_id"]]]
            block = block / np.linalg.norm(block, axis=1, keepdims=True)
            d = 1.0 - block @ block.T
            np.fill_diagonal(d, np.inf)
            values.loc[cohort["tcr_id"]] = np.exp(-d / SIGMA).mean(axis=1)
        components.append(rank(values))

    values = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
        pos = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [pos[t] for t in cohort["tcr_id"]]
        d = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
        np.fill_diagonal(d, np.inf)
        values.loc[cohort["tcr_id"]] = np.exp(-d / 20.0).mean(axis=1)
    components.append(rank(values))

    density = sum(components) / len(components)
    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    return (1 - STRUCTURE_WEIGHT) * density + STRUCTURE_WEIGHT * rank(structure["pae_tcr_pep_min"])


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scores = submitted_score(test)

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    normed = cached["vectors"] / np.linalg.norm(cached["vectors"], axis=1, keepdims=True)

    rows = []
    for epitope, cohort in labels.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = normed[[position[t] for t in ids]]
        distances = 1.0 - block @ block.T
        _, cluster = connected_components(
            sp.csr_matrix(distances <= CLUSTER_RADIUS), directed=False
        )
        sizes = pd.Series(cluster).value_counts()
        clustered = np.isin(cluster, sizes[sizes >= 2].index)

        y = cohort["pos"].to_numpy()
        s = cohort["tcr_id"].map(scores).to_numpy()
        pos_idx, neg_idx = np.where(y == 1)[0], np.where(y == 0)[0]
        total_pairs = len(pos_idx) * len(neg_idx)

        # Every positive-negative pair, classified and scored.
        sp_s, sn_s = s[pos_idx][:, None], s[neg_idx][None, :]
        wrong = (sp_s < sn_s).astype(float) + 0.5 * (sp_s == sn_s)

        pc, nc = clustered[pos_idx][:, None], clustered[neg_idx][None, :]
        same = (cluster[pos_idx][:, None] == cluster[neg_idx][None, :]) & pc & nc
        klass = np.full(wrong.shape, "", dtype=object)
        klass[pc & nc & same] = "same_cluster"
        klass[pc & nc & ~same] = "diff_cluster"
        klass[pc & ~nc] = "pos_clu_neg_sing"
        klass[~pc & nc] = "pos_sing_neg_clu"
        klass[~pc & ~nc] = "both_singleton"

        record = {"epitope": epitope, "auc": roc_auc_score(y, s), "pairs": total_pairs}
        for name in CLASSES:
            mask = klass == name
            record[f"{name}__share"] = mask.sum() / total_pairs
            record[f"{name}__loss"] = wrong[mask].sum() / total_pairs
            record[f"{name}__err"] = wrong[mask].mean() if mask.sum() else np.nan
        rows.append(record)

    frame = pd.DataFrame(rows)

    # Correctness check: a complete, exclusive partition must sum to 1 - AUC.
    reconstructed = frame[[f"{c}__loss" for c in CLASSES]].sum(axis=1)
    residual = (reconstructed - (1 - frame["auc"])).abs().max()
    if residual > 1e-9:
        raise AssertionError(f"partition does not reconstruct 1-AUC; residual {residual:.2e}")
    print(f"partition reconstructs 1 - AUC exactly (max residual {residual:.2e})\n")

    pd.set_option("display.width", 260)
    print("=== share of all positive-negative pairs falling in each class ===")
    share = frame.set_index("epitope")[[f"{c}__share" for c in CLASSES]]
    share.columns = CLASSES
    share.loc["MEAN"] = share.mean()
    print((share * 100).round(1).to_string())

    print("\n=== AUC loss contributed by each class (sums to 1 - AUC) ===")
    loss = frame.set_index("epitope")[[f"{c}__loss" for c in CLASSES]]
    loss.columns = CLASSES
    loss["1-AUC"] = 1 - frame.set_index("epitope")["auc"]
    loss.loc["MEAN"] = loss.mean()
    print(loss.round(4).to_string())

    print("\n=== error rate *within* each class (how hard is a pair of this kind) ===")
    err = frame.set_index("epitope")[[f"{c}__err" for c in CLASSES]]
    err.columns = CLASSES
    err.loc["MEAN"] = err.mean()
    print(err.round(3).to_string())

    mean_loss = loss.loc["MEAN", list(CLASSES)].sort_values(ascending=False)
    mean_share = share.loc["MEAN", list(CLASSES)]
    print("\n=== where the error actually lives ===")
    for name, value in mean_loss.items():
        print(
            f"  {name:<18} {value / mean_loss.sum() * 100:5.1f}% of all loss   "
            f"(on {mean_share[name] * 100:5.1f}% of pairs)"
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Pairwise AUC error decomposition\n\n"
        "Mean per-epitope AUC is the probability a validated TCR outranks a\n"
        "non-validating one from the same cohort, so 1 - AUC is a set of inverted\n"
        "pairs. Grouping them by cluster state gives a mutually exclusive, complete\n"
        "partition whose losses sum to 1 - AUC (asserted in the script).\n\n"
        "## Share of pairs\n\n" + (share * 100).round(1).to_markdown()
        + "\n\n## AUC loss contributed\n\n" + loss.round(4).to_markdown()
        + "\n\n## Error rate within class\n\n" + err.round(3).to_markdown() + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
