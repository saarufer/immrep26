#!/usr/bin/env python
"""Two checks on the headline estimate itself, rather than on any new feature.

**S06 -- leave-one-clonotype-family-out.** Leave-one-epitope-out models a new
pMHC, but it does not show that the score avoids leaning on clonotype families
that recur across the dataset. If performance depends on a family being
represented more than once, then the private epitopes -- where almost every
clonotype is a singleton (E27: record:clonotype 1.000) -- are a harder setting
than LOEO suggests, and 0.846 is optimistic. This is not comparable to a LOEO
number: it answers a different question, which is how much of the ranking
survives when each family is scored against a cohort stripped of its relatives.

**D03 -- alpha/beta discordance.** The alpha chain carries more of this signal
than beta (E6: 0.799 against 0.751), which inverts the usual expectation and has
no explanation yet. One candidate is alpha/beta mispairing in the source data: a
TCR whose beta looks typical of the epitope while its alpha looks random would
show a large rank disagreement between single-chain densities. That is testable
as a feature, and E3's rule applies -- it must vary within each private epitope.

Run: python scripts/25_family_holdout_and_chain_discordance.py
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
from immrep.metric import per_epitope_auc  # noqa: E402

SIGMA = 0.1
FAMILY_RADIUS = 0.12   # same radius that defines a cluster elsewhere
STRUCTURE_WEIGHT = 0.3


def load_space(name: str):
    cached = np.load(ROOT / f"data/interim/sceptr/{name}.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors, {t: i for i, t in enumerate(cached["tcr_id"])}


def rank_within(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    return frame.assign(r=frame.groupby("epitope_aa")["v"].rank(pct=True)).set_index("tcr_id")["r"]


def main() -> int:
    test = load_test()
    labels = build_training_frame()
    vectors, position = load_space("default")
    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")

    # ---------- S06 ----------
    print("=== S06: leave-one-clonotype-family-out ===")
    print("   each TCR is scored against a cohort with its own family removed,")
    print("   so no relative can supply its density\n")
    rows = []
    for epitope, cohort in test.groupby("epitope_aa"):
        labelled = labels[labels["epitope_aa"] == epitope]
        if labelled.empty:
            continue
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[position[t] for t in ids]]
        distances = 1.0 - block @ block.T
        _, family = connected_components(
            sp.csr_matrix(distances <= FAMILY_RADIUS), directed=False
        )
        sizes = pd.Series(family).value_counts()

        standard = np.exp(-np.where(np.eye(len(ids), dtype=bool), np.inf, distances) / SIGMA).mean(axis=1)

        # Family-masked: a TCR's own family contributes nothing to its density.
        masked = np.empty(len(ids))
        for i in range(len(ids)):
            keep = family != family[i]
            masked[i] = (
                np.exp(-distances[i, keep] / SIGMA).mean() if keep.any() else 0.0
            )

        y = labelled["pos"].to_numpy()
        order = {t: i for i, t in enumerate(ids)}
        pick = [order[t] for t in labelled["tcr_id"]]
        rows.append(
            {
                "epitope": epitope,
                "n": len(ids),
                "families_size>=2": int((sizes >= 2).sum()),
                "in_a_family": float(np.isin(family, sizes[sizes >= 2].index).mean()),
                "auc_standard": roc_auc_score(y, standard[pick]),
                "auc_family_masked": roc_auc_score(y, masked[pick]),
            }
        )
    s06 = pd.DataFrame(rows).set_index("epitope")
    s06["delta"] = s06["auc_family_masked"] - s06["auc_standard"]
    s06.loc["MEAN"] = s06.mean()
    pd.set_option("display.width", 250)
    print(s06.round(4).to_string())
    print(
        "\n  A large negative delta would mean the ranking is carried by family\n"
        "  recurrence, which the private epitopes cannot supply."
    )

    # ---------- D03 ----------
    print("\n=== D03: alpha/beta chain discordance ===")
    alpha_vectors, alpha_pos = load_space("a_sceptr")
    beta_vectors, beta_pos = load_space("b_sceptr")

    discord = pd.Series(index=test["tcr_id"], dtype=float)
    alpha_only = pd.Series(index=test["tcr_id"], dtype=float)
    beta_only = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].tolist()
        parts = {}
        for name, vecs, pos in (("a", alpha_vectors, alpha_pos), ("b", beta_vectors, beta_pos)):
            block = vecs[[pos[t] for t in ids]]
            distances = 1.0 - block @ block.T
            np.fill_diagonal(distances, np.inf)
            parts[name] = pd.Series(np.exp(-distances / SIGMA).mean(axis=1), index=ids).rank(pct=True)
        alpha_only.loc[ids] = parts["a"].to_numpy()
        beta_only.loc[ids] = parts["b"].to_numpy()
        # Signed disagreement: positive means alpha looks more typical than beta.
        discord.loc[ids] = (parts["a"] - parts["b"]).to_numpy()

    for name, series in (
        ("alpha-only density", alpha_only),
        ("beta-only density", beta_only),
        ("rank(alpha) - rank(beta)", discord),
        ("|rank(alpha) - rank(beta)|", discord.abs()),
    ):
        aucs = per_epitope_auc(labels, labels["tcr_id"].map(series).to_numpy())
        print(f"  {name:<28} mean={aucs.mean():.4f}   "
              + " ".join(f"{e[:5]}={v:.3f}" for e, v in aucs.items()))

    print("\n  within-epitope variance on the private epitopes (E3's rule):")
    private = test[test["split"] == "private"]
    audit = private.assign(discord=private["tcr_id"].map(discord))
    print(audit.groupby("epitope_aa")["discord"].agg(["std", "min", "max"]).round(3).to_string())

    print("\n  does it add to the submitted blend?")
    density = (
        rank_within(
            pd.Series(
                {t: v for t, v in zip(test["tcr_id"], np.zeros(len(test)))}, dtype=float
            ),
            test,
        )
        if False
        else None
    )
    # Rebuild the submitted score components.
    def cohort_density(space: str, sigma: float) -> pd.Series:
        vecs, pos = load_space(space)
        out = pd.Series(index=test["tcr_id"], dtype=float)
        for _, cohort in test.groupby("epitope_aa"):
            block = vecs[[pos[t] for t in cohort["tcr_id"]]]
            distances = 1.0 - block @ block.T
            np.fill_diagonal(distances, np.inf)
            out.loc[cohort["tcr_id"]] = np.exp(-distances / sigma).mean(axis=1)
        return out

    tcrdist = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
        pos = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows_idx = [pos[t] for t in cohort["tcr_id"]]
        distances = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows_idx, rows_idx)]
        np.fill_diagonal(distances, np.inf)
        tcrdist.loc[cohort["tcr_id"]] = np.exp(-distances / 20.0).mean(axis=1)

    blend = (
        rank_within(cohort_density("default", SIGMA), test)
        + rank_within(cohort_density("small", SIGMA), test)
        + rank_within(tcrdist, test)
    ) / 3
    submitted = (1 - STRUCTURE_WEIGHT) * blend + STRUCTURE_WEIGHT * rank_within(
        structure["pae_tcr_pep_min"], test
    )
    base = per_epitope_auc(labels, labels["tcr_id"].map(submitted).to_numpy()).mean()
    for weight in (0.1, 0.2):
        for name, series in (("discord", discord), ("|discord|", discord.abs())):
            combined = (1 - weight) * submitted + weight * rank_within(series, test)
            score = per_epitope_auc(labels, labels["tcr_id"].map(combined).to_numpy()).mean()
            print(f"    + {name:<11} w={weight}  {score:.4f}  gain {score - base:+.4f}")
    print(f"    submitted baseline {base:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
