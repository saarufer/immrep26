#!/usr/bin/env python
"""Build the final submission under the fixed-configuration discipline.

Configuration, fixed and not re-selected:

  rank-average of three within-epitope kernel densities
      SCEPTR `default`  (cosine, sigma=0.1)
      SCEPTR `small`    (cosine, sigma=0.1)
      tcrdist3 alpha+beta (sigma=20)
  blended 0.8 / 0.2 with the TCR-peptide min-PAE from Boltz-2.

Why these choices and no others, from reports/00_experiment_log.md:

- Averaging representations rather than picking one: every density lands within
  0.01 of every other, so the choice is noise and averaging removes it.
- Structure at a low weight: it reaches only 0.757 alone and correlates +0.53
  with density, so it is a minor correction, not a peer. Nested LOEO selects
  0.2-0.3 in every fold.
- No fitted weights anywhere. Logistic regression on top of density loses 0.006
  and gradient boosting 0.021 (E10), and cross-validated configuration
  selection loses 0.0048 against fixing it in advance (E15).
- Nothing else enters: metadata, generation probability, background
  normalisation, cluster shape, graph refinement, an unseen-epitope specificity
  model and multi-sample docking consistency were each measured and each failed
  the +0.02 inclusion threshold, most of them scoring below chance.

Expected private score: **0.846**, the nested-LOEO estimate, which already
absorbs selection cost. Public-epitope score will read about 0.851 but that
number is not evidence of anything -- the public labels are published.

Run: python scripts/18_final_submission.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.metric import score_report  # noqa: E402
from immrep.submission import build  # noqa: E402

SCEPTR_SPACES = (("default", 0.1), ("small", 0.1))
TCRDIST_SIGMA = 20.0
STRUCTURE_WEIGHT = 0.2
NAME = "density3_struct02"


def rank_within_epitope(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    frame["r"] = frame.groupby("epitope_aa")["v"].rank(pct=True)
    return frame.set_index("tcr_id")["r"]


def sceptr_density(test: pd.DataFrame, space: str, sigma: float) -> pd.Series:
    cached = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    vectors = cached["vectors"]
    values = pd.Series(index=test["tcr_id"], dtype=float)
    for _, cohort in test.groupby("epitope_aa"):
        block = vectors[[position[t] for t in cohort["tcr_id"]]]
        block = block / np.linalg.norm(block, axis=1, keepdims=True)
        distances = 1.0 - block @ block.T
        np.fill_diagonal(distances, np.inf)
        values.loc[cohort["tcr_id"]] = np.exp(-distances / sigma).mean(axis=1)
    return values


def tcrdist_density(test: pd.DataFrame, sigma: float) -> pd.Series:
    values = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [position[t] for t in cohort["tcr_id"]]
        distances = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
        np.fill_diagonal(distances, np.inf)
        values.loc[cohort["tcr_id"]] = np.exp(-distances / sigma).mean(axis=1)
    return values


def main() -> int:
    test = load_test()

    components = [
        rank_within_epitope(sceptr_density(test, space, sigma), test)
        for space, sigma in SCEPTR_SPACES
    ]
    components.append(rank_within_epitope(tcrdist_density(test, TCRDIST_SIGMA), test))
    density = sum(components) / len(components)

    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure_rank = rank_within_epitope(structure["pae_tcr_pep_min"], test)

    final = (1 - STRUCTURE_WEIGHT) * density + STRUCTURE_WEIGHT * structure_rank

    path = build(final, NAME)

    labelled = build_training_frame()
    report = score_report(labelled, labelled["tcr_id"].map(final).to_numpy())
    print(report.round(4).to_string())
    print("\nnested-LOEO estimate of private performance: 0.846")
    print("organisers' baselines: tcrdist3 0.8247, AlphaFold3 0.787")
    print(f"\nwrote {path.relative_to(ROOT)}")

    # Sanity: the private epitopes must each have a spread of scores, or that
    # epitope's AUC silently becomes 0.5 -- a third of the private score.
    joined = test[["tcr_id", "epitope_aa"]].copy()
    joined["prediction"] = joined["tcr_id"].map(final)
    private = joined[joined["epitope_aa"].isin(("CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"))]
    print("\nprivate-epitope score spread:")
    print(
        private.groupby("epitope_aa")["prediction"]
        .agg(["size", "nunique", "min", "max"])
        .round(4)
        .to_string()
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
