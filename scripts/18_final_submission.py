#!/usr/bin/env python
"""Build the final submission under the fixed-configuration discipline.

Configuration, fixed and not re-selected:

  rank-average of three within-epitope kernel densities
      SCEPTR `default`  (cosine, sigma=0.1)
      SCEPTR `small`    (cosine, sigma=0.1)
      tcrdist3 alpha+beta (sigma=20)
  blended 0.7 / 0.3 with the TCR-peptide min-PAE from Boltz-2.

Why these choices and no others, from reports/00_experiment_log.md:

- Averaging representations rather than picking one: every density lands within
  0.01 of every other, so the choice is noise and averaging removes it.
- Structure at a low weight: it reaches only 0.757 alone and correlates +0.53
  with density, so it is a minor correction, not a peer. Nested LOEO selects
  0.2-0.3 in every fold.
- Weight 0.3 rather than 0.2 on a dominance argument, not by optimising any one
  regime (E21). Thinning the public cohorts to private-like sparsity shows
  structure's contribution concentrates exactly where density fails: on the two
  epitopes that collapse under thinning its gain triples, from +0.006 to +0.020,
  while on the robust three it stays flat. Sweeping the weight, 0.3 ties 0.2 at
  full coverage (-0.0002) and beats it when sparse (+0.0042 at 45% coverage,
  +0.0008 at 20%), so it is weakly dominant and the private cohorts are the
  sparse case. Weight 0.5 is better still at 20% coverage but loses 0.013 at
  full coverage, so its risk is not symmetric. Note what that -0.0002 is and is
  not: it is the observed difference on the five tested epitopes, not a bound on
  what 0.3 could cost on an unseen one.
- No fitted weights anywhere. Logistic regression on top of density loses 0.006
  and gradient boosting 0.021 (E10), and cross-validated configuration
  selection loses 0.0048 against fixing it in advance (E15).

  To be precise about what that buys: the scoring function has no fitted
  parameters, but its architecture, features and blend weight were chosen in
  light of labelled validation data. Having no fitted parameters is not the same
  as having no selection exposure -- E15 measured that exposure at 0.0048 over
  323 candidates.
- Nothing else enters: metadata, generation probability, background
  normalisation, cluster shape, graph refinement, an unseen-epitope specificity
  model and multi-sample docking consistency were each measured and each failed
  the +0.02 inclusion threshold, most of them scoring below chance.

Expected private score: **0.80 central, 0.67-0.91 plausible** (E30), not the
0.846 nested-LOEO figure. LOEO holds out epitope identity but not cohort
structure, and the private epitopes have four times less clonotype-family
structure than the public five (7.2% against 28.8% of TCRs in a family). They
therefore lose most of the easy clustered-positive-against-isolated-negative
pairs that are 38.6% of public comparisons at a 2% error rate, and are scored
61-93% on both-singleton pairs instead. Reprojecting each class's measured error
rate onto the private compositions gives 0.799, and the spread comes from one
unknowable quantity: how separable those isolated TCRs are. Public-epitope score will read about 0.850 but that
number is not evidence of anything -- the public labels are published.

A matching public score also does not verify that the private rows are
correctly aligned: scoring covers only the 1,782 labelled public rows, so a
misalignment confined to the 602 private rows would be invisible. That is
checked separately, by asserting id uniqueness and order against test.csv, that
every tcr_id carries its own structure row, and that each complex's input YAML
contains that tcr_id's own chain sequences and the HLA platform matching its
allele.

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
STRUCTURE_WEIGHT = 0.3
NAME = "density3_struct03"


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
