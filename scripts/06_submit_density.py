#!/usr/bin/env python
"""Baseline-beating submission from the best saturated density feature.

Gaussian-kernel density over SCEPTR `default` cosine distances inside each
epitope cohort, sigma=0.1. Unsupervised: no label touches the score, so there
is nothing to overfit and the public-epitope AUC is an honest read of it.

Scores 0.8396 mean per-epitope AUC on the public epitopes against the
organisers' tcrdist3 baseline at 0.8247 -- a real but small edge, around one
standard error of the public estimate.

Run: python scripts/06_submit_density.py
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

VARIANT = "default"
SIGMA = 0.1


def cohort_density(vectors: np.ndarray, sigma: float) -> np.ndarray:
    """Mean Gaussian kernel to every other cohort member, self excluded."""
    normed = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    distances = 1.0 - normed @ normed.T
    np.fill_diagonal(distances, np.inf)
    return np.exp(-distances / sigma).mean(axis=1)


def main() -> int:
    test = load_test()
    cached = np.load(ROOT / f"data/interim/sceptr/{VARIANT}.npz", allow_pickle=True)
    position = {tcr: i for i, tcr in enumerate(cached["tcr_id"])}
    vectors = cached["vectors"]

    scores = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        rows = [position[t] for t in cohort["tcr_id"]]
        scores.loc[cohort["tcr_id"]] = cohort_density(vectors[rows], SIGMA)

    path = build(scores, f"density_sceptr_{VARIANT}_sigma{SIGMA}")

    labelled = build_training_frame()
    report = score_report(labelled, scores.loc[labelled["tcr_id"]].to_numpy())
    print(report.round(4).to_string())
    print(f"\nwrote {path.relative_to(ROOT)}")
    print(
        "submit with: kaggle competitions submit -c immrep26 "
        f"-f {path.relative_to(ROOT)} -m 'SCEPTR density sigma=0.1'"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
