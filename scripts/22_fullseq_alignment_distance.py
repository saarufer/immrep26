#!/usr/bin/env python
"""Density over whole-V-domain alignment distance, not CDR3 + V/J.

Prompted by a teammate's independent branch, which reports this distance adding
+0.0044 under leave-one-epitope-out and -- unlike candidate-pool expansion --
carrying no private-side discount, because it changes how two TCRs are compared
rather than how many TCRs there are to compare against.

Every distance tried here so far reads the CDR3 loops plus V/J identity:
tcrdist3 by construction, SCEPTR through its tokenised CDR + gene input. The
competition ships `full_seq_alpha_aa` / `full_seq_beta_aa`, the assembled
variable domains, and a straight alignment over those weights the framework
regions too. Whether that is signal or dilution is the question.

Needleman-Wunsch, BLOSUM62, gap open 11 / extend 1, via parasail. Similarity is
normalised so that sequence length cannot dominate:

    d(i,j) = 1 - S(i,j) / sqrt(S(i,i) * S(j,j))

Run: python scripts/22_fullseq_alignment_distance.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

OUTDIR = ROOT / "data/interim/fullseq"
GAP_OPEN, GAP_EXTEND = 11, 1


def cohort_distance(sequences: list[str]) -> np.ndarray:
    """Normalised global-alignment distance matrix for one chain of one cohort."""
    import parasail

    matrix = parasail.blosum62
    n = len(sequences)
    # Self-similarities first: the normaliser.
    self_score = np.array(
        [parasail.nw_scan_16(s, s, GAP_OPEN, GAP_EXTEND, matrix).score for s in sequences],
        dtype=np.float64,
    )
    distances = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            score = parasail.nw_scan_16(
                sequences[i], sequences[j], GAP_OPEN, GAP_EXTEND, matrix
            ).score
            denominator = np.sqrt(self_score[i] * self_score[j])
            value = 1.0 - score / denominator if denominator > 0 else 1.0
            distances[i, j] = distances[j, i] = value
    return distances


def main() -> int:
    test = load_test()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for epitope, cohort in test.groupby("epitope_aa", sort=True):
        out = OUTDIR / f"{epitope}.npz"
        if out.exists():
            print(f"  {epitope:<11} cached")
            continue
        cohort = cohort.reset_index(drop=True)
        alpha = cohort_distance(cohort["full_seq_alpha_aa"].tolist())
        beta = cohort_distance(cohort["full_seq_beta_aa"].tolist())
        np.savez_compressed(
            out, alpha=alpha, beta=beta, tcr_id=cohort["tcr_id"].to_numpy(dtype=object)
        )
        print(
            f"  {epitope:<11} n={len(cohort):4d}  "
            f"alpha[{alpha[alpha > 0].min():.3f},{alpha.max():.3f}] "
            f"beta[{beta[beta > 0].min():.3f},{beta.max():.3f}]"
        )
    print(f"\ncached in {OUTDIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
