#!/usr/bin/env python
"""Compute and cache per-epitope paired-abTCR tcrdist3 distance matrices.

Every density feature downstream is a reduction of these matrices, so they are
computed once and cached as .npz next to a row-order manifest.

Distances are computed *within* each epitope cohort, which is what the
organisers' `nndist` baseline does and the only thing the metric rewards: mean
per-epitope AUC never compares scores across epitopes.

Run: python scripts/03_tcrdist_matrices.py
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

OUTDIR = ROOT / "data/interim/tcrdist"


def to_tcrdist_gene(gene: str | float) -> str | float:
    """test.csv writes `TRAV14_DV4*01`; tcrdist3's reference uses `TRAV14/DV4*01`.

    All 2,636 rows map cleanly under this single substitution (verified: zero
    unmatched levels across TRAV/TRAJ/TRBV/TRBJ).
    """
    return gene.replace("_", "/") if isinstance(gene, str) else gene


def cohort_distances(cohort: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Paired alpha and beta distance matrices for one epitope cohort."""
    from tcrdist.repertoire import TCRrep

    frame = pd.DataFrame(
        {
            "cdr3_a_aa": cohort["cdr3_alpha_aa"].to_numpy(),
            "v_a_gene": cohort["TRAV_IMGT"].map(to_tcrdist_gene).to_numpy(),
            "j_a_gene": cohort["TRAJ_IMGT"].map(to_tcrdist_gene).to_numpy(),
            "cdr3_b_aa": cohort["cdr3_beta_aa"].to_numpy(),
            "v_b_gene": cohort["TRBV_IMGT"].map(to_tcrdist_gene).to_numpy(),
            "j_b_gene": cohort["TRBJ_IMGT"].map(to_tcrdist_gene).to_numpy(),
            "count": 1,
        }
    )
    # deduplicate=False keeps a 1:1 row correspondence with the cohort, so the
    # matrix index is the cohort index. Identical clonotypes would otherwise be
    # collapsed and the row order lost.
    rep = TCRrep(
        cell_df=frame,
        organism="human",
        chains=["alpha", "beta"],
        compute_distances=True,
        deduplicate=False,
    )
    return rep.pw_alpha.astype(np.float32), rep.pw_beta.astype(np.float32)


def main() -> int:
    test = load_test()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for epitope, cohort in test.groupby("epitope_aa", sort=True):
        cohort = cohort.reset_index(drop=True)
        out = OUTDIR / f"{epitope}.npz"
        if out.exists():
            print(f"  {epitope:<11} cached")
            continue
        alpha, beta = cohort_distances(cohort)
        if alpha.shape != (len(cohort), len(cohort)):
            raise ValueError(f"{epitope}: matrix shape {alpha.shape} != {len(cohort)}")
        np.savez_compressed(
            out,
            alpha=alpha,
            beta=beta,
            tcr_id=cohort["tcr_id"].to_numpy(dtype=object),
        )
        print(
            f"  {epitope:<11} n={len(cohort):4d} "
            f"alpha[{alpha.min():.0f},{alpha.max():.0f}] "
            f"beta[{beta.min():.0f},{beta.max():.0f}] -> {out.name}"
        )
    print(f"\ncached matrices in {OUTDIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
