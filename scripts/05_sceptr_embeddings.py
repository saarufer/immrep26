#!/usr/bin/env python
"""SCEPTR representations for every test TCR, and cohort distance matrices.

tcrdist3 density saturates around 0.83 mean per-epitope AUC, and generation
probability adds nothing orthogonal (scripts/04_pgen.py). The remaining lever on
a neighbourhood-based signal is a better metric space, so the same density
reductions are recomputed over a learned paired-chain representation.

Variants are kept separate rather than averaged: in tcrdist space the alpha
chain alone outranks the beta chain alone, so the single-chain models
(`a_sceptr`, `b_sceptr`) are worth measuring on their own.

Run: python scripts/05_sceptr_embeddings.py
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
from immrep.constants import ALL_EPITOPES  # noqa: E402

OUTDIR = ROOT / "data/interim/sceptr"

# `default` is the published paired-chain model, `large` its bigger sibling.
# The single-chain variants isolate each chain's contribution, which matters
# here because in tcrdist space the alpha chain alone outranks the beta chain
# alone. `cdr3_only` drops the V/J context, separating motif signal from gene
# usage; `finetuned` was tuned on labelled specificity data.
VARIANTS = ("default", "large", "small", "a_sceptr", "b_sceptr", "cdr3_only", "finetuned")


def to_sceptr_frame(test: pd.DataFrame) -> pd.DataFrame:
    """SCEPTR's prescribed columns. Gene names use the same slash form as
    tcrdist3, so the `_` -> `/` substitution from script 03 applies.

    SCEPTR's embedding table only covers functional V genes. Exactly one row in
    test.csv falls outside it (TRBV23-1, an ORF, in the private CINGVCWTV
    cohort). Dropping it is not an option because every row needs a prediction,
    so its V gene is set to None, which SCEPTR accepts and handles by embedding
    that chain from its CDR3 alone.
    """
    from libtcrlm.schema.tcr import TravGene, TrbvGene

    slash = lambda g: g.replace("_", "/") if isinstance(g, str) else g  # noqa: E731
    frame = pd.DataFrame(
        {
            "TRAV": test["TRAV_IMGT"].map(slash),
            "CDR3A": test["cdr3_alpha_aa"],
            "TRAJ": test["TRAJ_IMGT"].map(slash),
            "TRBV": test["TRBV_IMGT"].map(slash),
            "CDR3B": test["cdr3_beta_aa"],
            "TRBJ": test["TRBJ_IMGT"].map(slash),
        }
    )
    for col, vocab in (
        ("TRAV", {g.name for g in TravGene}),
        ("TRBV", {g.name for g in TrbvGene}),
    ):
        outside = ~frame[col].str.split("*").str[0].isin(vocab)
        if outside.any():
            print(
                f"  note: {col} outside SCEPTR vocabulary for "
                f"{int(outside.sum())} row(s) "
                f"({sorted(frame.loc[outside, col].unique())}); "
                "falling back to CDR3-only for that chain"
            )
            frame.loc[outside, col] = None
    return frame


def main() -> int:
    import sceptr.variant as variants

    test = load_test()
    frame = to_sceptr_frame(test)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for name in VARIANTS:
        out = OUTDIR / f"{name}.npz"
        if out.exists():
            print(f"  {name:<10} cached")
            continue
        model = getattr(variants, name)()
        vectors = model.calc_vector_representations(frame)
        if vectors.shape[0] != len(test):
            raise ValueError(f"{name}: got {vectors.shape[0]} vectors for {len(test)} TCRs")

        # Cosine distance within each epitope cohort. Cosine rather than
        # Euclidean because SCEPTR's published distance is cosine-based.
        normed = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
        cohorts: dict[str, np.ndarray] = {}
        ids: dict[str, np.ndarray] = {}
        for epitope in sorted(ALL_EPITOPES):
            mask = (test["epitope_aa"] == epitope).to_numpy()
            block = normed[mask]
            cohorts[epitope] = (1.0 - block @ block.T).astype(np.float32)
            ids[epitope] = test.loc[mask, "tcr_id"].to_numpy(dtype=object)

        np.savez_compressed(
            out,
            vectors=vectors.astype(np.float32),
            tcr_id=test["tcr_id"].to_numpy(dtype=object),
            **{f"D__{e}": m for e, m in cohorts.items()},
            **{f"ID__{e}": v for e, v in ids.items()},
        )
        print(f"  {name:<10} dim={vectors.shape[1]:4d} -> {out.name}")

    print(f"\ncached in {OUTDIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
