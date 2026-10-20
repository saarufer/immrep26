#!/usr/bin/env python
"""Generation probability (OLGA Pgen) for every test TCR, both chains.

Rationale: density features only rank a TCR relative to its neighbours, so they
say little about the 34-83% of each cohort that sits in no cluster. Pgen is
neighbour-free. A spuriously annotated TCR -- e.g. background carried through a
dextramer sort -- should look like a typical draw from the naive repertoire,
whereas an antigen-selected clonotype is expected to be rarer than chance. Pgen
also separates "convergent because antigen-driven" from "convergent because
easy to generate", which plain density conflates.

Run: python scripts/04_pgen.py
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

OUT = ROOT / "data/interim/pgen.csv"


def compute_pgen(sequences: pd.Series, chain_folder: str, recomb_type: str) -> pd.Series:
    """Pgen per unique CDR3, then broadcast back. Uniques save ~30% of calls."""
    from tcrdist.pgen import OlgaModel

    model = OlgaModel(chain_folder=chain_folder, recomb_type=recomb_type)
    uniques = sequences.dropna().unique().tolist()
    values = model.compute_aa_cdr3_pgens(uniques)
    lookup = dict(zip(uniques, values))
    return sequences.map(lookup)


def main() -> int:
    test = load_test()
    OUT.parent.mkdir(parents=True, exist_ok=True)

    frame = test[["tcr_id", "epitope_aa", "cdr3_alpha_aa", "cdr3_beta_aa"]].copy()
    frame["pgen_alpha"] = compute_pgen(frame["cdr3_alpha_aa"], "human_T_alpha", "VJ")
    frame["pgen_beta"] = compute_pgen(frame["cdr3_beta_aa"], "human_T_beta", "VDJ")

    # Pgen spans many orders of magnitude, so work in log space. Zeros occur
    # when OLGA cannot generate the sequence under its model; floor them one
    # decade below the smallest observed non-zero value rather than dropping
    # the row, since "ungeneratable" is itself informative.
    for chain in ("alpha", "beta"):
        col = f"pgen_{chain}"
        positive = frame.loc[frame[col] > 0, col]
        floor = positive.min() / 10 if len(positive) else 1e-30
        frame[f"log10_{col}"] = np.log10(frame[col].clip(lower=floor))
        frame[f"{col}_is_zero"] = (frame[col] <= 0).astype(int)

    # A paired TCR's joint probability under independent chain recombination.
    frame["log10_pgen_paired"] = frame["log10_pgen_alpha"] + frame["log10_pgen_beta"]

    frame.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}  rows={len(frame)}")
    print(
        frame[["log10_pgen_alpha", "log10_pgen_beta", "log10_pgen_paired"]]
        .describe()
        .round(2)
        .to_string()
    )
    print(
        "\nungeneratable (Pgen=0): "
        f"alpha={int(frame.pgen_alpha_is_zero.sum())} "
        f"beta={int(frame.pgen_beta_is_zero.sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
