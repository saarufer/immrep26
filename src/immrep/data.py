"""Loading and joining the IMMREP26 competition files with the TCRvdb labels."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .constants import N_TEST_ROWS, PRIVATE_EPITOPES, PUBLIC_EPITOPES

ROOT = Path(__file__).resolve().parents[2]
TEST_CSV = ROOT / "data/raw/test.csv"
TCRDIST_CSV = ROOT / "data/raw/tcrdist_predictions.csv"
TCRVDB_CSV = ROOT / "data/external/tcrvdb/ttd_gil_ylq_nlv_glc_01-10-2026_TCRvdb.csv"


def load_test() -> pd.DataFrame:
    """All 2,636 rows to be scored, with a `split` column added.

    Row order is preserved because submissions must keep it.
    """
    test = pd.read_csv(TEST_CSV)
    if len(test) != N_TEST_ROWS:
        raise ValueError(f"expected {N_TEST_ROWS} rows in test.csv, got {len(test)}")
    test["split"] = test["epitope_aa"].map(
        lambda e: "public" if e in PUBLIC_EPITOPES else "private"
    )
    unknown = set(test["epitope_aa"]) - set(PUBLIC_EPITOPES) - set(PRIVATE_EPITOPES)
    if unknown:
        raise ValueError(f"unexpected epitopes in test.csv: {sorted(unknown)}")
    return test


def load_tcrdist_baseline() -> pd.DataFrame:
    """The organisers' baseline: `prediction` is the negated within-epitope
    nearest-neighbour tcrdist3 distance (`nndist`)."""
    return pd.read_csv(TCRDIST_CSV)


def load_labels(*, drop_ambiguous: bool = True) -> pd.DataFrame:
    """TCRvdb functional-validation labels for the five public epitopes.

    Returns `pos` as 0/1. The 252 ambiguous rows carry no usable label and are
    excluded from the official metric, so they are dropped by default.

    TCRvdb is request-only and its terms forbid redistribution, so this file is
    gitignored; see data/external/README.md.
    """
    if not TCRVDB_CSV.exists():
        raise FileNotFoundError(
            f"{TCRVDB_CSV} missing. Request TCRvdb at https://tcrvdb.nkiai.nl"
        )
    labels = pd.read_csv(TCRVDB_CSV)
    one_hot = labels[["pos", "neg", "ambiguous"]].sum(axis=1)
    if not (one_hot == 1).all():
        raise ValueError("pos/neg/ambiguous are not mutually exclusive in TCRvdb")
    if drop_ambiguous:
        labels = labels[labels["ambiguous"] == 0].copy()
    labels["pos"] = labels["pos"].astype(int)
    return labels


def build_training_frame() -> pd.DataFrame:
    """Public-epitope rows of test.csv joined to their TCRvdb labels.

    This is the whole supervised training set: 1,782 labelled rows that live
    inside the test file itself. Joined on `tcr_id`, which TCRvdb shares with
    test.csv verbatim.
    """
    test = load_test()
    labels = load_labels()
    merged = test.merge(labels[["tcr_id", "pos"]], on="tcr_id", how="inner")
    expected = len(labels)
    if len(merged) != expected:
        raise ValueError(
            f"label join lost rows: {len(merged)} joined vs {expected} labels"
        )
    return merged
