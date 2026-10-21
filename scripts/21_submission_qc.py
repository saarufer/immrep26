#!/usr/bin/env python
"""Assertions that a matching public score cannot provide.

A submission predicted locally at 0.8498 scoring 0.849 confirms that the metric,
the label join and the public rows are all right. It confirms nothing about the
602 private rows: scoring covers only the 1,782 labelled public rows, so a
misalignment confined to the private epitopes would be invisible in that number.
Claiming otherwise was an error in an earlier report; these checks replace the
claim.

Beyond id-set equality, the decisive check is provenance: for a sample weighted
towards the private rows, confirm that each complex's input YAML holds that
tcr_id's *own* chain sequences and the HLA platform matching its own allele.
That catches the failure that id-set equality cannot -- right identifier, wrong
structure behind it.

Run: python scripts/21_submission_qc.py [submission.csv]
"""

from __future__ import annotations

import collections
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import N_TEST_ROWS, PRIVATE_EPITOPES  # noqa: E402
from immrep.data import load_test  # noqa: E402

DEFAULT_SUBMISSION = ROOT / "submissions/density3_struct03.csv"
INPUTS = ROOT / "data/interim/boltz/inputs"
SAMPLE_PER_SPLIT = 40
SEED = 0


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SUBMISSION
    test = load_test()
    submission = pd.read_csv(path)
    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv")

    failures: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}{'  ' + detail if detail else ''}")
        if not ok:
            failures.append(name)

    print(f"submission: {path.relative_to(ROOT)}\n")
    print("=== submission integrity ===")
    check("tcr_id unique", bool(submission["tcr_id"].is_unique))
    check(
        "row order identical to test.csv",
        submission["tcr_id"].tolist() == test["tcr_id"].tolist(),
    )
    check("every prediction finite", bool(np.isfinite(submission["prediction"]).all()))
    check(f"{N_TEST_ROWS} rows", len(submission) == N_TEST_ROWS, f"got {len(submission)}")

    print("\n=== features joined by id, not by row position ===")
    check(
        "structure covers every test tcr_id",
        set(test["tcr_id"]) == set(structure["tcr_id"]),
        f"{len(set(test['tcr_id']) - set(structure['tcr_id']))} missing",
    )
    check("structure tcr_id unique", bool(structure["tcr_id"].is_unique))

    print("\n=== provenance: does each complex hold its own tcr_id's sequences? ===")
    rng = np.random.default_rng(SEED)
    sample = np.concatenate(
        [
            rng.choice(
                test.index[test["split"] == split].to_numpy(),
                min(SAMPLE_PER_SPLIT, int((test["split"] == split).sum())),
                replace=False,
            )
            for split in ("private", "public")
        ]
    )
    mismatches = []
    platforms: dict[str, set[str]] = collections.defaultdict(set)
    for i in sample:
        row = test.loc[i]
        sequences = re.findall(
            r"sequence: (\S+)", (INPUTS / f"{row.tcr_id}.yaml").read_text()
        )
        if len(sequences) != 4:
            mismatches.append((row.tcr_id, f"{len(sequences)} chains"))
            continue
        if sequences[0] != row.full_seq_alpha_aa:
            mismatches.append((row.tcr_id, "alpha chain"))
        elif sequences[1] != row.full_seq_beta_aa:
            mismatches.append((row.tcr_id, "beta chain"))
        elif sequences[3] != row.epitope_aa:
            mismatches.append((row.tcr_id, f"peptide {sequences[3]} != {row.epitope_aa}"))
        platforms[row.hla_long].add(sequences[2])

    check(
        f"{len(sample)} sampled complexes carry their own sequences "
        f"({int((test.loc[sample, 'split'] == 'private').sum())} private)",
        not mismatches,
        str(mismatches[:3]),
    )
    check(
        "one MHC platform sequence per allele",
        all(len(v) == 1 for v in platforms.values()),
        str({k: len(v) for k, v in platforms.items()}),
    )
    check(
        "distinct alleles use distinct platforms",
        len({next(iter(v)) for v in platforms.values()}) == len(platforms),
    )

    print("\n=== private rows specifically ===")
    private = test[test["epitope_aa"].isin(PRIVATE_EPITOPES)].merge(
        structure, on="tcr_id", how="left"
    )
    check(
        "no private row missing structure features",
        bool(private["pae_tcr_pep_min"].notna().all()),
        f"{int(private['pae_tcr_pep_min'].isna().sum())} missing",
    )
    joined = test[["tcr_id", "epitope_aa"]].merge(submission, on="tcr_id")
    flat = [
        epitope
        for epitope, group in joined[joined["epitope_aa"].isin(PRIVATE_EPITOPES)].groupby(
            "epitope_aa"
        )
        if group["prediction"].nunique() < 10
    ]
    check(
        "every private epitope has a spread of scores",
        not flat,
        f"near-constant: {flat}" if flat else "",
    )

    print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILURES: ' + str(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
