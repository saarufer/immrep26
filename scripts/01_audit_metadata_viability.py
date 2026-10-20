#!/usr/bin/env python
"""Audit which VDJdb metadata fields can possibly rank TCRs *within* each epitope.

Motivation
----------
IMMREP26 is scored per split: the five published epitopes land on the public
leaderboard, the three held-back epitopes decide the final ranking. Any feature
is therefore only useful if it varies *within* a private epitope. A field that
is constant across all 243 CINGVCWTV rows contributes exactly zero ranking
information there, no matter how predictive it looks in cross-validation on the
public epitopes.

This script quantifies that per epitope, so the feature set can be restricted to
fields that survive on the private side. It writes a markdown table to
reports/ and prints the verdict.

Run: python scripts/01_audit_metadata_viability.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import PRIVATE_EPITOPES, PUBLIC_EPITOPES, VDJDB_VERSION  # noqa: E402

VDJDB_FULL = ROOT / "data/external/vdjdb" / f"vdjdb-{VDJDB_VERSION}" / "vdjdb_full.txt"
OUT = ROOT / "reports/01_metadata_viability.md"

# Fields a "predict annotation reliability" model would reach for first.
CANDIDATE_FIELDS = [
    "vdjdb.score",
    "reference.id",
    "method.identification",
    "method.verification",
    "method.singlecell",
    "method.sequencing",
    "method.frequency",
    "meta.cell.subset",
    "meta.subject.cohort",
    "meta.subject.id",
    "meta.tissue",
    "meta.donor.MHC",
    "mhc.a",
]

# A field with two levels split 132:3 is formally "variable" but carries almost
# no usable signal, so we also report the share held by the dominant level.
SKEW_CEILING = 0.95


def load_paired_vdjdb(path: Path) -> pd.DataFrame:
    """Human entries with both chains resolved, one row per paired clonotype."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    paired = df[
        (df["species"] == "HomoSapiens")
        & df["cdr3.alpha"].notna()
        & df["cdr3.beta"].notna()
    ]
    return paired.drop_duplicates(subset=["cdr3.alpha", "cdr3.beta"])


def audit(paired: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    private, public = list(PRIVATE_EPITOPES), list(PUBLIC_EPITOPES)
    cohorts = {ep: paired[paired["antigen.epitope"] == ep] for ep in private + public}

    levels, skew = [], []
    for field in CANDIDATE_FIELDS:
        lv = {"feature": field}
        sk = {"feature": field}
        for ep, cohort in cohorts.items():
            values = cohort[field].dropna()
            lv[ep] = values.nunique()
            sk[ep] = (
                values.value_counts(normalize=True).iloc[0] if len(values) else float("nan")
            )
        levels.append(lv)
        skew.append(sk)

    lv = pd.DataFrame(levels).set_index("feature")
    sk = pd.DataFrame(skew).set_index("feature")
    # Usable on the private side only if it varies, and not 132:3-style, in all three.
    lv["private_min_levels"] = lv[private].min(axis=1)
    lv["private_max_skew"] = sk[private].max(axis=1)
    lv["usable_private"] = (lv["private_min_levels"] >= 2) & (
        lv["private_max_skew"] < SKEW_CEILING
    )
    return lv, sk


def main() -> int:
    if not VDJDB_FULL.exists():
        sys.exit(f"missing {VDJDB_FULL}; run scripts/00_fetch_external.sh first")

    paired = load_paired_vdjdb(VDJDB_FULL)
    levels, skew = audit(paired)

    private = list(PRIVATE_EPITOPES)
    usable = sorted(levels.index[levels["usable_private"]])
    dead = sorted(levels.index[~levels["usable_private"]])

    cols = private + list(PUBLIC_EPITOPES)
    lines = [
        "# Metadata viability inside private epitopes",
        "",
        f"VDJdb `{VDJDB_VERSION}`, {len(paired):,} paired human clonotypes.",
        "Cells are the number of distinct non-null levels within that epitope's cohort.",
        f"A field counts as usable only if it has >=2 levels in all three private "
        f"epitopes and its dominant level holds <{SKEW_CEILING:.0%} of rows.",
        "",
        "| feature | " + " | ".join(cols) + " | max skew (private) | usable |",
        "|---|" + "---|" * (len(cols) + 2),
    ]
    for field, row in levels.iterrows():
        cells = " | ".join(str(int(row[c])) for c in cols)
        lines.append(
            f"| `{field}` | {cells} | {row['private_max_skew']:.2f} | "
            f"{'yes' if row['usable_private'] else '**no**'} |"
        )
    lines += [
        "",
        f"**Usable on all three private epitopes:** {', '.join(f'`{u}`' for u in usable) or 'none'}",
        "",
        f"**Dead on at least one private epitope:** {', '.join(f'`{d}`' for d in dead)}",
        "",
        "## Consequence",
        "",
        "Annotation-provenance features cannot rank the private epitopes, so "
        "leave-one-epitope-out scores on the public epitopes overstate what they "
        "will deliver on the final ranking. Treat them as an ablation / diagnostic, "
        "not as a submission component, and keep the ranking model on paired-abTCR "
        "sequence and structure, which transfer.",
        "",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")

    print(levels[cols + ["private_min_levels", "private_max_skew", "usable_private"]].to_string())
    print(f"\nusable on all 3 private epitopes: {usable or 'NONE'}")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
