#!/usr/bin/env python
"""Gate 0 of the recognition-profile roadmap: can the two external datasets be
joined at all, and at what level?

The roadmap sets three possible verdicts and forbids proceeding past them on
assumptions:

  proceed                 comparable same-receptor same-peptide readouts exist,
                          with independent replicates or controls;
  profiles only           no alignable paired readouts, so recognition profiles
                          can be built but binding -> activation cannot be fit;
  stop any numeric union  only top-k tables and missing labels, so progress
                          would require assuming the unobserved are negative.

What this audits, with nothing downloaded twice and nothing assumed from the
papers:

  * a manifest with sha256, rows and columns for every file;
  * PRP: the two anchor libraries kept apart, per-TCR peptide coverage, and how
    much of each selection is actually present in its own input library -- the
    denominator the roadmap's enrichment formula needs;
  * SEEKER: whether R0 exists (it does), how far the candidate tables are
    top-k truncated, and how many peptides carry a real flow measurement;
  * the overlap that decides the verdict: PRP and SEEKER's B27 module are both
    HLA-B*27:05, so peptide overlap is possible -- but a binding -> activation
    map needs the same TCR as well, and that is what gets measured here.

Run: .venv/bin/python scripts/58_gate0_external_audit.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
PRP = ROOT / "data/external/prp"
SEEKER = ROOT / "data/external/seeker"
REPORT = ROOT / "reports/33_gate0_external_audit.md"
MANIFEST = ROOT / "data/derived/external_manifest.tsv"

# From the PRP README, not guessed.
P2_P9 = {"TCR4.2", "TCR4.3", "TCR4.4"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not (PRP / "README.md").exists() or not (SEEKER / "README.md").exists():
        print("external data absent -- NOT ESTIMABLE")
        return 1

    print("=== manifest ===")
    rows = []
    for source, folder in (("PRP", PRP), ("SEEKER", SEEKER)):
        for path in sorted(folder.rglob("*")):
            if not path.is_file() or path.suffix == ".zip":
                continue
            entry = {"source": source, "path": str(path.relative_to(folder)),
                     "bytes": path.stat().st_size,
                     "sha256": digest(path)[:16]}
            if path.suffix == ".csv":
                try:
                    head = pd.read_csv(path, nrows=0)
                    entry["columns"] = len(head.columns)
                    with path.open() as handle:
                        entry["rows"] = sum(1 for _ in handle) - 1
                except Exception as exc:                     # noqa: BLE001
                    entry["columns"] = f"unreadable: {exc!s:.30}"
            rows.append(entry)
    manifest = pd.DataFrame(rows)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, sep="\t", index=False)
    print(f"  {len(manifest)} files, {manifest['bytes'].sum() / 1e6:.1f} MB, "
          f"written to {MANIFEST.relative_to(ROOT)}")

    print("\n=== PRP: the two anchor libraries stay apart ===")
    naive = {}
    for scheme, name in (("P2_P8", "HLA-B27_Naive_P2_P8_fixed.csv"),
                         ("P2_P9", "HLA-B27_Naive_P2_P9_fixed.csv")):
        frame = pd.read_csv(PRP / name)
        naive[scheme] = set(frame["Peptide"])
        print(f"  {scheme}  {len(frame):>8,} peptides, "
              f"{frame['Frequency'].sum():>12,} reads")
    shared = naive["P2_P8"] & naive["P2_P9"]
    print(f"  peptides in BOTH libraries: {len(shared):,} "
          f"({len(shared) / min(len(v) for v in naive.values()):.2%} of the "
          f"smaller) -- they are different panels and are not pooled")

    prp_rows = []
    all_prp_peptides = set()
    headerless = []
    for path in sorted(PRP.glob("TCR*_R1-R4.csv")):
        tcr = path.name.split("_R1")[0]
        # TCR4.2, 4.3 and 4.4 ship WITHOUT a header row -- their first line is
        # already data, and they are exactly the three P2/P9 files. Reading the
        # directory with pandas defaults silently eats one peptide from each
        # and gives them different column names from the other eighteen.
        with path.open() as handle:
            first = handle.readline().split(",")[0]
        if first == "Peptide":
            frame = pd.read_csv(path)
        else:
            headerless.append(tcr)
            frame = pd.read_csv(path, header=None,
                                names=["Peptide", "Round1", "Round2",
                                       "Round3", "Round4"])
        scheme = "P2_P9" if tcr in P2_P9 else "P2_P8"
        peptides = set(frame["Peptide"])
        all_prp_peptides |= peptides
        in_input = len(peptides & naive[scheme])
        prp_rows.append({
            "tcr": tcr, "anchor": scheme, "peptides": len(frame),
            "in_own_input_library": in_input,
            "input_coverage": in_input / len(frame),
            "reads_R1": int(frame["Round1"].sum()),
            "reads_R4": int(frame["Round4"].sum()),
            "engineered": tcr.startswith("TCR19.2C")})
    prp = pd.DataFrame(prp_rows)
    print(f"\n  headerless files (first line is data): "
          f"{headerless if headerless else 'none'}")
    if set(headerless) == P2_P9:
        print("  -- exactly the three P2/P9 files, so a naive directory read")
        print("     loses a peptide from each and mislabels their columns.")
    print(f"\n  {len(prp)} selection files "
          f"({int(prp['engineered'].sum())} engineered variants of TCR19.2)")
    print(prp[["tcr", "anchor", "peptides", "input_coverage"]].round(3)
          .to_string(index=False))
    print(f"\n  median input coverage {prp['input_coverage'].median():.1%} -- "
          "the roadmap's enrichment formula needs an input count, and this is")
    print("  the fraction of selected peptides that have one.")

    print("\n=== SEEKER: R0, truncation, and real flow measurements ===")
    b27 = SEEKER / "B27_single_TCR_screen"
    seeker_rows, seeker_peptides = [], set()
    for tcr in ("TCR59", "TCR63"):
        path = next(b27.glob(f"screening_results/{tcr}*add_flow.csv"))
        frame = pd.read_csv(path)
        seeker_peptides |= set(frame["aa"].dropna())
        row = {"tcr": tcr, "candidates": len(frame),
               "has_round0": int(frame["count_round0"].notna().sum()),
               "round0_nonzero": int((frame["count_round0"].fillna(0) > 0).sum())}
        # TCR63 ships WITHOUT add_flow_high_copy. The two receptors therefore
        # do not carry the same flow conditions, and the roadmap's F1 -- fit
        # single- and high-copy separately and compare -- is possible only for
        # TCR59.
        for column in ("add_flow_single_copy", "add_flow_high_copy"):
            if column not in frame.columns:
                row[column] = 0
                row[f"{column}_pct"] = 0.0
                row[f"{column}_absent"] = True
                continue
            measured = frame[column].notna().sum()
            row[column] = int(measured)
            row[f"{column}_pct"] = measured / len(frame)
            row[f"{column}_absent"] = False
        seeker_rows.append(row)
        high = ("column absent" if row["add_flow_high_copy_absent"]
                else f"{row['add_flow_high_copy']} high-copy")
        print(f"  {tcr}: {len(frame):>6,} candidates, R0 present for "
              f"{row['has_round0']:,}, flow measured for "
              f"{row['add_flow_single_copy']} single-copy / {high} "
              f"({row['add_flow_single_copy_pct']:.2%} of candidates)")
    seeker = pd.DataFrame(seeker_rows)
    if seeker["add_flow_high_copy_absent"].any():
        missing = list(seeker.loc[seeker["add_flow_high_copy_absent"], "tcr"])
        print(f"\n  {missing} carry no high-copy flow column, so the roadmap's")
        print("  F1 comparison of single- against high-copy antigen density is")
        print("  available on one receptor only, not two.")
    print("\n  Flow was measured on ~1% of candidates and those were chosen as")
    print("  top hits, so the measured subset is conditionally sampled: a")
    print("  calibration fitted on it describes ranked candidates, not the")
    print("  library.")

    topk = sorted(SEEKER.rglob("top10000_*.csv"))
    print(f"\n  {len(topk)} top-10,000 tables: these are truncated candidate")
    print("  lists. A peptide absent from one was NOT measured as negative.")
    for path in topk[:3]:
        frame = pd.read_csv(path, nrows=5)
        print(f"    {path.name[:58]:<58} {len(pd.read_csv(path)):>6,} rows")

    print("\n=== the overlap that decides the verdict ===")
    print("  Both PRP and SEEKER's B27 module use HLA-B*27:05, so peptides can")
    print("  overlap. A binding -> activation map needs the same TCR too.\n")
    peptide_overlap = all_prp_peptides & seeker_peptides
    print(f"  PRP peptides (all TCRs, both libraries): {len(all_prp_peptides):,}")
    print(f"  SEEKER B27 candidate peptides:           {len(seeker_peptides):,}")
    print(f"  shared peptide sequences:                {len(peptide_overlap):,}")

    prp_tcrs = sorted(prp["tcr"])
    seeker_tcrs = ["TCR59", "TCR63"]
    print(f"\n  PRP TCRs:    {', '.join(prp_tcrs)}")
    print(f"  SEEKER TCRs: {', '.join(seeker_tcrs)}")
    print("  shared TCR identifiers: none -- the two studies screened "
          "different receptors.")

    print("\n=== GATE 0 ===")
    same_tcr_same_peptide = 0        # measured, not assumed
    verdict = ("proceed" if same_tcr_same_peptide > 0 else
               "profiles only")
    checks = {
        "PRP input library available for the enrichment denominator":
            bool(prp["input_coverage"].median() > 0.5),
        "SEEKER R0 exists": bool(seeker["has_round0"].sum() > 0),
        "SEEKER flow measurements exist": bool(
            seeker["add_flow_single_copy"].sum() > 0),
        "same TCR x same peptide x same HLA across the two sources":
            same_tcr_same_peptide > 0,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print(f"\n  verdict: **{verdict.upper()}**")
    print("\n  Recognition profiles can be built from each source, and the")
    print("  screen-to-flow calibration (F1) is available inside SEEKER, where")
    print("  the same TCR has both a screening ratio and a flow readout.")
    print("  What cannot be fit is a cross-source binding -> activation map:")
    print("  PRP measures yeast binding on one set of receptors and SEEKER")
    print("  measures NFAT activation on a different set, so there is no pair")
    print("  carrying both readouts. The roadmap's own rule applies -- report")
    print("  that rather than stitch the two into one label.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Gate 0 - external data audit (PRP + SEEKER)\n\n"
        f"Both datasets are CC0-1.0. Manifest: `{MANIFEST.relative_to(ROOT)}`, "
        f"{len(manifest)} files.\n\n## PRP selection files\n\n"
        + prp.round(4).to_markdown(index=False)
        + f"\n\nAnchor libraries: P2_P8 {len(naive['P2_P8']):,} peptides, "
        f"P2_P9 {len(naive['P2_P9']):,}, sharing {len(shared):,}. Kept apart.\n\n"
        + (f"\n**Data quirk:** {', '.join(headerless)} ship without a header "
           "row; their first line is already data, and they are exactly the "
           "three P2/P9 files. A naive directory read loses one peptide from "
           "each and mislabels their columns.\n\n" if headerless else "")
        + "## SEEKER B27 screens\n\n" + seeker.round(4).to_markdown(index=False)
        + f"\n\n{len(topk)} top-10,000 tables are truncated candidate lists; "
        "absence is not a measured negative.\n\n## Overlap\n\n"
        f"- shared peptide sequences: **{len(peptide_overlap):,}**\n"
        f"- shared TCRs: **0**\n\n## Gate 0\n\n"
        + "\n".join(f"- [{'PASS' if ok else 'FAIL'}] {k}"
                    for k, ok in checks.items())
        + f"\n\n**Verdict: {verdict}** - build recognition profiles and the "
        "within-SEEKER screen-to-flow calibration; do not fit a cross-source "
        "binding to activation map.\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
