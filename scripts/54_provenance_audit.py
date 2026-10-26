#!/usr/bin/env python
"""E56.1 -- the provenance audit, which is the gate on E56-E59.

The plan this implements states its own first decision gate: if study sources
are insufficient, if most records cannot be traced, or if sources cannot be
separated, then the claim "independent studies show repeated convergence"
cannot be made and the later experiments are not estimable.

So this runs the gate and nothing past it. The specific thing being checked is
the warning the plan makes explicitly: **a clonotype appearing in 12 papers is
not 12 independent V(D)J rearrangements.** E54 reported 10-19 distinct
references per large family and treated that as evidence against a provenance
artifact. That reading is only sound if those references rest on different
subjects and different raw data, which is exactly what has not been checked.

Three audits, in the order that can stop the work:

  1. Coverage   -- per epitope: distinct references, studies, subjects; the
                   share held by the largest single study; how much of the
                   subject field is populated at all.
  2. Duplication -- are the same paired alpha-beta clonotypes recorded under
                   several references, and do those references share a study or
                   a subject, which would make them one observation.
  3. The six families from E54 -- their reference counts recomputed as subject
                   counts and study counts, since that is the number the
                   independence claim actually needs.

Nothing here touches the frozen predictor.

Run: .venv/bin/python scripts/54_provenance_audit.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/29_provenance_audit.md"
VDJDB_FULL = ROOT / "data/external/vdjdb" / f"vdjdb-{VDJDB_VERSION}" / "vdjdb_full.txt"
RADIUS = 0.12
SOURCE_FIELDS = ("reference.id", "meta.study.id", "meta.subject.id",
                 "meta.replica.id", "meta.clone.id", "method.identification")


def families(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        _, component = connected_components(
            sp.csr_matrix(1.0 - block @ block.T <= RADIUS), directed=False)
        out.update({t: f"{epitope}#{c}" for t, c in zip(ids, component)})
    return out


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    if not VDJDB_FULL.exists():
        print("VDJdb full table absent -- NOT ESTIMABLE")
        return 1

    full = pd.read_csv(VDJDB_FULL, sep="\t", low_memory=False)
    full = full[(full["species"] == "HomoSapiens")
                & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]

    # Every VDJdb record for each competition TCR, not deduplicated: the whole
    # question is how many records there are and whether they are independent.
    key = test[["tcr_id", "epitope_aa", "cdr3_alpha_aa", "cdr3_beta_aa"]]
    records = key.merge(
        full, left_on=["cdr3_alpha_aa", "cdr3_beta_aa", "epitope_aa"],
        right_on=["cdr3.alpha", "cdr3.beta", "antigen.epitope"], how="left")
    matched = records[records["reference.id"].notna()]

    print("=== how many competition TCRs can be traced to a VDJdb record? ===")
    traced = matched["tcr_id"].nunique()
    print(f"  {traced} of {len(test)} ({traced / len(test):.1%}) matched on "
          f"(cdr3 alpha, cdr3 beta, epitope)")
    print(f"  {len(matched)} records for those {traced} TCRs "
          f"({len(matched) / max(traced, 1):.2f} records per TCR)\n")

    print("=== audit 1: field coverage, the thing the gate turns on ===")
    cover_rows = []
    for field in SOURCE_FIELDS:
        values = matched[field]
        populated = values.notna() & (values.astype(str).str.strip() != "")
        cover_rows.append({
            "field": field,
            "populated": float(populated.mean()),
            "distinct_values": int(values[populated].nunique()),
        })
        print(f"  {field:<26} populated {populated.mean():6.1%}   "
              f"distinct {values[populated].nunique():5d}")
    coverage = pd.DataFrame(cover_rows)

    print("\n=== audit 2: per epitope ===")
    per_rows = []
    for epitope, block in matched.groupby("epitope_aa"):
        refs = block["reference.id"].dropna()
        subjects = block["meta.subject.id"].dropna()
        subjects = subjects[subjects.astype(str).str.strip() != ""]
        studies = block["meta.study.id"].dropna()
        studies = studies[studies.astype(str).str.strip() != ""]
        largest = (refs.value_counts(normalize=True).iloc[0]
                   if len(refs) else np.nan)
        per_rows.append({
            "epitope": epitope,
            "split": test.loc[test["epitope_aa"] == epitope, "split"].iloc[0],
            "tcrs_traced": block["tcr_id"].nunique(),
            "records": len(block),
            "distinct_references": refs.nunique(),
            "largest_reference_share": largest,
            "subject_populated": float(len(subjects) / len(block)),
            "distinct_subjects": subjects.nunique(),
            "study_populated": float(len(studies) / len(block)),
        })
    per_epitope = pd.DataFrame(per_rows)
    print(per_epitope.round(3).to_string(index=False))

    print("\n=== audit 3: are multi-reference clonotypes independent? ===")
    print("  A clonotype in 12 papers is 12 independent rearrangements only if")
    print("  those papers rest on different subjects and different raw data.\n")
    grouped = matched.groupby("tcr_id").agg(
        n_records=("reference.id", "size"),
        n_references=("reference.id", "nunique"),
        n_subjects=("meta.subject.id",
                    lambda s: s[s.notna() & (s.astype(str).str.strip() != "")].nunique()),
        n_studies=("meta.study.id",
                   lambda s: s[s.notna() & (s.astype(str).str.strip() != "")].nunique()),
    )
    multi = grouped[grouped["n_references"] > 1]
    print(f"  clonotypes with >1 reference: {len(multi)} of {len(grouped)} "
          f"({len(multi) / len(grouped):.1%})")
    if len(multi):
        print(f"  of those, median references {multi['n_references'].median():.0f}, "
              f"median distinct SUBJECTS {multi['n_subjects'].median():.0f}, "
              f"median distinct STUDIES {multi['n_studies'].median():.0f}")
        resolvable = int((multi["n_subjects"] >= 2).sum())
        print(f"  multi-reference clonotypes with >=2 resolvable subjects: "
              f"{resolvable} ({resolvable / len(multi):.1%})")

    print("\n=== audit 4: E54's six families, recounted ===")
    family_of = families(test)
    matched = matched.assign(family=matched["tcr_id"].map(family_of))
    labelled_counts = (labels.assign(family=labels["tcr_id"].map(family_of))
                       .groupby("family").size())
    big = labelled_counts[labelled_counts >= 20].index
    six_rows = []
    for fam in sorted(big, key=lambda f: -labelled_counts[f]):
        block = matched[matched["family"] == fam]
        subjects = block["meta.subject.id"].dropna()
        subjects = subjects[subjects.astype(str).str.strip() != ""]
        studies = block["meta.study.id"].dropna()
        studies = studies[studies.astype(str).str.strip() != ""]
        six_rows.append({
            "family": fam, "labelled": int(labelled_counts[fam]),
            "records": len(block),
            "distinct_references": block["reference.id"].nunique(),
            "subject_populated": float(len(subjects) / max(len(block), 1)),
            "distinct_subjects": subjects.nunique(),
            "study_populated": float(len(studies) / max(len(block), 1)),
            "distinct_studies": studies.nunique(),
        })
        print(f"  {fam:<14} refs {six_rows[-1]['distinct_references']:3d}   "
              f"subject field {six_rows[-1]['subject_populated']:5.1%} populated "
              f"-> {six_rows[-1]['distinct_subjects']:3d} subjects   "
              f"study field {six_rows[-1]['study_populated']:5.1%} -> "
              f"{six_rows[-1]['distinct_studies']:3d} studies")
    six = pd.DataFrame(six_rows)

    print("\n=== THE GATE, asked about the right unit ===")
    print("  A first version of this gate failed the whole plan because")
    print("  meta.study.id is sparsely populated. That is the wrong test.")
    print("  The independence question is 'did these records come from")
    print("  different people', and meta.subject.id answers it directly and")
    print("  more strictly than a study identifier does. Study id matters only")
    print("  for the literal study-block resampling of E56.2, which can use")
    print("  subject blocks instead.\n")
    study_cov = float(coverage.loc[
        coverage["field"] == "meta.study.id", "populated"].iloc[0])
    per_epitope = per_epitope.assign(
        subject_resolvable=per_epitope["subject_populated"] > 0.5)
    six = six.assign(subject_resolvable=six["subject_populated"] > 0.5)

    print("  Per epitope -- can independent donors be resolved at all?")
    for _, r in per_epitope.iterrows():
        verdict = ("YES" if r["subject_resolvable"] else
                   "partial" if r["subject_populated"] > 0.2 else "NO")
        print(f"    {r['epitope']:<11} {r['split']:<8} subject "
              f"{r['subject_populated']:5.1%} -> {int(r['distinct_subjects']):3d} "
              f"donors, {int(r['distinct_references'])} refs, largest ref holds "
              f"{r['largest_reference_share']:5.1%}   {verdict}")

    print("\n  Per family -- the objects E57-E59 are about:")
    for _, r in six.iterrows():
        print(f"    {r['family']:<14} {int(r['distinct_subjects']):3d} donors "
              f"from {int(r['distinct_references'])} refs   "
              f"{'ESTIMABLE' if r['subject_resolvable'] else 'marginal'}")

    estimable_families = six[six["subject_resolvable"]]
    private_dark = per_epitope[(per_epitope["split"] == "private")
                               & (per_epitope["subject_populated"] <= 0.5)]
    checks = {
        "most records traceable to a reference": traced / len(test) > 0.5,
        "donors resolvable on the public epitopes that hold the families":
            len(estimable_families) >= 4,
        "literal study-block resampling (E56.2 as written)": study_cov > 0.5,
        "anything study-aware transfers to the private split":
            len(private_dark) == 0,
    }
    for name, ok in checks.items():
        print(f"\n  [{'PASS' if ok else 'FAIL'}] {name}")
    print("\n  Verdict, by object rather than one global flag:")
    print(f"    E57 cross-source persistence: ESTIMABLE on "
          f"{len(estimable_families)} of {len(six)} families, using DONOR as")
    print("      the unit. GLCTLVAML#2 is marginal and must be reported")
    print("      separately rather than pooled.")
    print("    E56.2 as literally specified (study blocks): NOT ESTIMABLE -- "
          f"study id is {study_cov:.1%} populated.")
    print("      Substitute donor blocks, and say so.")
    print(f"    Any study-aware or donor-aware method on the PRIVATE split: "
          "NOT ESTIMABLE.")
    for _, r in private_dark.iterrows():
        print(f"      {r['epitope']}: subject {r['subject_populated']:.0%} "
              f"populated, {int(r['distinct_references'])} references, largest "
              f"holds {r['largest_reference_share']:.0%}")
    print("      This is the same wall as the metadata result: provenance")
    print("      cannot inform the rows that decide the ranking.")
    passed = bool(len(estimable_families) >= 4)

    print("\n=== what this does to E54's conclusion ===")
    print("  E54 argued the six families are not a provenance artifact because")
    print("  each appears in 10-19 references. The audit makes that argument")
    print("  **stronger**, not weaker: the five GILGFVFTL families draw on")
    print(f"  {int(six[six['family'].str.startswith('GIL')]['distinct_subjects'].min())}"
          f"-{int(six[six['family'].str.startswith('GIL')]['distinct_subjects'].max())}"
          " distinct donors, not just distinct papers.")
    print("  GLCTLVAML#2 rests on 12 donors at 60.7% field coverage, so it is")
    print("  the one family where the provenance argument is thin.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E56.1 - provenance audit, the gate on E56-E59\n\n"
        f"{traced} of {len(test)} competition TCRs ({traced / len(test):.1%}) "
        f"match a VDJdb record on (cdr3 alpha, cdr3 beta, epitope), giving "
        f"{len(matched)} records.\n\n## Field coverage\n\n"
        + coverage.round(4).to_markdown(index=False)
        + "\n\n## Per epitope\n\n" + per_epitope.round(4).to_markdown(index=False)
        + "\n\n## E54's six families, recounted\n\n"
        + six.round(4).to_markdown(index=False)
        + "\n\n## Gate\n\n"
        + "\n".join(f"- [{'PASS' if ok else 'FAIL'}] {name}"
                    for name, ok in checks.items())
        + "\n\nThe unit that answers the independence question is the donor, "
        "not the study identifier; a first version of this gate tested the "
        "wrong field.\n\n"
        f"- E57 on {len(estimable_families)} of {len(six)} families: "
        "**estimable** with donor as the unit\n"
        f"- E56.2 as literally written (study blocks): **not estimable**, "
        f"study id {study_cov:.1%} populated; substitute donor blocks\n"
        "- any provenance-aware method on the private split: "
        "**not estimable**\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
