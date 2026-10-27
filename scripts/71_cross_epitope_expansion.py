#!/usr/bin/env python
"""E80 -- cross-epitope reference expansion constrained by peptide similarity.

Every pool experiment here (E12, E27, E33, E35) used same-epitope TCRs. This
uses TCRs annotated to different but similar peptides.

Run: .venv/bin/python scripts/71_cross_epitope_expansion.py
"""
from __future__ import annotations
import importlib.util, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/46_cross_epitope_expansion.md"
MAX_HAMMING, K = 4, 8


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 99


def main() -> int:
    spec = importlib.util.spec_from_file_location(
        "pool", ROOT / "scripts/33_pool_negative_control.py")
    pool_mod = importlib.util.module_from_spec(spec)
    sys.modules["pool"] = pool_mod
    spec.loader.exec_module(pool_mod)
    import sceptr
    sceptr.disable_hardware_acceleration()
    test, labels = load_test(), build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    full = pd.read_csv(ROOT / f"data/external/vdjdb/vdjdb-{VDJDB_VERSION}/vdjdb_full.txt",
                       sep="\t", low_memory=False)
    full = full[(full["species"] == "HomoSapiens")
                & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]

    print("=== availability: how many similar-peptide records exist? ===")
    rows = []
    for epitope, cohort in test.groupby("epitope_aa"):
        allele = "A\\*02" if cohort["hla_short"].iloc[0] == "A2" else "A\\*01"
        same_class = full[full["mhc.a"].astype(str).str.contains(allele)]
        near = sorted({p for p in same_class["antigen.epitope"].unique()
                       if isinstance(p, str) and p != epitope
                       and hamming(epitope, p) <= MAX_HAMMING})
        records = int(same_class["antigen.epitope"].isin(near).sum())
        rows.append({"target": epitope, "split": cohort["split"].iloc[0],
                     "similar_peptides": len(near), "records": records,
                     "nearest": min(near, key=lambda p: hamming(epitope, p))
                     if near else "-"})
        print(f"  {epitope:<11} {cohort['split'].iloc[0]:<8} {len(near):>3} peptides, "
              f"{records:>5} records, nearest {rows[-1]['nearest']}")
    availability = pd.DataFrame(rows)
    private = availability[availability["split"] == "private"]["records"].sum()
    print(f"\n  the three private epitopes can gain {private} records in total")

    print("\n=== NLVPMVATV: the only evaluable case ===")
    target = "NLVPMVATV"
    cohort = test[test["epitope_aa"] == target].reset_index(drop=True)
    query = pd.DataFrame({"cdr3a": cohort["cdr3_alpha_aa"], "cdr3b": cohort["cdr3_beta_aa"],
                          "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
                          "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"],
                          "ep": target})
    qv = pool_mod.embed(query)
    a2 = full[full["mhc.a"].astype(str).str.contains("A\\*02")]
    near = sorted({p for p in a2["antigen.epitope"].unique()
                   if isinstance(p, str) and p != target
                   and hamming(target, p) <= MAX_HAMMING})
    print("  nearest: " + ", ".join(
        f"{p}(d={hamming(target, p)})"
        for p in sorted(near, key=lambda p: hamming(target, p))[:4]))
    keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
    extra = a2[a2["antigen.epitope"].isin(near)].drop_duplicates(
        ["cdr3.alpha", "cdr3.beta"])
    extra = extra[[(a, b) not in keys for a, b in
                   zip(extra["cdr3.alpha"], extra["cdr3.beta"])]]
    print(f"  {int(a2['antigen.epitope'].isin(near).sum())} records collapse to "
          f"**{len(extra)}** new paired clonotypes after deduplication")
    xv = pool_mod.embed(pd.DataFrame({
        "cdr3a": extra["cdr3.alpha"], "cdr3b": extra["cdr3.beta"],
        "trav": extra["v.alpha"], "traj": extra["j.alpha"],
        "trbv": extra["v.beta"], "trbj": extra["j.beta"], "ep": "NEAR"}))
    block = labels[labels["epitope_aa"] == target]

    def score(reference):
        d = 1.0 - qv @ reference.T
        d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
        series = pd.Series(-np.sort(d, axis=1)[:, :K].mean(axis=1),
                           index=cohort["tcr_id"])
        return roc_auc_score(block["pos"], block["tcr_id"].map(series))
    before, after = score(qv), score(np.vstack([qv, xv]))
    print(f"  cohort alone {before:.4f} -> with similar-peptide TCRs "
          f"{after:.4f}  ({after - before:+.4f})")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("# E80 - cross-epitope expansion by peptide similarity\n\n"
                      + availability.to_markdown(index=False)
                      + f"\n\nNLVPMVATV: {len(extra)} new clonotypes, "
                      f"{before:.4f} -> {after:.4f}.\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
