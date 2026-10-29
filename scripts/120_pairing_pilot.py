#!/usr/bin/env python
"""E140 -- Pilot: does cell-level alphabeta pairing ambiguity explain revalidation failure?

The small pre-experiment a collaborator proposed, and the only form of their
E139 that can be completed before the deadline. E139's table-level test already
showed that `clonotype_aa_screened` equals the stored pair in 2,636/2,636 rows
and that chain-sharing multiplicity has the *wrong* sign (+0.3089) because it
is convergent recurrence. What that test could not see is the deep form: a cell
expressing two productive alpha chains where only one reached VDJdb.

The 10x Genomics "A New Way of Exploring Immunity" application note is the one
public source that exposes it. Per barcode it gives every productive contig
(so a second alpha chain is directly visible) and a dextramer call with UMI
counts (so capture strength is visible too, which is their E140). Four donors,
**189,512 barcodes**, covering GILGFVFTL and GLCTLVAML -- two public epitopes,
so labels exist.

Three questions, in order:

  1. how many receptors actually carry cell-level pairing ambiguity?
  2. does that ambiguity predict revalidation FAILURE, as hypothesised?
  3. does it survive holding D fixed?

**Pre-registered, before the run.** Their hypothesis predicts dual-alpha cells
fail revalidation more, so P(valid | dual-alpha) - P(valid | single-alpha)
should be clearly negative; I set the bar at <= -0.10 to count as support.
My own prediction is **-0.05 to +0.05**, i.e. null, because E139 found the
observable chain-sharing signal is recurrence rather than error. The cell-level
measurement is genuinely different from the table-level one, so this is a real
test and not a foregone conclusion.

For capture strength (their E140) the prediction is the opposite: supporting
cell count and dextramer UMI should associate *positively* with revalidation,
but be absorbed by D, because both are recurrence proxies and E29 established D
already exploits clonal recurrence (masking families costs -0.0282).

Stratified evaluation, not residual regression: these are small-count discrete
quantities and E112's tie artefact has now fired four times, most recently
producing four spurious residuals above the project record in E139.

Run: .venv/bin/python scripts/120_pairing_pilot.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

TENX = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
            "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/tenx")
REPORT = ROOT / "reports/98_pairing_pilot.md"
DEX = {"GILGFVFTL": "A0201_GILGFVFTL_Flu-MP_Influenza",
       "GLCTLVAML": "A0201_GLCTLVAML_BMLF1_EBV"}
SUPPORT_BAR = -0.10


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    lab = scored.set_index("tcr_id")["pos"]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    print("=== loading the 10x application note ===")
    CT, BM = [], []
    for d in (1, 2, 3, 4):
        c = pd.read_csv(TENX / f"d{d}_all_contig_annotations.csv", low_memory=False)
        c["donor"] = d
        c["bc"] = f"d{d}:" + c["barcode"].astype(str)
        CT.append(c[["bc", "donor", "is_cell", "high_confidence", "chain",
                     "productive", "cdr3", "umis"]])
        keep = ["barcode"] + [v for v in DEX.values()
                              if v in pd.read_csv(
                                  TENX / f"d{d}_binarized_matrix.csv",
                                  nrows=0).columns]
        keep += [v + "_binder" for v in DEX.values()]
        b = pd.read_csv(TENX / f"d{d}_binarized_matrix.csv",
                        usecols=lambda x: x in keep, low_memory=False)
        b["bc"] = f"d{d}:" + b["barcode"].astype(str)
        BM.append(b)
    C = pd.concat(CT, ignore_index=True)
    B = pd.concat(BM, ignore_index=True)
    print(f"  contigs  {len(C):,}   barcodes {C.bc.nunique():,}")
    print(f"  dextramer rows {len(B):,}")

    def truthy(s):
        return s.astype(str).str.lower().isin(["true", "yes", "1"])

    C = C[truthy(C.is_cell) & truthy(C.high_confidence) & truthy(C.productive)]
    C = C[C.cdr3.notna() & (C.cdr3 != "None")]
    print(f"  productive high-confidence contigs in cells: {len(C):,}"
          f"  cells {C.bc.nunique():,}")

    # ---- per-cell chain inventory ----
    inv = (C.groupby(["bc", "chain"])
             .agg(n=("cdr3", "nunique"),
                  cdr3s=("cdr3", lambda x: "|".join(sorted(set(x)))),
                  umi=("umis", "sum")).reset_index())
    piv_n = inv.pivot(index="bc", columns="chain", values="n").fillna(0)
    piv_s = inv.pivot(index="bc", columns="chain", values="cdr3s")
    cells = pd.DataFrame({"n_TRA": piv_n.get("TRA", 0),
                          "n_TRB": piv_n.get("TRB", 0),
                          "TRA": piv_s.get("TRA"), "TRB": piv_s.get("TRB")})
    cells = cells[(cells.n_TRA >= 1) & (cells.n_TRB >= 1)]
    print(f"\n=== 1. how common is cell-level pairing ambiguity? ===")
    print(f"  cells with >=1 productive TRA and TRB: {len(cells):,}")
    tot = len(cells)
    for k in (1, 2, 3):
        na = int((cells.n_TRA == k).sum())
        nb = int((cells.n_TRB == k).sum())
        print(f"    {k} productive TRA: {na:7,} ({na/tot*100:5.2f}%)"
              f"     {k} productive TRB: {nb:7,} ({nb/tot*100:5.2f}%)")
    print(f"    >=2 TRA (dual alpha): {int((cells.n_TRA>=2).sum()):,} "
          f"({(cells.n_TRA>=2).mean()*100:.2f}%)")
    print(f"    >=2 TRB:              {int((cells.n_TRB>=2).sum()):,} "
          f"({(cells.n_TRB>=2).mean()*100:.2f}%)")

    # ---- match competition rows to source cells, per epitope ----
    B = B.set_index("bc")
    rows = []
    print("\n=== 2. matching competition rows to their source cells ===")
    for ep, col in DEX.items():
        bcol = col + "_binder"
        if bcol not in B.columns:
            print(f"  {ep}: no binder column, skipped")
            continue
        pos_bc = set(B.index[truthy(B[bcol])])
        sub = cells[cells.index.isin(pos_bc)]
        comp = test[test.epitope_aa == ep]
        # explode the multi-chain cells so a dual-alpha cell matches on either
        rec = []
        for bc, r in sub.iterrows():
            for a in str(r.TRA).split("|"):
                for b in str(r.TRB).split("|"):
                    rec.append((a, b, bc, r.n_TRA, r.n_TRB))
        M = pd.DataFrame(rec, columns=["a", "b", "bc", "n_TRA", "n_TRB"])
        key = comp["cdr3_alpha_aa"] + "|" + comp["cdr3_beta_aa"]
        M["k"] = M.a + "|" + M.b
        g = M.groupby("k").agg(cells=("bc", "nunique"),
                               max_nTRA=("n_TRA", "max"),
                               max_nTRB=("n_TRB", "max"))
        hit = comp.assign(k=key.to_numpy()).merge(
            g, left_on="k", right_index=True, how="inner")
        print(f"  {ep:<12}dextramer-positive cells {len(sub):6,}"
              f"   competition rows matched {len(hit):4d} of {len(comp)}")
        for _, r in hit.iterrows():
            rows.append([r.tcr_id, ep, int(r.cells), int(r.max_nTRA),
                         int(r.max_nTRB)])
    P = pd.DataFrame(rows, columns=["tcr_id", "epitope", "n_cells",
                                    "max_nTRA", "max_nTRB"])
    P = P.drop_duplicates("tcr_id")
    P["y"] = lab.reindex(P.tcr_id).to_numpy()
    P = P.dropna(subset=["y"])
    P["y"] = P.y.astype(int)
    print(f"\n  matched rows with a visible label: {len(P)}")
    if len(P) < 20:
        print("  STOP: too few matched rows to test anything.")
        return 0

    print(f"\n=== 3. does dual-alpha predict revalidation FAILURE? ===")
    print(f"  pre-registered: hypothesis supported only if diff <= {SUPPORT_BAR:+.2f};"
          f" my prediction -0.05 to +0.05")
    dual = P.max_nTRA >= 2
    print(f"  {'group':<28}{'n':>5}{'P(valid)':>10}")
    print(f"  {'dual alpha (>=2 TRA)':<28}{int(dual.sum()):5d}"
          f"{P.loc[dual,'y'].mean() if dual.any() else np.nan:10.3f}")
    print(f"  {'single alpha':<28}{int((~dual).sum()):5d}"
          f"{P.loc[~dual,'y'].mean():10.3f}")
    diff = (P.loc[dual, "y"].mean() - P.loc[~dual, "y"].mean()) if dual.any() else np.nan
    print(f"  difference: {diff:+.4f}   -> "
          f"{'SUPPORTS the hypothesis' if diff <= SUPPORT_BAR else 'does NOT support it'}")
    dualb = P.max_nTRB >= 2
    if dualb.any():
        print(f"  (dual beta, for symmetry: n={int(dualb.sum())}, "
              f"P(valid) {P.loc[dualb,'y'].mean():.3f} against "
              f"{P.loc[~dualb,'y'].mean():.3f})")

    print(f"\n=== 4. capture strength (their E140) ===")
    for c in ("n_cells", "max_nTRA", "max_nTRB"):
        if P[c].nunique() < 2:
            continue
        a = roc_auc_score(P.y, P[c])
        print(f"  {c:<14}AUC {a:.4f}   median {P[c].median():.1f}  "
              f"range {P[c].min()}-{P[c].max()}")

    print(f"\n=== 5. does any of it survive holding D fixed? ===")
    P["d"] = dr.reindex(P.tcr_id).to_numpy()
    print("  (stratified, not residual regression -- E112's tie artefact has")
    print("   now fired four times on small-count features)")
    for c in ("n_cells", "max_nTRA"):
        aucs = []
        for ep, g in P.groupby("epitope"):
            if len(g) < 20:
                continue
            q = pd.qcut(g.d.rank(method="first"),
                        min(4, max(2, len(g)//10)), labels=False)
            for k in sorted(set(q)):
                s = g[q == k]
                if s.y.nunique() < 2 or s[c].nunique() < 2:
                    continue
                aucs.append(roc_auc_score(s.y, s[c]))
        if aucs:
            print(f"  {c:<14}{len(aucs):2d} strata   mean AUC {np.mean(aucs):.4f}"
                  f"   above 0.5: {sum(1 for a in aucs if a>0.5)}/{len(aucs)}")
        else:
            print(f"  {c:<14}no usable strata")

    lines = ["# E140 — Pilot: cell-level αβ pairing ambiguity", "",
             "Generated by `scripts/120_pairing_pilot.py`.", "",
             "Source: 10x Genomics *A New Way of Exploring Immunity*, 4 donors,",
             f"{C.bc.nunique():,} cells with productive high-confidence contigs.", "",
             "## 1. How common is cell-level pairing ambiguity?", "",
             f"- cells with ≥1 productive TRA and TRB: **{len(cells):,}**",
             f"- **≥2 productive TRA (dual α): {int((cells.n_TRA>=2).sum()):,} "
             f"({(cells.n_TRA>=2).mean()*100:.2f}%)**",
             f"- ≥2 productive TRB: {int((cells.n_TRB>=2).sum()):,} "
             f"({(cells.n_TRB>=2).mean()*100:.2f}%)", "",
             "## 2–3. Does it predict revalidation failure?", "",
             f"Matched competition rows with a label: **{len(P)}**.", "",
             "| group | n | P(valid) |", "|---|---:|---:|",
             f"| dual α (≥2 TRA) | {int(dual.sum())} | "
             f"{P.loc[dual,'y'].mean() if dual.any() else float('nan'):.3f} |",
             f"| single α | {int((~dual).sum())} | {P.loc[~dual,'y'].mean():.3f} |",
             "", f"Difference **{diff:+.4f}**, against a pre-registered support "
             f"bar of ≤ {SUPPORT_BAR:+.2f}.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
