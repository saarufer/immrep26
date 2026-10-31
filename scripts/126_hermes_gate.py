#!/usr/bin/env python
"""E148 part 2 -- gating the HERMES compatibility score.

Reads `data/interim/hermes_compat.csv` produced by scripts/125 and puts dLL
through the standing gates, with the controls the external review asked for.

Pre-registered in scripts/125 before the run: residual-on-D **0.52-0.57**,
nested **-0.002 to +0.002**.

Run: .venv/bin/python scripts/126_hermes_gate.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/103_hermes_compatibility.md"
PREDICTED_RESID = (0.52, 0.57)
PREDICTED_NESTED = (-0.002, 0.002)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    lab = scored.set_index("tcr_id")["pos"]
    pub = sorted(scored["epitope_aa"].unique())
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    H = pd.read_csv(ROOT / "data/interim/hermes_compat.csv").set_index("tcr_id")
    print(f"=== {len(H)} rows, {int(H.dll.notna().sum())} finite ===")
    cov = test["tcr_id"].isin(H.index[H.dll.notna()])
    print(f"  coverage of the 2,636 competition rows: {int(cov.sum())}")
    for e, g in test.groupby("epitope_aa"):
        c = g["tcr_id"].isin(H.index[H.dll.notna()])
        print(f"  {e:<12}{g['split'].iloc[0]:<9}{int(c.sum()):4d}/{len(g):4d}")

    F = {}
    for c in ("dll", "dll_mean", "ll_full", "ll_pmhc"):
        if c in H:
            F[c] = H[c].reindex(test["tcr_id"]).astype(float)
    # controls the review asked for
    Q = ROOT / "data/interim/pdockq_ipsae.csv"
    CTRL = {}
    if Q.exists():
        P = pd.read_csv(Q).set_index("tcr_id")
        for c in ("n_contacts", "iface_plddt", "pdockq1"):
            if c in P:
                CTRL[c] = P[c].reindex(test["tcr_id"]).astype(float)
    E = ROOT / "data/interim/interface_energetics.csv"
    if E.exists():
        Eg = pd.read_csv(E).set_index("tcr_id")
        for c in ("bsa", "polar_frac"):
            if c in Eg:
                CTRL[c] = Eg[c].reindex(test["tcr_id"]).astype(float)

    print(f"\n=== is dLL just contact count, burial or confidence? ===")
    print(f"  (the review's control; within-epitope Spearman)")
    for c, s in CTRL.items():
        r = float(np.mean([spearmanr(F["dll"].reindex(g["tcr_id"]),
                                     s.reindex(g["tcr_id"]),
                                     nan_policy="omit").statistic
                           for _, g in test.groupby("epitope_aa")]))
        print(f"  rho(dll, {c:<14}) = {r:+.3f}")

    print(f"\n=== gates; floor 0.5098, record 0.5609, "
          f"pre-registered {PREDICTED_RESID[0]:.2f}-{PREDICTED_RESID[1]:.2f} ===")
    print(f"  {'quantity':<14}{'macro':>8}{'typ rho':>9}{'priv':>6}"
          f"{'resid':>8}{'leak':>7}{'valid':>7}")
    rows = []
    for nm, s in F.items():
        v = s.copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.median()).to_numpy()
        v = v.fillna(v.median())
        tg = typicality_gate(test, scored, v)
        pv, _ = private_variance_gate(test, v)
        mm = hyp.mean_auc(labels, v)
        res = pd.Series(index=v.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in v.index]
            xx = dr.reindex(i).rank(pct=True).to_numpy()
            yy = v.reindex(i).rank(pct=True).to_numpy()
            r = yy - np.polyval(np.polyfit(xx, yy, 1), xx)
            res.loc[i] = r
            leak.append(abs(spearmanr(r, xx).statistic))
        rr = hyp.mean_auc(labels, res)
        resid, mx = max(rr, 1 - rr), max(leak)
        print(f"  {nm:<14}{max(mm,1-mm):8.4f}{tg['rho_with_typicality']:+9.3f}"
              f"{'pass' if pv else 'FAIL':>6}{resid:8.4f}{mx:7.3f}"
              f"{'yes' if mx<0.1 else 'NO':>7}")
        rows.append([nm, max(mm, 1 - mm), resid, mx, pv])

    print("\n=== leak-free: within D-quintile strata ===")
    for nm, s in F.items():
        a = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in lab.index
                 and np.isfinite(s.get(t, np.nan))]
            if len(i) < 40:
                continue
            q = pd.qcut(dr.reindex(i).rank(method="first"), 4, labels=False)
            for k in range(4):
                sel = [t for t, qq in zip(i, q) if qq == k]
                yy = lab.reindex(sel).to_numpy()
                xx = s.reindex(sel).to_numpy()
                if len(set(yy)) < 2 or len(set(xx)) < 2:
                    continue
                a.append(roc_auc_score(yy, xx))
        if a:
            from scipy.stats import binomtest
            above = sum(1 for v in a if v > 0.5)
            p = binomtest(above, len(a), 0.5, alternative="greater").pvalue
            print(f"  {nm:<14}{len(a):2d} strata  mean AUC {np.mean(a):.4f}"
                  f"  above 0.5: {above}/{len(a)}  binom p {p:.3f}")

    print("\n=== blended into D, nested LOEO ===")
    base = hyp.mean_auc(labels, dr)
    WS = [0.0, 0.05, 0.10, 0.15, 0.20]
    print(f"  D alone {base:.4f}   pre-registered nested "
          f"{PREDICTED_NESTED[0]:+.3f} to {PREDICTED_NESTED[1]:+.3f}")
    best_nested = None
    for nm, s in F.items():
        v = s.copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.median()).to_numpy()
        v = v.fillna(v.median())
        for sg in (1.0, -1.0):
            sr = hyp.rank_within(sg * v, test)
            per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
            vals = {w: float(np.mean([per[w][e] for e in pub])) for w in WS}
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pick][e])
            nest = float(np.mean(held))
            bw = max(WS, key=lambda w: vals[w])
            if sg > 0 or vals[bw] - base > 0:
                print(f"  {nm:<14}sign {int(sg):+d} w={bw:.2f}  oracle "
                      f"{vals[bw]:.4f} ({vals[bw]-base:+.4f})  nested "
                      f"{nest:.4f} ({nest-base:+.4f})")
            if best_nested is None or nest - base > best_nested:
                best_nested = nest - base

    print(f"\n=== verdict ===")
    rr = [r for r in rows if r[0] == "dll"][0]
    print(f"  dll residual {rr[2]:.4f}  -> pre-registration "
          f"{'HIT' if PREDICTED_RESID[0] <= rr[2] <= PREDICTED_RESID[1] else 'MISSED'}")
    print(f"  best nested {best_nested:+.4f}  -> pre-registration "
          f"{'HIT' if PREDICTED_NESTED[0] <= best_nested <= PREDICTED_NESTED[1] else 'MISSED'}")

    lines = ["# E148 — HERMES interface compatibility", "",
             "Generated by `scripts/125_hermes_compatibility.py` and",
             "`scripts/126_hermes_gate.py`.", "",
             "`dLL = Σ_peptide [ log P(aa | full TCR–pMHC) − log P(aa | pMHC only) ]`,",
             "HERMES `hermes_bp_000` 10-member ensemble, ensembled at logits level,",
             "hydrogens placed by `reduce` 4.17 built from source (matching the",
             "released preprocessing), seed 0.", "",
             f"Pre-registered: residual {PREDICTED_RESID[0]:.2f}–{PREDICTED_RESID[1]:.2f}, "
             f"nested {PREDICTED_NESTED[0]:+.3f} to {PREDICTED_NESTED[1]:+.3f}.", "",
             "| quantity | macro | residual on D | leak | valid | private var |",
             "|---|---:|---:|---:|---|---|"]
    for nm, mm, resid, mx, pv in rows:
        lines.append(f"| `{nm}` | {mm:.4f} | **{resid:.4f}** | {mx:.3f} "
                     f"| {'yes' if mx<0.1 else 'NO'} "
                     f"| {'pass' if pv else 'FAIL'} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
