#!/usr/bin/env python
"""E145 -- pDockQ and ipSAE: the structural quality-tier channel.

From a literature sweep: "Structural quality-tier assessment for TCR-pMHC
functional enrichment" (PMC13447355) trains a random forest on 1,325
AlphaFold3 models to assign a quality tier, and reports that **excluding
low-quality models raises precision to 0.72 on TCRvdb** -- validated on
GLCTLVAML and YLQPRTFLL, two of our own public epitopes. Its features are
pLDDT, ipTM, iPAE, iPDE, **ipSAE** and **pDockQ v1/v2**.

Against the admission standard:

  Q1  partly new. `data/interim/boltz_features.csv` already holds 35 columns of
      ipTM/pLDDT/PAE, and that channel is measured: 0.7387 on task A but
      **0.5362** on the real task, with `iptm_mhc_pep` -- an interface holding
      no TCR -- at 0.647 and `iptm` at **-0.0193** as a residual. But pDockQ
      and ipSAE are *different composites*: pDockQ weights interface pLDDT by
      the contact count, and ipSAE restricts pTM to interface pairs under a PAE
      cutoff. Neither is in our feature set.
  Q2  **passes, and rarely** -- it needs only structures, and we hold PAE
      matrices, per-residue pLDDT and PDB files for all 2,636 rows including
      all 602 private ones (15,816 files = 2,636 x 6 models).
  Q3  passes. The tier model is trained on *structural* quality (DockQ against
      reference structures), not on binding labels, so it is independent of the
      VDJdb annotation process this task audits -- unlike E121's Hi-TpH.
  Q4  information.

**Pre-registered residual-on-D: 0.51-0.55**, i.e. near the 0.5098 floor, with
the reasoning on record -- these are composites of quantities already measured
at -0.0193. The genuine difference is the *threshold* framing, so a binarised
tier-style variant is tested alongside the continuous scores.

Run: .venv/bin/python scripts/123_pdockq_ipsae.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

PRED = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
REPORT = ROOT / "reports/101_pdockq_ipsae.md"
CACHE = ROOT / "data/interim/pdockq_ipsae.csv"
CUT = 8.0          # interface contact cutoff, Angstrom
PAE_CUT = 15.0     # ipSAE PAE cutoff
D0 = 10.0          # pDockQ2 PAE scale
PREDICTED = (0.51, 0.55)
TCR, PMHC = ("A", "B"), ("M", "P")


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def parse_pdb(p: Path):
    xyz, ch, res = [], [], []
    for line in p.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            ch.append(line[21])
            res.append(int(line[22:26]))
    return np.array(xyz), np.array(ch), np.array(res)


def metrics(pdb: Path, pae: np.ndarray, plddt: np.ndarray):
    xyz, ch, _ = parse_pdb(pdb)
    n = len(xyz)
    if n == 0 or pae.shape[0] < n or len(plddt) < n:
        return None
    pae = pae[:n, :n]
    pl = np.asarray(plddt[:n], dtype=float)
    if pl.max() <= 1.5:
        pl = pl * 100.0
    t = np.isin(ch, TCR)
    m = np.isin(ch, PMHC)
    if t.sum() == 0 or m.sum() == 0:
        return None
    ti, mi = np.where(t)[0], np.where(m)[0]
    # interface residue pairs within CUT
    tree = cKDTree(xyz[mi])
    pairs = tree.query_ball_point(xyz[ti], r=CUT)
    rows, cols = [], []
    for a, lst in zip(ti, pairs):
        for b in lst:
            rows.append(a)
            cols.append(mi[b])
    if not rows:
        return None
    rows = np.array(rows)
    cols = np.array(cols)
    ncon = len(rows)
    iface = np.unique(np.concatenate([rows, cols]))
    ipl = float(pl[iface].mean())
    # pDockQ v1 (Bryant 2022 constants)
    x = ipl * np.log10(max(ncon, 1))
    pdq1 = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018
    # pDockQ2: interface pLDDT weighted by PAE-derived pair confidence
    pw = 1.0 / (1.0 + (pae[rows, cols] / D0) ** 2)
    x2 = ipl * float(pw.mean())
    pdq2 = 0.724 / (1 + np.exp(-0.052 * (x2 - 152.611))) + 0.018
    # ipSAE-style: pTM restricted to interface pairs under a PAE cutoff
    keep = pae[rows, cols] < PAE_CUT
    nk = max(int(keep.sum()), 1)
    d0 = max(1.24 * (nk - 15) ** (1 / 3) - 1.8, 0.5) if nk > 15 else 0.5
    ptm_terms = 1.0 / (1.0 + (pae[rows, cols][keep] / d0) ** 2)
    ipsae = float(ptm_terms.mean()) if keep.any() else 0.0
    return dict(n_contacts=ncon, iface_plddt=ipl, pdockq1=float(pdq1),
                pdockq2=float(pdq2), ipsae=ipsae,
                iface_pae_mean=float(pae[rows, cols].mean()),
                frac_pae_confident=float(keep.mean()))


def main() -> int:
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    lab = scored.set_index("tcr_id")["pos"]
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    if CACHE.exists():
        M = pd.read_csv(CACHE).set_index("tcr_id")
        print(f"=== cached: {len(M)} rows ===")
    else:
        print(f"=== computing pDockQ / ipSAE from {PRED} ===")
        rec = {}
        dirs = sorted(d for d in PRED.iterdir() if d.is_dir())
        print(f"  {len(dirs)} prediction directories")
        for i, d in enumerate(dirs):
            pdbs = sorted(d.glob("*_model_0.pdb"))
            paes = sorted(d.glob("pae_*_model_0.npz"))
            plds = sorted(d.glob("plddt_*_model_0.npz"))
            if not (pdbs and paes and plds):
                continue
            try:
                pae = np.load(paes[0])["pae"]
                pl = np.load(plds[0])["plddt"]
            except Exception:
                continue
            r = metrics(pdbs[0], pae, pl)
            if r:
                rec[d.name] = r
            if (i + 1) % 400 == 0:
                print(f"  {i+1}/{len(dirs)}  parsed {len(rec)}")
        M = pd.DataFrame(rec).T
        M.index.name = "tcr_id"
        M.to_csv(CACHE)
        print(f"  wrote {CACHE.relative_to(ROOT)}  ({len(M)} rows)")

    cov = test["tcr_id"].isin(M.index)
    print(f"\n=== coverage: {int(cov.sum())} of {len(test)} rows ===")
    for e, g in test.groupby("epitope_aa"):
        c = g["tcr_id"].isin(M.index)
        print(f"  {e:<12}{g['split'].iloc[0]:<9}{int(c.sum()):4d}/{len(g):4d}")

    cols = ["pdockq1", "pdockq2", "ipsae", "iface_plddt", "n_contacts",
            "iface_pae_mean", "frac_pae_confident"]
    F = {c: M[c].reindex(test["tcr_id"]).astype(float) for c in cols if c in M}
    # tier-style binarisation, the paper's actual framing
    if "pdockq1" in F:
        q = F["pdockq1"]
        F["tier_not_LQ (pdockq1>median)"] = (q > q.median()).astype(float)

    print(f"\n=== gates; floor 0.5098, record 0.5609, "
          f"pre-registered {PREDICTED[0]:.2f}-{PREDICTED[1]:.2f} ===")
    print(f"  {'quantity':<30}{'macro':>8}{'typ rho':>9}{'priv':>6}"
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
        print(f"  {nm:<30}{max(mm,1-mm):8.4f}{tg['rho_with_typicality']:+9.3f}"
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
            print(f"  {nm:<30}{len(a):2d} strata  mean AUC {np.mean(a):.4f}"
                  f"  above 0.5: {sum(1 for v in a if v>0.5)}/{len(a)}")

    print("\n=== blended into D, nested LOEO ===")
    pub = sorted(scored["epitope_aa"].unique())
    base = hyp.mean_auc(labels, dr)
    WS = [0.0, 0.05, 0.10, 0.15, 0.20]
    print(f"  D alone {base:.4f}")
    for nm, s in F.items():
        v = s.copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.median()).to_numpy()
        v = v.fillna(v.median())
        best = None
        for sg in (1.0, -1.0):
            sr = hyp.rank_within(sg * v, test)
            per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
            vals = {w: float(np.mean([per[w][e] for e in pub])) for w in WS}
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pick][e])
            bw = max(WS, key=lambda w: vals[w])
            if best is None or vals[bw] > best[1]:
                best = (sg, vals[bw], float(np.mean(held)), bw)
        sg, ov, nest, bw = best
        print(f"  {nm:<30}sign {int(sg):+d} w={bw:.2f}  oracle {ov:.4f} "
              f"({ov-base:+.4f})  nested {nest:.4f} ({nest-base:+.4f})")

    lines = ["# E145 — pDockQ and ipSAE, the structural quality-tier channel", "",
             "Generated by `scripts/123_pdockq_ipsae.py`.", "",
             f"Pre-registered residual-on-D: {PREDICTED[0]:.2f}–{PREDICTED[1]:.2f}.",
             "", "| quantity | macro | residual on D | leak | valid | private var |",
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
