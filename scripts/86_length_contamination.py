#!/usr/bin/env python
"""Is CDR3 length contaminating the within-epitope distance, and can that be used?

E95 established that CDR3 length separates the *cohorts* — mean CDR3α runs
11.55 on LLWNGPMAV to 14.20 on GILGFVFTL, η² 0.171 across the six HLA-A2
epitopes. That is a between-epitope fact and the competition scores
within-epitope, so it looked like a dead end. It is not, because length also
varies **inside** a cohort (SD 1.7–2.2 residues), and the density that D runs
on is a cosine distance in an isotropic 64-dimensional space. If a large share
of that distance is simply "these two CDR3s are different lengths", then
typicality is partly **length commonness**: a TCR looks typical because most of
its cohort happens to share its length.

Two consequences, pulling opposite ways, and both are tested:

  1. **Length commonness as its own feature.** It is a within-cohort quantity,
     so unlike E95's between-cohort lengths it is not structurally barred, and
     unlike every provenance variable it varies inside every cohort.
  2. **Length-controlled density.** If length is contamination, removing it
     should sharpen what remains. Two removals: standardising each TCR's
     density against same-length TCRs, and residualising the distance matrix
     itself against |Δlength| before any density is computed.

Standing gates apply: `private_variance_gate` (E3), `typicality_gate` (0.6378),
the pre-registered **+0.02** inclusion bar, and nested LOEO against E15's
−0.0048 selection cost.

Run: .venv/bin/python scripts/86_length_contamination.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/63_length_contamination.md"
SIGMA = 0.1
TOPK = 8
BAR = 0.02


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    test = load_test().copy()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    V = cached["vectors"]
    V = V / np.linalg.norm(V, axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(cached["tcr_id"])}
    test["la"] = test["cdr3_alpha_aa"].str.len().astype(float)
    test["lb"] = test["cdr3_beta_aa"].str.len().astype(float)
    epitopes = sorted(test["epitope_aa"].unique())

    feats: dict[str, dict] = {}

    def put(name, mapping):
        feats[name] = mapping

    # ---------- build every variant in one pass over the cohorts ----------
    contamination = []
    for e, g in test.groupby("epitope_aa"):
        idx = np.array([po[t] for t in g["tcr_id"]])
        X = V[idx]
        D = 1.0 - X @ X.T
        n = len(g)
        la, lb = g["la"].to_numpy(), g["lb"].to_numpy()
        dla = np.abs(la[:, None] - la[None, :])
        dlb = np.abs(lb[:, None] - lb[None, :])
        iu = np.triu_indices(n, 1)
        A = np.column_stack([dla[iu], dlb[iu], np.ones(len(iu[0]))])
        beta, *_ = np.linalg.lstsq(A, D[iu], rcond=None)
        pred = A @ beta
        r2 = 1 - ((D[iu] - pred) ** 2).sum() / ((D[iu] - D[iu].mean()) ** 2).sum()
        contamination.append([e, g["split"].mode().iloc[0], n, float(r2)])

        # length-residualised distance: subtract the |Δlen|-explained component
        full = np.column_stack([dla.ravel(), dlb.ravel(), np.ones(n * n)])
        Dres = D - (full @ beta).reshape(n, n)
        Dres = Dres - Dres.min() + 1e-6
        np.fill_diagonal(Dres, np.inf)

        Dd = D.copy()
        np.fill_diagonal(Dd, np.inf)
        dens = -np.sort(Dd, axis=1)[:, :TOPK].mean(axis=1)
        dens_res = -np.sort(Dres, axis=1)[:, :TOPK].mean(axis=1)

        # length commonness, three forms
        pa = g["la"].map(g["la"].value_counts(normalize=True)).to_numpy()
        pb = g["lb"].map(g["lb"].value_counts(normalize=True)).to_numpy()
        joint = (g.groupby(["la", "lb"])["tcr_id"].transform("size") / n).to_numpy()
        # smoothed: adjacent lengths are similar, so a Gaussian over length
        def smooth(vals, bw=1.0):
            d = np.abs(vals[:, None] - vals[None, :])
            return np.exp(-(d ** 2) / (2 * bw ** 2)).mean(axis=1)
        sm = smooth(la) * smooth(lb)

        # density standardised within a length stratum
        key = pd.Series([f"{int(a)}_{int(b)}" for a, b in zip(la, lb)])
        dser = pd.Series(dens)
        mu = dser.groupby(key.values).transform("mean").to_numpy()
        sd = dser.groupby(key.values).transform("std").fillna(0).to_numpy()
        strat = np.where(sd > 1e-9, (dens - mu) / (sd + 1e-9), 0.0)
        # smoother version: regress density on the two lengths, keep residual
        B = np.column_stack([la, lb, la ** 2, lb ** 2, la * lb, np.ones(n)])
        cf, *_ = np.linalg.lstsq(B, dens, rcond=None)
        dens_lenresid = dens - B @ cf

        ids = g["tcr_id"].to_numpy()
        for nm, arr in (("density (reference)", dens),
                        ("lenfreq product", pa * pb),
                        ("lenfreq joint", joint),
                        ("lenfreq smoothed", sm),
                        ("density on length-residualised distance", dens_res),
                        ("density z within length stratum", strat),
                        ("density minus quadratic length fit", dens_lenresid)):
            feats.setdefault(nm, {}).update(dict(zip(ids, arr)))

    print("=== 1. how much of the within-epitope distance is |Δlength|? ===")
    print(f"  {'epitope':<12}{'split':<9}{'n':>5}{'R2':>8}")
    for c in sorted(contamination, key=lambda r: -r[3]):
        print(f"  {c[0]:<12}{c[1]:<9}{c[2]:5d}{c[3]:8.3f}")
    print(f"  {'MEAN':<12}{'':9}{'':5}"
          f"{np.mean([c[3] for c in contamination]):8.3f}")

    # ---------- 2. every variant, scored ----------
    print("\n=== 2. each variant on task A, with the standing gates ===")
    print(f"  {'feature':<42}{'macro':>8}{'rho w/ typ':>12}{'typ gate':>10}"
          f"{'priv gate':>11}")
    series = {k: pd.Series(v) for k, v in feats.items()}
    rows = []
    for nm, s in series.items():
        g = typicality_gate(test, scored, s)
        pv, pvdetail = private_variance_gate(test, s)
        macro = hyp.mean_auc(labels, s)
        rho = g["rho_with_typicality"]
        print(f"  {nm:<42}{macro:8.4f}{rho:+12.3f}"
              f"{'pass' if abs(rho) < 0.6378 else 'FAIL':>10}"
              f"{'pass' if pv else 'FAIL':>11}")
        rows.append([nm, macro, rho, abs(rho) < 0.6378, pv, pvdetail])

    # ---------- 3. blend the survivors into D ----------
    print("\n=== 3. blending into D ===")
    d = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d, test)
    base = hyp.mean_auc(labels, dr)
    print(f"  D alone {base:.4f}")
    blend_rows = []
    for nm, s in series.items():
        if nm == "density (reference)":
            continue
        sr = hyp.rank_within(s, test)
        ws, vals = [0.05, 0.10, 0.15, 0.20, 0.30], []
        for w in ws:
            vals.append(hyp.mean_auc(labels, (1 - w) * dr + w * sr))
        bi = int(np.argmax(vals))
        # nested LOEO over the weight grid including 0
        held = []
        for e in epitopes:
            tr = [x for x in epitopes if x != e]
            def sc(w, eps):
                b = (1 - w) * dr + w * sr
                pe = hyp.per_epitope(labels, b)
                return float(np.mean([pe[t] for t in eps if t in pe]))
            bw = max([0.0] + ws, key=lambda w: sc(w, tr))
            pe = hyp.per_epitope(labels, (1 - bw) * dr + bw * sr)
            if e in pe:
                held.append(pe[e])
        nested = float(np.mean(held))
        print(f"  {nm:<42} oracle w={ws[bi]:.2f} {vals[bi]:.4f} "
              f"({vals[bi]-base:+.4f})   nested {nested:.4f} ({nested-base:+.4f})")
        blend_rows.append([nm, ws[bi], vals[bi], vals[bi] - base, nested,
                           nested - base])

    # ---------- report ----------
    lines = ["# E96 — CDR3 length contamination of the within-epitope distance",
             "", "Generated by `scripts/86_length_contamination.py`.", "",
             "## 1. How much of the within-epitope cosine distance is |Δlength|?",
             "", "Linear fit of pairwise cosine distance on |Δ CDR3α| and "
             "|Δ CDR3β|, within each cohort.", "",
             "| epitope | split | n | R² |", "|---|---|---:|---:|"]
    for c in sorted(contamination, key=lambda r: -r[3]):
        lines.append(f"| {c[0]} | {c[1]} | {c[2]} | **{c[3]:.3f}** |")
    lines += ["", f"Mean **{np.mean([c[3] for c in contamination]):.3f}**.", "",
              "## 2. Each variant on task A", "",
              "| feature | macro AUC | ρ with typicality | typicality gate | private gate |",
              "|---|---:|---:|---|---|"]
    for r in rows:
        lines.append(f"| {r[0]} | **{r[1]:.4f}** | {r[2]:+.3f} "
                     f"| {'pass' if r[3] else 'FAIL'} | {'pass' if r[4] else 'FAIL'} ({r[5]}) |")
    lines += ["", "## 3. Blending into D", "",
              f"D alone **{base:.4f}**. Bar is +{BAR:.2f}; E15's selection cost "
              "is −0.0048.", "",
              "| feature | oracle w | oracle macro | vs D | nested LOEO | vs D |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in blend_rows:
        lines.append(f"| {r[0]} | {r[1]:.2f} | {r[2]:.4f} | {r[3]:+.4f} "
                     f"| {r[4]:.4f} | **{r[5]:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
