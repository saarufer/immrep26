#!/usr/bin/env python
"""Hubness: a different functional of the same distance matrix.

E98 concluded that the ceiling is not caused by the candidate features
overlapping — three orthogonal mechanisms summed at 98% of additive and still
reached only +0.0009. So the next thing worth trying has to be outside the
family, and the family is narrower than it looks: every density variant this
project has built is the same functional, *mean distance to my k nearest
neighbours*, computed in a different space.

Hubness is a different functional of the identical matrix, and a well known
pathology of k-NN in high dimensions (Radovanović et al. 2010). For a point x:

  out-degree view (what D uses)  mean distance from x to its k nearest
  in-degree view  N_k(x)         how many points have x among *their* k nearest

These are correlated but not monotone. A **hub** sits in everyone's neighbour
list while its own neighbours are far; an **antihub** is in nobody's list. In
a TCR cohort the hub is the generic receptor that is vaguely close to
everything and specifically close to nothing — which is a different thing from
being in a sparse neighbourhood, and plausibly a different thing for
revalidation.

Three standard corrections are tested as drop-in replacements for the distance
that density is computed on:

  **local scaling (NICDM)**  d'(x,y) = d(x,y) / sqrt(r_k(x) · r_k(y))
  **mutual proximity**       MP(x,y) = P(d(x,·) > d(x,y)) · P(d(y,·) > d(x,y))
  **mutual k-NN**            keep only edges where each is in the other's top-k

Local scaling is worth stating plainly: it is **a per-TCR adaptive bandwidth**.
E97 established that the bandwidth should vary by cohort and bought +0.0007
inside D for it; this is the same idea at the level of the individual receptor,
and it has never been tried here.

The standing gates apply, plus the +0.02 bar and nested LOEO against E15's
−0.0048. E29 is the relevant warning: masking each TCR's own family costs
−0.0282, so clonal recurrence *is* signal, and any correction that flattens
families should be expected to pay for it.

Run: .venv/bin/python scripts/89_hubness.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/66_hubness.md"
K = 8
SIGMA = 0.1


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    epitopes = sorted(test["epitope_aa"].unique())
    pub = sorted(scored["epitope_aa"].unique())

    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}

    feats: dict[str, dict] = {}
    hub_stats = []
    for e, g in test.groupby("epitope_aa"):
        ids = g["tcr_id"].to_numpy()
        X = W[[po[t] for t in ids]]
        D = 1.0 - X @ X.T
        n = len(ids)
        np.fill_diagonal(D, np.inf)
        order = np.argsort(D, axis=1)
        knn = order[:, :K]

        # --- the two degree views ---
        dens = -np.sort(D, axis=1)[:, :K].mean(axis=1)
        indeg = np.bincount(knn.ravel(), minlength=n).astype(float)

        # --- local scaling (NICDM): per-TCR adaptive bandwidth ---
        rk = np.sort(D, axis=1)[:, K - 1]
        Dls = D / np.sqrt(np.outer(rk, rk) + 1e-12)
        np.fill_diagonal(Dls, np.inf)
        dens_ls = -np.sort(Dls, axis=1)[:, :K].mean(axis=1)

        # --- mutual proximity: P(both find each other far) ---
        finite = np.where(np.isinf(D), np.nan, D)
        rank_in_row = np.argsort(np.argsort(np.nan_to_num(finite, nan=1e9),
                                            axis=1), axis=1) / (n - 1)
        MP = 1.0 - (1 - rank_in_row) * (1 - rank_in_row.T)
        np.fill_diagonal(MP, np.inf)
        dens_mp = -np.sort(MP, axis=1)[:, :K].mean(axis=1)

        # --- mutual kNN: only reciprocal edges count ---
        inknn = np.zeros((n, n), dtype=bool)
        np.put_along_axis(inknn, knn, True, axis=1)
        mutual = inknn & inknn.T
        Dmk = np.where(mutual, D, np.inf)
        cnt = mutual.sum(axis=1)
        with np.errstate(invalid="ignore"):
            dens_mk = np.where(cnt > 0,
                               -np.where(mutual, D, 0).sum(axis=1)
                               / np.maximum(cnt, 1), -1.0)
        dens_mk = dens_mk + 0.01 * cnt      # reward reciprocity count

        sk = float(pd.Series(indeg).skew())
        hub_stats.append([e, g["split"].mode().iloc[0], n, sk,
                          float((indeg == 0).mean()), float(indeg.max()),
                          float(spearmanr(dens, indeg).statistic)])

        for nm, arr in (("density (reference)", dens),
                        ("in-degree N_k", indeg),
                        ("density on local scaling (NICDM)", dens_ls),
                        ("density on mutual proximity", dens_mp),
                        ("mutual-kNN density", dens_mk),
                        ("in-degree minus density rank",
                         pd.Series(indeg).rank(pct=True).to_numpy()
                         - pd.Series(dens).rank(pct=True).to_numpy())):
            feats.setdefault(nm, {}).update(dict(zip(ids, arr)))

    print("=== 1. is there hubness here at all? ===")
    print(f"  {'epitope':<12}{'split':<9}{'n':>5}{'N_k skew':>10}"
          f"{'antihub %':>11}{'max N_k':>9}{'rho(dens,N_k)':>15}")
    for h in hub_stats:
        print(f"  {h[0]:<12}{h[1]:<9}{h[2]:5d}{h[3]:10.2f}{h[4]*100:11.1f}"
              f"{h[5]:9.0f}{h[6]:+15.3f}")

    print("\n=== 2. each variant on task A ===")
    series = {k: pd.Series(v) for k, v in feats.items()}
    print(f"  {'feature':<36}{'macro':>8}{'rho w/ typ':>12}{'typ':>6}{'priv':>7}")
    rows = []
    for nm, s in series.items():
        g = typicality_gate(test, scored, s)
        pv, pvd = private_variance_gate(test, s)
        m = hyp.mean_auc(labels, s)
        rho = g["rho_with_typicality"]
        print(f"  {nm:<36}{m:8.4f}{rho:+12.3f}"
              f"{'pass' if abs(rho) < 0.6378 else 'FAIL':>6}"
              f"{'pass' if pv else 'FAIL':>7}")
        rows.append([nm, m, rho, abs(rho) < 0.6378, pv])

    # how decorrelated is each from the density it would replace?
    print("\n=== 3. decorrelation from the reference density ===")
    ref = series["density (reference)"]
    dec = []
    for nm, s in series.items():
        if nm == "density (reference)":
            continue
        r = float(np.mean([spearmanr(s.reindex(g["tcr_id"]),
                                     ref.reindex(g["tcr_id"])).statistic
                           for _, g in test.groupby("epitope_aa")]))
        print(f"  {nm:<36}{r:+8.3f}")
        dec.append([nm, r])

    print("\n=== 4. blended into D ===")
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)
    base = hyp.mean_auc(labels, dr)
    print(f"  D alone {base:.4f}")
    WS = [0.05, 0.10, 0.15, 0.20, 0.30]
    blend = []
    for nm, s in series.items():
        if nm == "density (reference)":
            continue
        sr = hyp.rank_within(s, test)
        per_w = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
        per_w[0.0] = hyp.per_epitope(labels, dr)
        vals = {w: float(np.mean([per_w[w][e] for e in pub])) for w in per_w}
        bw = max(WS, key=lambda w: vals[w])
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            pick = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
            held.append(per_w[pick][e])
        nested = float(np.mean(held))
        print(f"  {nm:<36} oracle w={bw:.2f} {vals[bw]:.4f} ({vals[bw]-base:+.4f})"
              f"   nested {nested:.4f} ({nested-base:+.4f})")
        blend.append([nm, bw, vals[bw], vals[bw] - base, nested, nested - base])

    lines = ["# E99 — Hubness: a different functional of the same distance matrix",
             "", "Generated by `scripts/89_hubness.py`.", "",
             "## 1. Is there hubness here?", "",
             "| epitope | split | n | N_k skew | antihub % | max N_k | ρ(density, N_k) |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for h in hub_stats:
        lines.append(f"| {h[0]} | {h[1]} | {h[2]} | {h[3]:.2f} | {h[4]*100:.1f} "
                     f"| {h[5]:.0f} | {h[6]:+.3f} |")
    lines += ["", "## 2. Each variant on task A", "",
              "| feature | macro | ρ with typicality | typicality gate | private gate |",
              "|---|---:|---:|---|---|"]
    for r in rows:
        lines.append(f"| {r[0]} | **{r[1]:.4f}** | {r[2]:+.3f} "
                     f"| {'pass' if r[3] else 'FAIL'} | {'pass' if r[4] else 'FAIL'} |")
    lines += ["", "## 3. Decorrelation from the density it would replace", "",
              "| feature | within-epitope ρ with density |", "|---|---:|"]
    for nm, r in dec:
        lines.append(f"| {nm} | {r:+.3f} |")
    lines += ["", "## 4. Blended into D", "", f"D alone **{base:.4f}**.", "",
              "| feature | oracle w | oracle | vs D | nested LOEO | vs D |",
              "|---|---:|---:|---:|---:|---:|"]
    for b in blend:
        lines.append(f"| {b[0]} | {b[1]:.2f} | {b[2]:.4f} | {b[3]:+.4f} "
                     f"| {b[4]:.4f} | **{b[5]:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
