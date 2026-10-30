#!/usr/bin/env python
"""E142 -- Averaging the hyperparameter instead of selecting it.

Every candidate since E96 tried to add information, and E131 explained why they
all failed: clearing a measured +0.005 needs a TRUE effect of +0.0086 and the
largest channel ever found here is +0.0016. That bound applies to *information*.

This experiment does not add information. It removes an estimation error.

D's four components are combined with **uniform** weights,
`(density3 + m30 + vdomain + structure) / 4`, so the combination layer already
averages. But the hyperparameters *inside* the components were **selected** on
the five public epitopes -- density3's per-space bandwidths were picked from a
nine-point grid, as were m30's and the V-domain's top-k. E15 measured the cost
of selecting on five epitopes at **-0.0048**, which is three times the largest
feature gain the project ever found.

In a regime where the candidates being selected among are statistically
indistinguishable -- which E134 established is exactly this regime, all seven
top submissions inside +-1 null sd -- the argmax of K noisy estimates regresses
on new data by roughly the selection bias, while the *mean* of the K does not.
So averaging should beat selecting, and the gain needs no new channel.

The test is symmetric and uses E118's better estimator. On each held-out
pseudo-epitope fold:

  SELECT   pick the bandwidth per space that maximises the other 17 folds,
           exactly reproducing the procedure D's configuration came from
  AVERAGE  take the mean rank over the whole grid, no selection at all

Both are scored on the held-out fold, so neither sees it.

**Pre-registered: +0.001 to +0.004**, the reasoning being that E15's -0.0048
is the recoverable ceiling and recovery is partial. A result at or below zero
means the bandwidth selection was never costing anything.

Run: .venv/bin/python scripts/121_average_not_select.py
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

REPORT = ROOT / "reports/99_average_not_select.md"
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([5., 8., 12., 20., 30., 50., 80., 120., 200.])}
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}   # what D actually uses
MIN_N = 30
PREDICTED = (0.001, 0.004)


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
    lab = labels.dropna(subset=["pos"]).copy()
    lab["pos"] = lab["pos"].astype(int)
    y = lab.set_index("tcr_id")["pos"]
    eps = sorted(test["epitope_aa"].unique())

    def rank(s):
        return hyp.rank_within(s.reindex(test["tcr_id"]), test)

    # ---- distance matrices per space per cohort ----
    print("=== loading distance matrices ===")
    Dm = {}
    for spc in ("default", "small"):
        c = np.load(ROOT / f"data/interim/sceptr/{spc}.npz", allow_pickle=True)
        W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
        po = {t: i for i, t in enumerate(c["tcr_id"])}
        Dm[spc] = {}
        for e, g in test.groupby("epitope_aa"):
            idx = [po[t] for t in g["tcr_id"]]
            M = 1.0 - W[idx] @ W[idx].T
            np.fill_diagonal(M, np.inf)
            Dm[spc][e] = (M, g["tcr_id"].to_numpy())
    Dm["tcrdist"] = {}
    for e, g in test.groupby("epitope_aa"):
        cc = np.load(ROOT / f"data/interim/tcrdist/{e}.npz", allow_pickle=True)
        pp = {t: i for i, t in enumerate(cc["tcr_id"])}
        r = [pp[t] for t in g["tcr_id"]]
        M = (cc["alpha"] + cc["beta"]).astype(float)[np.ix_(r, r)]
        np.fill_diagonal(M, np.inf)
        Dm["tcrdist"][e] = (M, g["tcr_id"].to_numpy())
    print(f"  three spaces x {len(eps)} cohorts")

    # ---- density at every bandwidth, pre-ranked ----
    print("\n=== building density at every bandwidth on the grid ===")
    DENS = {}
    for spc in GRID:
        for sg in GRID[spc]:
            v = pd.Series(index=test["tcr_id"], dtype=float)
            for e in eps:
                M, ids = Dm[spc][e]
                v.loc[ids] = np.exp(-M / sg).mean(axis=1)
            DENS[(spc, float(sg))] = rank(v)
        print(f"  {spc}: {len(GRID[spc])} bandwidths")

    # ---- the other three components, held fixed ----
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv"
                        ).set_index("tcr_id")
    structure = rank(boltz["pae_tcr_pep_min"])
    vdomain = rank(rival.vdomain_topk(test))
    sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                      ).set_index("tcr_id")["prediction"]
    frozen = decomp.submitted_score(test)
    d3_fixed = sum(DENS[(s, FIXED[s])] for s in GRID) / 3
    # recover D's m30 term by subtraction, which E123 verified to 1e-16
    m30 = 4 * rank(sub) - rank(d3_fixed) - vdomain - structure
    m30 = rank(m30)

    def assemble(d3):
        return rank((rank(d3) + m30 + vdomain + structure) / 4)

    # ---- the 18 pseudo-epitope folds (E118's estimator) ----
    m = test.merge(lab[["tcr_id", "pos"]], on="tcr_id")
    m["_k"] = m["cdr3_alpha_aa"] + "|" + m["cdr3_beta_aa"] + "|" + m["epitope_aa"]
    f = pd.read_csv(ROOT / "data/external/vdjdb/vdjdb-2023-06-01/vdjdb_full.txt",
                    sep="\t", low_memory=False)
    f = f[(f["species"] == "HomoSapiens") & f["cdr3.alpha"].notna()
          & f["cdr3.beta"].notna()].copy()
    f["_k"] = f["cdr3.alpha"] + "|" + f["cdr3.beta"] + "|" + f["antigen.epitope"]
    ref = (f[f["_k"].isin(set(m["_k"]))].dropna(subset=["reference.id"])
           .drop_duplicates("_k").set_index("_k")["reference.id"])
    m["ref"] = m["_k"].map(ref)
    folds = {}
    for (e, k), g in m.dropna(subset=["ref"]).groupby(["epitope_aa", "ref"]):
        if len(g) >= MIN_N and g["pos"].nunique() > 1:
            folds[(e, str(k))] = g["tcr_id"].tolist()
    print(f"\n=== {len(folds)} pseudo-epitope folds, "
          f"{sum(len(v) for v in folds.values())} TCRs ===")

    def fauc(s, ids):
        yy = y.reindex(ids).to_numpy()
        return roc_auc_score(yy, s.reindex(ids).to_numpy()) \
            if len(set(yy)) > 1 else np.nan

    # cache assembled predictions: one per (space, bandwidth) triple is too many,
    # so vary one space at a time around the fixed point, plus the full average
    cand = {}
    for spc in GRID:
        for sg in GRID[spc]:
            d3 = sum(DENS[(s, sg if s == spc else FIXED[s])] for s in GRID) / 3
            cand[(spc, float(sg))] = assemble(d3)
    d3_avg = sum(DENS[k] for k in DENS) / len(DENS)
    AVG = assemble(d3_avg)
    SEL = assemble(d3_fixed)

    print("\n=== 1. the two objects, scored on the public macro ===")
    per_sel = hyp.per_epitope(labels, SEL)
    per_avg = hyp.per_epitope(labels, AVG)
    pub = sorted(lab["epitope_aa"].unique())
    print(f"  {'epitope':<12}{'D (selected sigma)':>20}{'averaged sigma':>17}{'diff':>9}")
    for e in pub:
        print(f"  {e:<12}{per_sel[e]:20.4f}{per_avg[e]:17.4f}"
              f"{per_avg[e]-per_sel[e]:+9.4f}")
    ms = float(np.mean([per_sel[e] for e in pub]))
    ma = float(np.mean([per_avg[e] for e in pub]))
    print(f"  {'MACRO':<12}{ms:20.4f}{ma:17.4f}{ma-ms:+9.4f}")
    print("  (public is IN-SAMPLE for the selected sigma and out-of-sample for")
    print("   the average, so this comparison is biased AGAINST averaging)")

    print("\n=== 2. the symmetric out-of-fold test, 18 folds ===")
    print(f"  pre-registered: {PREDICTED[0]:+.3f} to {PREDICTED[1]:+.3f}")
    base = {k: fauc(SEL, v) for k, v in folds.items()}
    ks = [k for k in folds if not np.isnan(base[k])]
    cache = {c: {k: fauc(cand[c], folds[k]) for k in ks} for c in cand}
    avgf = {k: fauc(AVG, folds[k]) for k in ks}
    held_sel, held_avg = [], []
    for k in ks:
        trn = [x for x in ks if x != k]
        # SELECT: best single bandwidth config on the other 17 folds
        best = max(cand, key=lambda c: np.nanmean([cache[c][t] for t in trn]))
        held_sel.append(cache[best][k])
        held_avg.append(avgf[k])
    hs, ha = np.array(held_sel), np.array(held_avg)
    d = ha - hs
    print(f"  SELECT  (bandwidth chosen on 17 folds): {hs.mean():.4f}")
    print(f"  AVERAGE (no selection at all)         : {ha.mean():.4f}")
    print(f"  difference: {d.mean():+.4f} +- {d.std(ddof=1)/np.sqrt(len(d)):.4f}"
          f"   folds where averaging wins: {int((d>0).sum())}/{len(d)}")
    lo = d.mean() - 1.96*d.std(ddof=1)/np.sqrt(len(d))
    hi = d.mean() + 1.96*d.std(ddof=1)/np.sqrt(len(d))
    print(f"  95% CI [{lo:+.4f}, {hi:+.4f}]   "
          f"{'EXCLUDES' if lo>0 or hi<0 else 'includes'} zero")
    print(f"  prediction {PREDICTED[0]:+.3f}..{PREDICTED[1]:+.3f}: "
          f"{'HIT' if PREDICTED[0] <= d.mean() <= PREDICTED[1] else 'MISSED'}")

    print("\n=== 3. how much does the bandwidth matter at all? ===")
    allm = {c: np.nanmean([cache[c][k] for k in ks]) for c in cand}
    best_c = max(allm, key=allm.get)
    worst_c = min(allm, key=allm.get)
    print(f"  best  config {best_c}  {allm[best_c]:.4f}")
    print(f"  worst config {worst_c}  {allm[worst_c]:.4f}")
    print(f"  spread across the whole grid: {allm[best_c]-allm[worst_c]:.4f}")
    print(f"  D's own fixed point scores {np.nanmean([base[k] for k in ks]):.4f}")
    print(f"  the average scores          {ha.mean():.4f}")

    lines = ["# E142 — Averaging the hyperparameter instead of selecting it", "",
             "Generated by `scripts/121_average_not_select.py`.", "",
             "This adds no information; it removes an estimation error. D's four",
             "components already carry uniform weights, but the bandwidths inside",
             "them were selected on five epitopes, and E15 measured that selection",
             "cost at **−0.0048**.", "",
             f"Pre-registered: {PREDICTED[0]:+.3f} to {PREDICTED[1]:+.3f}.", "",
             "## Public macro (biased against averaging)", "",
             "| epitope | D (selected σ) | averaged σ | diff |",
             "|---|---:|---:|---:|"]
    for e in pub:
        lines.append(f"| {e} | {per_sel[e]:.4f} | {per_avg[e]:.4f} "
                     f"| {per_avg[e]-per_sel[e]:+.4f} |")
    lines += [f"| **MACRO** | **{ms:.4f}** | **{ma:.4f}** | **{ma-ms:+.4f}** |", "",
              "## Symmetric out-of-fold, 18 pseudo-epitope folds", "",
              f"- SELECT (bandwidth chosen on the other 17): **{hs.mean():.4f}**",
              f"- AVERAGE (no selection): **{ha.mean():.4f}**",
              f"- difference **{d.mean():+.4f}**, 95% CI [{lo:+.4f}, {hi:+.4f}], "
              f"averaging wins {int((d>0).sum())}/{len(d)} folds", "",
              f"Spread across the entire bandwidth grid: "
              f"{allm[best_c]-allm[worst_c]:.4f}.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
