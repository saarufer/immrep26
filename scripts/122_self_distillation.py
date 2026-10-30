#!/usr/bin/env python
"""E144 -- Transductive self-distillation of D.

The last untried route. E131 bounded adding information (+0.0086 needed,
+0.0016 found), E142 bounded bandwidth averaging (flat surface), E143 found
robust aggregation EV-neutral. This is the third member of the
remove-estimation-error family and it operates at a layer the other two do not.

**The mechanism, stated so it can fail.** A student trained on a teacher's
outputs cannot acquire information the teacher lacks -- by construction it only
approximates D. Its one route to beating D is **self-distillation**: D's
density terms are k-NN estimates over small cohorts (CINGVCWTV 243 rows with a
pool of 11, LTDEMIAQY 130 with a pool of 5), so they carry estimation noise,
and a smoother hypothesis class fitted to D's own ranking can average that
noise out. If D's ranking is already smooth in sequence space, the student can
only lose.

This is transductive and legitimately so: D's ranks exist for all 2,636 rows
including all 602 private ones, and **the student never sees a true label**.
Every student prediction is cross-fitted, so no row is scored by a model that
was trained on it.

Evaluation is against the true public labels and E118's 18 pseudo-epitope
folds, neither of which entered training.

**Pre-registered:** student alone **-0.01 to -0.05** on the public macro, since
approximation can only lose; blend with D **-0.002 to +0.002**. A blend above
+0.002 would mean the smoothing is real.

Run: .venv/bin/python scripts/122_self_distillation.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/100_self_distillation.md"
SEED = 0
NFOLD = 5
MIN_N = 30
PREDICTED_STUDENT = (-0.05, -0.01)
PREDICTED_BLEND = (-0.002, 0.002)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    test = load_test().reset_index(drop=True)
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    lab = labels.dropna(subset=["pos"]).copy()
    lab["pos"] = lab["pos"].astype(int)
    y = lab.set_index("tcr_id")["pos"]
    pub = sorted(lab["epitope_aa"].unique())
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    # ---- the teacher's target: D's within-epitope percentile rank ----
    t = pd.Series(index=test["tcr_id"], dtype=float)
    for e, g in test.groupby("epitope_aa"):
        t.loc[g["tcr_id"]] = dr.reindex(g["tcr_id"]).rank(pct=True).to_numpy()

    # ---- the student's features: sequence only, available everywhere ----
    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {x: i for i, x in enumerate(c["tcr_id"])}
    X = W[[po[x] for x in test["tcr_id"]]]
    extra = np.column_stack([
        test["cdr3_alpha_aa"].fillna("").str.len().to_numpy(),
        test["cdr3_beta_aa"].fillna("").str.len().to_numpy()])
    # epitope as one-hot, so the student can carry cohort offsets
    oh = pd.get_dummies(test["epitope_aa"]).to_numpy().astype(float)
    X = np.hstack([X, extra, oh])
    print(f"=== student features: {X.shape[1]} dims "
          f"({W.shape[1]} SCEPTR + 2 lengths + {oh.shape[1]} epitope) ===")
    print(f"  target: D's within-epitope percentile rank, no true labels used")

    # ---- cross-fitted students ----
    tv = t.to_numpy()
    MODELS = {
        "ridge": lambda: RidgeCV(alphas=np.logspace(-3, 3, 25)),
        "gbm": lambda: HistGradientBoostingRegressor(
            max_iter=300, learning_rate=0.05, max_depth=4,
            random_state=SEED),
    }
    S = {}
    for nm, mk in MODELS.items():
        pred = np.zeros(len(test))
        for tr, te in KFold(NFOLD, shuffle=True, random_state=SEED).split(X):
            m = mk()
            m.fit(X[tr], tv[tr])
            pred[te] = m.predict(X[te])
        S[nm] = pd.Series(pred, index=test["tcr_id"])
        rho = np.mean([pd.Series(pred, index=test["tcr_id"]
                                 ).reindex(g["tcr_id"]).corr(
            t.reindex(g["tcr_id"]), method="spearman")
            for _, g in test.groupby("epitope_aa")])
        print(f"  {nm:<6} cross-fitted, within-epitope rho with the teacher "
              f"{rho:+.3f}")

    # ---- 18 pseudo-epitope folds ----
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

    def fauc(s, ids):
        yy = y.reindex(ids).to_numpy()
        return roc_auc_score(yy, s.reindex(ids).to_numpy()) \
            if len(set(yy)) > 1 else np.nan

    def macro(s):
        p = hyp.per_epitope(labels, hyp.rank_within(s, test))
        return float(np.mean([p[e] for e in pub]))

    base_m = macro(dr)
    base_f = np.nanmean([fauc(dr, v) for v in folds.values()])
    print(f"\n=== D, the teacher: public macro {base_m:.4f}   "
          f"18-fold {base_f:.4f} ===")
    print(f"  pre-registered student {PREDICTED_STUDENT[0]:+.3f}.."
          f"{PREDICTED_STUDENT[1]:+.3f}, blend "
          f"{PREDICTED_BLEND[0]:+.3f}..{PREDICTED_BLEND[1]:+.3f}")

    print("\n=== 1. the student alone ===")
    print(f"  {'model':<8}{'public macro':>14}{'vs D':>9}{'18-fold':>9}{'vs D':>9}")
    rows = []
    for nm, s in S.items():
        mm, ff = macro(s), np.nanmean([fauc(hyp.rank_within(s, test), v)
                                       for v in folds.values()])
        print(f"  {nm:<8}{mm:14.4f}{mm-base_m:+9.4f}{ff:9.4f}{ff-base_f:+9.4f}")
        rows.append([nm, "alone", mm - base_m, ff - base_f])

    print("\n=== 2. blended with D ===")
    WS = [0.05, 0.10, 0.20, 0.30, 0.50]
    print(f"  {'model':<8}{'w':>6}{'public macro':>14}{'vs D':>9}"
          f"{'18-fold':>9}{'vs D':>9}")
    for nm, s in S.items():
        sr = hyp.rank_within(s, test)
        for w in WS:
            b = (1 - w) * dr + w * sr
            mm = macro(b)
            ff = np.nanmean([fauc(hyp.rank_within(b, test), v)
                             for v in folds.values()])
            print(f"  {nm:<8}{w:6.2f}{mm:14.4f}{mm-base_m:+9.4f}"
                  f"{ff:9.4f}{ff-base_f:+9.4f}")
            rows.append([nm, f"w={w:.2f}", mm - base_m, ff - base_f])

    print("\n=== 3. where it would matter: the pool-starved cohorts ===")
    print("  D's density is noisiest where the reference pool is smallest, so")
    print("  that is where smoothing should pay. Public proxies by pool size:")
    pool = {"GILGFVFTL": 3005, "GLCTLVAML": 459, "NLVPMVATV": 9881,
            "TTDPSFLGRY": 25, "YLQPRTFLL": 54}
    best = max(((r[2], r) for r in rows if r[1] != "alone"), key=lambda x: x[0])[1]
    nm, wlab = best[0], best[1]
    w = float(wlab.split("=")[1])
    b = (1 - w) * dr + w * hyp.rank_within(S[nm], test)
    pD = hyp.per_epitope(labels, dr)
    pB = hyp.per_epitope(labels, hyp.rank_within(b, test))
    print(f"  using the best blend ({nm}, {wlab}):")
    print(f"  {'epitope':<12}{'pool':>7}{'D':>9}{'blend':>9}{'diff':>9}")
    for e in sorted(pub, key=lambda x: pool[x]):
        print(f"  {e:<12}{pool[e]:7d}{pD[e]:9.4f}{pB[e]:9.4f}{pB[e]-pD[e]:+9.4f}")

    best_blend = max(r[2] for r in rows if r[1] != "alone")
    worst_alone = min(r[2] for r in rows if r[1] == "alone")
    print(f"\n=== verdict ===")
    print(f"  student alone, worst: {worst_alone:+.4f}  -> prediction "
          f"{PREDICTED_STUDENT[0]:+.3f}..{PREDICTED_STUDENT[1]:+.3f} "
          f"{'HIT' if PREDICTED_STUDENT[0] <= worst_alone <= PREDICTED_STUDENT[1] else 'MISSED'}")
    print(f"  best blend:           {best_blend:+.4f}  -> prediction "
          f"{PREDICTED_BLEND[0]:+.3f}..{PREDICTED_BLEND[1]:+.3f} "
          f"{'HIT' if PREDICTED_BLEND[0] <= best_blend <= PREDICTED_BLEND[1] else 'MISSED'}")

    lines = ["# E144 — Transductive self-distillation of D", "",
             "Generated by `scripts/122_self_distillation.py`.", "",
             "The student learns D's within-epitope percentile rank from sequence",
             "features only, cross-fitted over 5 folds, and **never sees a true",
             "label**. It cannot acquire information D lacks; its only route to",
             "beating D is smoothing D's own k-NN estimation noise.", "",
             f"D: public macro **{base_m:.4f}**, 18-fold **{base_f:.4f}**.", "",
             "| model | setting | public vs D | 18-fold vs D |",
             "|---|---|---:|---:|"]
    for nm, st, dm, df in rows:
        lines.append(f"| `{nm}` | {st} | {dm:+.4f} | {df:+.4f} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
