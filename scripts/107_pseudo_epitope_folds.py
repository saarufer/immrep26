#!/usr/bin/env python
"""Sixteen folds instead of five: can more units of selection resolve the gain?

E117 located the bottleneck precisely. TCR-level sampling noise is not the
problem — the paired family bootstrap puts the frozen configuration's advantage
at +0.0074 ± 0.0040. The problem is that **the fold is the unit of selection
and there are only five of them**, so a nested search produces up to +0.0203 on
shuffled labels (E106).

More real epitopes would fix it and are not obtainable:

  - **IMMREP25** is a different task — "previously unpublished TCR-epitope
    *binding* data" against this year's "whether pairs reported in public
    databases can be *validated*". Its epitopes carry no revalidation labels.
  - **A VDJdb replication proxy**, using E90's finding that multi-study
    annotations validate at odds ratio 17, needs both single- and
    multi-reference TCRs within one epitope. Of 21 VDJdb epitopes with ≥80
    paired-chain TCRs, exactly **one** has ≥15 on each side.
  - New revalidation labels require wet-lab work.

What is obtainable is **pseudo-epitopes**: split each cohort by the study that
contributed it. E90 established that studies differ sharply — same-epitope
validation rates run 0.084 to 0.923 — so a (epitope, study) cell is a cohort
with its own donors, protocol and sampling. Requiring ≥40 TCRs with both labels
gives **16 folds covering 1,528 TCRs**, against the current five.

**What this does and does not buy.** Pseudo-epitopes share peptides, so they
are not independent the way real epitopes are, and a configuration exploiting
peptide-specific structure would look better than it is. That confound is
testable rather than assumed: the gain is reported separately for folds held
out *within* an epitope that still has other folds, against folds whose whole
epitope is held out. E36 and E109 both found the structural signal is not
peptide-specific, so the two should agree; if they do not, the extra folds are
not usable.

Run: .venv/bin/python scripts/107_pseudo_epitope_folds.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/84_pseudo_epitope_folds.md"
MIN_N = 40
NULL_REPS = 40
EPS = 1e-9
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    eps_all = sorted(test["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    # ---------- the components, exactly as E104 defined them ----------
    Dm = {}
    for spc in ("default", "small"):
        cc = np.load(ROOT / f"data/interim/sceptr/{spc}.npz", allow_pickle=True)
        WW = cc["vectors"] / np.linalg.norm(cc["vectors"], axis=1, keepdims=True)
        pp = {t: i for i, t in enumerate(cc["tcr_id"])}
        Dm[spc] = {}
        for e, g in test.groupby("epitope_aa"):
            Xx = WW[[pp[t] for t in g["tcr_id"]]]
            Dd = 1.0 - Xx @ Xx.T
            np.fill_diagonal(Dd, np.inf)
            Dm[spc][e] = (Dd, g["tcr_id"].to_numpy())
    Dm["tcrdist"] = {}
    for e, g in test.groupby("epitope_aa"):
        cc = np.load(ROOT / f"data/interim/tcrdist/{e}.npz", allow_pickle=True)
        pp = {t: i for i, t in enumerate(cc["tcr_id"])}
        r = [pp[t] for t in g["tcr_id"]]
        Dd = (cc["alpha"] + cc["beta"]).astype(float)[np.ix_(r, r)]
        np.fill_diagonal(Dd, np.inf)
        Dm["tcrdist"][e] = (Dd, g["tcr_id"].to_numpy())
    nn1 = {s: {e: float(np.median(np.sort(Dm[s][e][0], axis=1)[:, 0]))
               for e in eps_all} for s in Dm}

    def d3(p):
        out = []
        for s in ("default", "small", "tcrdist"):
            gm = float(np.exp(np.mean(np.log(list(nn1[s].values())))))
            sg = {e: float(GRID[s][int(np.argmin(np.abs(
                np.log(GRID[s]) - np.log(max(FIXED[s] * (gm / nn1[s][e]) ** p,
                                             1e-12)))))]) for e in eps_all}
            v = pd.Series(index=test["tcr_id"], dtype=float)
            for e in eps_all:
                Dd, idl = Dm[s][e]
                v.loc[idl] = np.exp(-Dd / sg[e]).mean(axis=1)
            out.append(rank(v))
        return sum(out) / 3

    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    pae = rank(boltz["pae_tcr_pep_min"])
    fz = np.load(ROOT / "data/interim/contact_fingerprints.npz", allow_pickle=True)
    rfp = fz["fingerprint"].astype(float)
    wh = {t: i for i, t in enumerate(fz["tcr_id"])}
    nzf = rfp / (np.linalg.norm(rfp, axis=1, keepdims=True) + 1e-9)
    vv, lf, tn = {}, {}, {}
    for e, g in test.groupby("epitope_aa"):
        idl = [t for t in g["tcr_id"] if t in wh]
        B = nzf[[wh[t] for t in idl]]
        dd = 1.0 - B @ B.T
        np.fill_diagonal(dd, np.inf)
        for i, t in enumerate(idl):
            vv[t] = -float(np.sort(dd[i])[:8].mean())
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        for t, v in zip(g["tcr_id"], (la.map(la.value_counts(normalize=True))
                                      * lb.map(lb.value_counts(normalize=True))
                                      ).to_numpy()):
            lf[t] = v
        idx = np.array([po[t] for t in g["tcr_id"]])
        D = 1.0 - W[idx] @ W[idx].T
        np.fill_diagonal(D, np.inf)
        sd = np.sort(D, axis=1)
        tn.update(dict(zip(g["tcr_id"], 1.0 / np.log(np.maximum(
            np.maximum(sd[:, 1], EPS) / np.maximum(sd[:, 0], EPS), 1 + EPS)))))
    fpd = rank(pd.Series(vv).reindex(test["tcr_id"]))
    lfr, tnr = rank(pd.Series(lf)), rank(-pd.Series(tn))
    Eg = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    pf = Eg["polar_frac"].reindex(test["tcr_id"])
    pfr = rank(-pf.fillna(pf.median()))
    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])
    d3c = {p: d3(p) for p in (0.0, 0.5, 1.0)}
    rest = dr - d3(0.0) / 4 - pae / 4

    def build(p, sw, lw, tw, pw):
        s = rest + d3c[p] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        for w, r in ((lw, lfr), (tw, tnr), (pw, pfr)):
            if w:
                s = (1 - w) * s + w * r
        return s

    # ---------- the pseudo-epitope folds ----------
    lab = labels.dropna(subset=["pos"]).copy()
    lab["pos"] = lab["pos"].astype(int)
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
    print(f"{len(folds)} pseudo-epitope folds, "
          f"{sum(len(v) for v in folds.values())} TCRs")
    y = lab.set_index("tcr_id")["pos"]

    def fold_auc(s, ids):
        yy = y.reindex(ids).to_numpy()
        if len(np.unique(yy)) < 2:
            return np.nan
        return roc_auc_score(yy, s.reindex(ids).to_numpy())

    P, S, L, T, PW = [0., .5, 1.], [0., .15, .3], [0., .05], [0., .05, .10], \
        [0., .05, .10, .15]
    grid = {k: build(*k) for k in product(P, S, L, T, PW)}
    per = {k: {fk: fold_auc(v, ids) for fk, ids in folds.items()}
           for k, v in grid.items()}
    dper = {fk: fold_auc(dr, ids) for fk, ids in folds.items()}

    def nested(fold_keys, lab_y):
        held, picks = [], []
        for fk in fold_keys:
            tr = [x for x in fold_keys if x != fk]
            best = max(per, key=lambda k: np.nanmean([per[k][t] for t in tr]))
            picks.append(best)
            held.append(per[best][fk] - dper[fk])
        return float(np.nanmean(held)), len(set(picks))

    keys = list(folds)
    g16, u16 = nested(keys, y)
    print(f"\n=== nested selection over {len(keys)} pseudo-folds ===")
    print(f"  gain over D: {g16:+.4f}   distinct configurations chosen: {u16}")

    # within-epitope versus whole-epitope holdout, the confound check
    multi = [e for e in eps_all if sum(1 for k in keys if k[0] == e) > 1]
    wi, we = [], []
    for fk in keys:
        tr = [x for x in keys if x != fk]
        b = max(per, key=lambda k: np.nanmean([per[k][t] for t in tr]))
        d = per[b][fk] - dper[fk]
        (wi if fk[0] in multi else we).append(d)
        tr2 = [x for x in keys if x[0] != fk[0]]
        if tr2:
            b2 = max(per, key=lambda k: np.nanmean([per[k][t] for t in tr2]))
            we.append(per[b2][fk] - dper[fk])
    print(f"  fold held out, other folds of its epitope still in training: "
          f"{np.nanmean(wi):+.4f}  (n={len(wi)})")
    print(f"  whole epitope held out:                                     "
          f"{np.nanmean(we):+.4f}  (n={len(we)})")

    # ---------- the null, with 16 folds ----------
    print(f"\n=== permutation null with {len(keys)} folds ({NULL_REPS} reps) ===")
    rng = np.random.default_rng(3)
    nulls = []
    for _ in range(NULL_REPS):
        yp = y.copy()
        for e in eps_all:
            ids = [t for t in y.index if test.set_index('tcr_id').loc[t, 'epitope_aa'] == e] \
                if False else None
        ypv = y.copy()
        for e, g in m.groupby("epitope_aa"):
            ids = [t for t in g["tcr_id"] if t in ypv.index]
            ypv.loc[ids] = rng.permutation(ypv.reindex(ids).to_numpy())

        def fa(s, ids):
            yy = ypv.reindex(ids).to_numpy()
            if len(np.unique(yy)) < 2:
                return np.nan
            return roc_auc_score(yy, s.reindex(ids).to_numpy())
        pn = {k: {fk: fa(v, ids) for fk, ids in folds.items()}
              for k, v in grid.items()}
        dn = {fk: fa(dr, ids) for fk, ids in folds.items()}
        held = []
        for fk in keys:
            tr = [x for x in keys if x != fk]
            b = max(pn, key=lambda k: np.nanmean([pn[k][t] for t in tr]))
            held.append(pn[b][fk] - dn[fk])
        nulls.append(float(np.nanmean(held)))
    nulls = np.array(nulls)
    p = (int((nulls >= g16).sum()) + 1) / (len(nulls) + 1)
    print(f"  null mean {nulls.mean():+.4f}  sd {nulls.std():.4f}  "
          f"max {nulls.max():+.4f}")
    print(f"  observed {g16:+.4f}  ->  permutation p = {p:.3f}")
    print(f"  for comparison, E106 with 5 folds: null sd 0.0059-0.0198, p 0.082-0.200")

    lines = ["# E118 — Sixteen folds instead of five", "",
             "Generated by `scripts/107_pseudo_epitope_folds.py`.", "",
             "More real epitopes are not obtainable: IMMREP25 is a *binding*",
             "task not a *validation* task, and of 21 VDJdb epitopes with ≥80",
             "paired-chain TCRs exactly **one** has enough single- and",
             "multi-reference TCRs for E90's replication proxy. What is",
             "obtainable is **(epitope, study) pseudo-cohorts** — E90 measured",
             "same-epitope validation rates of 0.084 to 0.923 across studies, so",
             "these are cohorts with their own donors, protocol and sampling.", "",
             f"**{len(folds)} folds covering "
             f"{sum(len(v) for v in folds.values())} TCRs**, against five.", "",
             "## Result", "",
             f"- nested selection over {len(keys)} folds: **{g16:+.4f}** over D, "
             f"{u16} distinct configurations chosen",
             f"- fold held out with its epitope's other folds still training: "
             f"**{np.nanmean(wi):+.4f}** (n={len(wi)})",
             f"- whole epitope held out: **{np.nanmean(we):+.4f}** (n={len(we)})",
             "", "## The null, recomputed at this fold count", "",
             f"| | value |", "|---|---:|",
             f"| null mean | {nulls.mean():+.4f} |",
             f"| null sd | **{nulls.std():.4f}** |",
             f"| null max | {nulls.max():+.4f} |",
             f"| observed | {g16:+.4f} |",
             f"| **permutation p** | **{p:.3f}** |", "",
             "E106 with five folds: null sd 0.0059–0.0198, p 0.082–0.200.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
