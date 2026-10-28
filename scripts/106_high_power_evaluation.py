#!/usr/bin/env python
"""Is the +0.0035 real? A higher-powered test than five folds can give.

E106 established the bottleneck: a nested search over five epitopes produces up
to +0.0203 on shuffled labels, so the stacked gain at +0.0035 to +0.0103 cannot
be distinguished from selection optimism, permutation p 0.200. Every experiment
since has run into that wall.

The wall is not about the features. It is that the outer estimate is a mean
over **five numbers**, and five is too few. But five is only the fold count for
*selection*; the sampling error of a **fixed** configuration can be measured
with far more resolution, because there are 1,782 scorable TCRs rather than
five epitopes.

So the two sources of uncertainty are separated and measured apart:

  **Selection optimism** — already measured by E106's label permutation, which
  re-runs the whole search on shuffled labels. Null sd 0.0059–0.0198 depending
  on how honestly the sign step is nested.

  **Sampling error of a fixed configuration** — never measured. Take E104's
  five-component configuration with its weights *frozen*, and resample the
  evaluation set. The paired difference against D has a sampling distribution,
  and its width says whether the gain would survive a different draw of TCRs.

Three resampling schemes, in increasing conservatism:

  1. **TCR bootstrap within epitope** — the naive one, treats TCRs as
     independent.
  2. **Donor bootstrap** — TCRs from one donor share a rearrangement event and
     are not independent (E57 resolved 158 donors on GILGFVFTL, 43 on
     NLVPMVATV). Resampling donors respects the clustering and is the honest
     version where donors are known.
  3. **Family bootstrap** — resample single-linkage clonal families instead,
     which needs no donor metadata and so covers all eight cohorts including
     the private ones.

A fixed configuration carries no selection, so these intervals answer a
different and cleaner question than E106's: *given this configuration, is it
better than D on data like this?*

Run: .venv/bin/python scripts/106_high_power_evaluation.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/83_high_power_evaluation.md"
N_BOOT = 2000
FAMILY_RADIUS = 0.12
EPS = 1e-9
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}
# E104's configuration, frozen
P, SW, LW, TW, PW = 0.5, 0.3, 0.05, 0.05, 0.15


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
    cross = load_module("cross", "scripts/55_cross_donor_persistence.py")
    eps_all = sorted(test["epitope_aa"].unique())
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    # ---------- rebuild the frozen configuration ----------
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
    cfg = dr - d3(0.0) / 4 - pae / 4 + d3(P) / 4 + ((1 - SW) * pae + SW * fpd) / 4
    cfg = (1 - LW) * cfg + LW * lfr
    cfg = (1 - TW) * cfg + TW * tnr
    cfg = (1 - PW) * cfg + PW * pfr

    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    y = scored.set_index("tcr_id")["pos"]
    point = (float(np.mean([hyp.per_epitope(labels, cfg)[e] for e in pub]))
             - float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub])))
    print(f"frozen E104 configuration vs D, on all data: {point:+.4f}\n")

    # ---------- grouping structures ----------
    donors = cross.donor_sets(test)
    ids_all = test["tcr_id"].to_numpy()
    Xa = W[[po[t] for t in ids_all]]
    Df = 1.0 - Xa @ Xa.T
    np.fill_diagonal(Df, 0.0)
    _, fam = connected_components(sp.csr_matrix(Df <= FAMILY_RADIUS),
                                  directed=False)
    fam_of = dict(zip(ids_all, fam))
    don_of = {t: (sorted(donors[t])[0] if donors[t] else f"__solo_{t}")
              for t in ids_all}

    groups = {}
    for e in pub:
        ids = [t for t in test[test["epitope_aa"] == e]["tcr_id"] if t in y.index]
        groups[e] = {
            "tcr": {t: t for t in ids},
            "donor": {t: don_of[t] for t in ids},
            "family": {t: f"f{fam_of[t]}" for t in ids},
        }
    print("=== how many independent units does each scheme give? ===")
    print(f"  {'epitope':<12}{'TCRs':>7}{'donors':>9}{'families':>10}"
          f"{'donor resolved':>16}")
    unit_rows = []
    for e in pub:
        g = groups[e]
        nd = len({v for v in g["donor"].values() if not v.startswith("__solo")})
        res = np.mean([not v.startswith("__solo") for v in g["donor"].values()])
        print(f"  {e:<12}{len(g['tcr']):7d}{len(set(g['donor'].values())):9d}"
              f"{len(set(g['family'].values())):10d}{res:16.0%}")
        unit_rows.append([e, len(g["tcr"]), len(set(g["donor"].values())),
                          len(set(g["family"].values())), float(res), nd])

    # ---------- paired bootstrap ----------
    def boot(scheme, n=N_BOOT, seed=0):
        rng = np.random.default_rng(seed)
        pre = {}
        for e in pub:
            gmap = groups[e][scheme]
            buckets = {}
            for t, k in gmap.items():
                buckets.setdefault(k, []).append(t)
            keys = list(buckets)
            pre[e] = (keys, buckets,
                      y.reindex(list(gmap)).to_numpy(),
                      dr.reindex(list(gmap)).to_numpy(),
                      cfg.reindex(list(gmap)).to_numpy(),
                      {t: i for i, t in enumerate(gmap)})
        diffs = []
        for _ in range(n):
            per = []
            for e in pub:
                keys, buckets, yy, a, b, pos = pre[e]
                pick = rng.choice(len(keys), size=len(keys), replace=True)
                idx = [pos[t] for k in pick for t in buckets[keys[k]]]
                yb = yy[idx]
                if len(np.unique(yb)) < 2:
                    continue
                per.append(roc_auc_score(yb, b[idx]) - roc_auc_score(yb, a[idx]))
            if per:
                diffs.append(float(np.mean(per)))
        return np.array(diffs)

    print(f"\n=== paired bootstrap of the frozen configuration vs D "
          f"({N_BOOT} replicates) ===")
    print(f"  {'resampling unit':<20}{'mean':>9}{'sd':>9}{'2.5%':>9}{'97.5%':>9}"
          f"{'P(gain>0)':>12}")
    boot_rows = []
    for scheme, nm in (("tcr", "TCRs (naive)"), ("donor", "donors"),
                       ("family", "clonal families")):
        d = boot(scheme)
        lo, hi = np.percentile(d, [2.5, 97.5])
        pg = float((d > 0).mean())
        print(f"  {nm:<20}{d.mean():+9.4f}{d.std():9.4f}{lo:+9.4f}{hi:+9.4f}"
              f"{pg:12.1%}")
        boot_rows.append([nm, float(d.mean()), float(d.std()), float(lo),
                          float(hi), pg])

    lines = ["# E117 — Is the gain real? Separating selection optimism from "
             "sampling error", "",
             "Generated by `scripts/106_high_power_evaluation.py`. E106 measured",
             "**selection optimism** by permuting labels and re-running the whole",
             "search. This measures the other half — the **sampling error of a",
             "frozen configuration** — which five folds cannot resolve but 1,782",
             "scorable TCRs can.", "",
             f"E104's five-component configuration with weights frozen scores",
             f"**{point:+.4f}** against D on all data.", "",
             "## Independent units available per scheme", "",
             "| epitope | TCRs | donors | families | donor resolved |",
             "|---|---:|---:|---:|---:|"]
    for r in unit_rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]:.0%} |")
    lines += ["", f"## Paired bootstrap, {N_BOOT} replicates", "",
              "Resampling TCRs treats them as independent; donors and clonal",
              "families respect the clustering and are the honest versions.", "",
              "| resampling unit | mean | sd | 2.5% | 97.5% | P(gain > 0) |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in boot_rows:
        lines.append(f"| {r[0]} | {r[1]:+.4f} | {r[2]:.4f} | {r[3]:+.4f} "
                     f"| {r[4]:+.4f} | **{r[5]:.1%}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
