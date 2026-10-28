#!/usr/bin/env python
"""Code review of E96-E105: the sign-selection leak, quantified.

A review of the exploration arc found one substantive defect repeated in four
scripts. Every candidate term was oriented by

    signed = s if hyp.mean_auc(labels, s) >= 0.5 else -s

which decides the direction of each feature using the labels of **all five**
epitopes, including the one being held out. It is one bit of leakage per
feature, eight features in E105's stack, and it sits outside the nested loop so
the nesting cannot catch it.

Three versions of E105's nine-component stack are run here:

  **leaky**    the published procedure, sign from all five epitopes
  **honest**   sign decided on the four training epitopes only, inside the fold
  **signed**   no orientation step at all; the weight grid is symmetric about
               zero so the inner search picks the direction itself

The third is the right design: it removes the decision rather than relocating
it, at the cost of doubling the inner search space, which the permutation null
is re-measured to price.

Run: .venv/bin/python scripts/96_review_sign_leak.py
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

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/73_review_sign_leak.md"
EPS = 1e-9
FAMILY_RADIUS = 0.12
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}
W_POS = [0.0, 0.05, 0.10, 0.15]
W_SYM = [-0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15]
NULL_REPS = 24


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def twonn(d):
    q = np.partition(d, 1, axis=1)[:, :2]
    return 1.0 / np.log(np.maximum(np.maximum(q[:, 1], EPS)
                                   / np.maximum(q[:, 0], EPS), 1 + EPS))


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    eps_all = sorted(test["epitope_aa"].unique())
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    stack = load_module("stack", "scripts/95_full_stack_curve.py")  # reuse parsing
    del stack  # only imported to assert it stays importable

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
    vv, lf, tn, mk = {}, {}, {}, {}
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
        tn.update(dict(zip(g["tcr_id"], twonn(D))))
        n = len(idx)
        knn = np.argsort(D, axis=1)[:, :8]
        ink = np.zeros((n, n), bool)
        np.put_along_axis(ink, knn, True, axis=1)
        mu = ink & ink.T
        cnt = mu.sum(axis=1)
        # NOTE: the published mutual-kNN term carries an undeclared constant,
        # `+ 0.01 * cnt`, mixing a distance with a count on a hand-picked
        # scale. It is reproduced here so the comparison is like for like, and
        # flagged in the report rather than silently changed.
        mk.update(dict(zip(g["tcr_id"],
                           np.where(cnt > 0, -np.where(mu, D, 0).sum(axis=1)
                                    / np.maximum(cnt, 1), -1.0) + .01 * cnt)))
    fpd = rank(pd.Series(vv).reindex(test["tcr_id"]))

    ids_all = test["tcr_id"].to_numpy()
    ep_all = test["epitope_aa"].to_numpy()
    Xa = W[[po[t] for t in ids_all]]
    Dfull = 1.0 - Xa @ Xa.T
    np.fill_diagonal(Dfull, 0.0)
    _, fam = connected_components(sp.csr_matrix(Dfull <= FAMILY_RADIUS),
                                  directed=False)
    members = {e: np.where(ep_all == e)[0] for e in eps_all}
    rng = np.random.default_rng(0)
    L = np.zeros((len(ids_all), len(eps_all)))
    for _ in range(25):
        pool = {e: rng.choice(members[e], size=min(110, len(members[e])),
                              replace=False) for e in eps_all}
        for j, e in enumerate(eps_all):
            ref = pool[e]
            d = Dfull[:, ref].copy()
            d[np.arange(len(ids_all))[:, None] == ref[None, :]] = np.inf
            d[fam[:, None] == fam[ref][None, :]] = np.inf
            L[:, j] += twonn(d)
    L /= 25
    oi = np.array([eps_all.index(e) for e in ep_all])
    oth = L.copy()
    oth[np.arange(len(ids_all)), oi] = np.nan
    xl = (pd.Series(dict(zip(ids_all, L[np.arange(len(ids_all)), oi]
                             - np.nanmean(oth, axis=1))))
          .replace([np.inf, -np.inf], np.nan))

    E = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    prof = pd.read_csv(ROOT / "data/interim/germline_docking.csv").set_index("tcr_id")
    gcols = {"a": [x for x in prof.columns if x.startswith(("a_A1b", "a_A2b"))],
             "b": [x for x in prof.columns if x.startswith(("b_A1b", "b_A2b"))]}
    conf = {}
    for tag, gc in (("a", "TRAV_IMGT"), ("b", "TRBV_IMGT")):
        M = prof[gcols[tag]].to_numpy(float)
        M = M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-9)
        ix = {t: i for i, t in enumerate(prof.index)}
        gs = test.set_index("tcr_id")[gc].reindex(prof.index)
        out = np.full(len(prof), np.nan)
        for gn, sub in gs.groupby(gs):
            ii = np.array([ix[t] for t in sub.index])
            if len(ii) < 5:
                continue
            for k in ii:
                o = ii[ii != k]
                md = M[o].mean(0)
                md /= np.linalg.norm(md) + 1e-9
                out[k] = float(M[k] @ md)
        conf[tag] = pd.Series(out, index=prof.index)

    def raw(s):
        s = pd.Series(s).reindex(test["tcr_id"]).replace([np.inf, -np.inf], np.nan)
        return rank(s.fillna(s.median()))

    RAW = [("polar_frac", raw(E["polar_frac"])),
           ("charge_product", raw(E["charge_product"])),
           ("germline conformity", raw((conf["a"] + conf["b"]) / 2)),
           ("lenfreq", raw(pd.Series(lf))),
           ("packing", raw(E["packing"])),
           ("TwoNN LID", raw(pd.Series(tn))),
           ("cross-epitope LID", raw(xl)),
           ("mutual-kNN density", raw(pd.Series(mk)))]

    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))
    rest = dr - d3(0.0) / 4 - pae / 4
    d3c = {p: d3(p) for p in (0.0, 0.5, 1.0)}
    SW = [0.0, 0.15, 0.3]

    def score(p, sw, ws, signs):
        s = rest + d3c[p] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        for w, sg, (_, r) in zip(ws, signs, RAW):
            if w:
                s = (1 - abs(w)) * s + abs(w) * (r if sg > 0 else (1 - r))
        return s

    def search(train, lab, wgrid, signs_fixed):
        p_, sw_ = 0.0, 0.0
        ws = [0.0] * len(RAW)
        signs = list(signs_fixed) if signs_fixed is not None else [1] * len(RAW)

        def val(a, b, c_, sg):
            per = hyp.per_epitope(lab, score(a, b, c_, sg))
            return float(np.mean([per[t] for t in train if t in per]))
        p_ = max(d3c, key=lambda q: val(q, sw_, ws, signs))
        sw_ = max(SW, key=lambda q: val(p_, q, ws, signs))
        for _ in range(2):
            for j in range(len(RAW)):
                best, bv = (ws[j], signs[j]), -1.0
                for w in wgrid:
                    cur = ws[:]
                    cur[j] = abs(w)
                    sg = signs[:]
                    sg[j] = 1 if w >= 0 else -1
                    v = val(p_, sw_, cur, sg)
                    if v > bv:
                        bv, best = v, (abs(w), sg[j])
                ws[j], signs[j] = best
        return p_, sw_, tuple(ws), tuple(signs)

    all_signs = tuple(1 if hyp.mean_auc(labels, pd.Series(
        dict(zip(test["tcr_id"], r.reindex(test["tcr_id"]))))) >= 0.5 else -1
        for _, r in RAW)
    print(f"D {base:.4f}")
    print(f"signs from all five epitopes (the leak): {all_signs}\n")

    def run(mode, lab):
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            if mode == "leaky":
                cfg = search(tr, lab, W_POS, all_signs)
            elif mode == "honest":
                sg = tuple(1 if float(np.mean([hyp.per_epitope(lab, r)[t]
                                               for t in tr])) >= 0.5 else -1
                           for _, r in RAW)
                cfg = search(tr, lab, W_POS, sg)
            else:
                cfg = search(tr, lab, W_SYM, None)
            held.append(hyp.per_epitope(lab, score(*cfg))[e])
        return float(np.mean(held))

    print(f"  {'version':<34}{'nested':>9}{'vs D':>9}")
    res = {}
    for mode, desc in (("leaky", "leaky: sign from all five (published)"),
                       ("honest", "honest: sign from the four training"),
                       ("signed", "signed: symmetric weights, no sign step")):
        v = run(mode, labels)
        res[mode] = v
        print(f"  {desc:<34}{v:9.4f}{v-base:+9.4f}")

    print(f"\n  cost of removing the leak: "
          f"{res['honest']-res['leaky']:+.4f} (honest) "
          f"{res['signed']-res['leaky']:+.4f} (signed)")

    print(f"\n=== permutation null for the honest versions ({NULL_REPS} reps) ===")
    rng2 = np.random.default_rng(11)
    nulls = {"honest": [], "signed": []}
    for rep in range(NULL_REPS):
        lp = labels.copy()
        for e in pub:
            m = lp["epitope_aa"] == e
            lp.loc[m, "pos"] = rng2.permutation(lp.loc[m, "pos"].to_numpy())
        b = float(np.mean([hyp.per_epitope(lp, dr)[e] for e in pub]))
        for mode in ("honest", "signed"):
            nulls[mode].append(run(mode, lp) - b)
        if (rep + 1) % 6 == 0:
            print(f"    {rep+1:2d} done   honest sd {np.std(nulls['honest']):.4f} "
                  f"max {max(nulls['honest']):+.4f}   signed sd "
                  f"{np.std(nulls['signed']):.4f} max {max(nulls['signed']):+.4f}")
    lines = ["# E106 — Code review of E96–E105: the sign-selection leak", "",
             "Generated by `scripts/96_review_sign_leak.py`.", "",
             "Four scripts oriented each candidate with",
             "`signed = s if mean_auc(labels, s) >= 0.5 else -s`, which uses the",
             "labels of all five epitopes including the held-out one. One bit per",
             "feature, eight features, outside the nested loop.", "",
             "| version | nested | vs D |", "|---|---:|---:|"]
    for mode, desc in (("leaky", "leaky (published)"),
                       ("honest", "sign from training folds only"),
                       ("signed", "symmetric weights, no sign step")):
        lines.append(f"| {desc} | {res[mode]:.4f} | **{res[mode]-base:+.4f}** |")
    lines += ["", "## Permutation null, re-measured", "",
              "| version | null mean | null sd | null max | observed | p |",
              "|---|---:|---:|---:|---:|---:|"]
    for mode in ("honest", "signed"):
        n = np.array(nulls[mode])
        obs = res[mode] - base
        p = (int((n >= obs).sum()) + 1) / (len(n) + 1)
        print(f"  {mode:<8} null mean {n.mean():+.4f} sd {n.std():.4f} "
              f"max {n.max():+.4f}   observed {obs:+.4f}   p {p:.3f}")
        lines.append(f"| {mode} | {n.mean():+.4f} | {n.std():.4f} | {n.max():+.4f} "
                     f"| {obs:+.4f} | **{p:.3f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
