#!/usr/bin/env python
"""What if every component that ever showed a gain is stacked at once?

The stack has been built incrementally — three components in E98, four in
E100, five in E104 — and each step was tested against the step before it. That
leaves the user's question unanswered: several further candidates produced a
positive number at some point and were then dropped when a better one arrived.
Stacked together, do they add?

The honest form of the question is not one number but a curve. Each extra
component brings whatever residual signal it has **and** an extra parameter to
select, and E98 measured selection at −0.0011 for three parameters while E104
measured −0.0016 for five. If the per-component gain is falling while the
per-parameter cost is rising, the curve turns over, and where it turns is the
answer.

Nine candidates, ordered by the residual-on-D diagnostic E103 introduced:

  polar_frac (E104)            0.5586    interface N/O share
  structure fingerprint (E83)  0.5570    contact-fingerprint neighbourhood
  charge_product (E104)        0.5542    interface charge complementarity
  germline conformity (E103)   0.5458    V-gene docking conformity
  lenfreq (E96)                0.5422    CDR3 length-pair frequency
  packing (E104)               0.5494    partner atoms per interface atom
  TwoNN LID (E100)             0.5355    local intrinsic dimension
  cross-epitope LID (E101)     --        own minus other-cohort dimension
  mutual-kNN density (E99)     --        reciprocal-neighbour density

plus the per-epitope bandwidth (E97), which is a reparameterisation rather than
an added term. Every weight is re-selected inside each fold; the curve is the
fully nested macro as components are added in that order.

Run: .venv/bin/python scripts/95_full_stack_curve.py
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

REPORT = ROOT / "reports/72_full_stack_curve.md"
EPS = 1e-9
FAMILY_RADIUS = 0.12
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}
WGRID = [0.0, 0.05, 0.10, 0.15]


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

    # ---------- distances and the bandwidth-parameterised density ----------
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

    # ---------- every candidate term ----------
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
        mk.update(dict(zip(g["tcr_id"],
                           np.where(cnt > 0, -np.where(mu, D, 0).sum(axis=1)
                                    / np.maximum(cnt, 1), -1.0) + .01 * cnt)))
    fpd = rank(pd.Series(vv).reindex(test["tcr_id"]))

    # cross-epitope LID
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
    xl = pd.Series(dict(zip(ids_all, L[np.arange(len(ids_all)), oi]
                            - np.nanmean(oth, axis=1))))
    xl = xl.replace([np.inf, -np.inf], np.nan)

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

    def prep(s, invert_if_needed=True):
        s = pd.Series(s).reindex(test["tcr_id"]).replace([np.inf, -np.inf], np.nan)
        s = s.fillna(s.median())
        if invert_if_needed and hyp.mean_auc(labels, s) < 0.5:
            s = -s
        return rank(s)

    TERMS = [
        ("polar_frac (E104)", prep(E["polar_frac"])),
        ("charge_product (E104)", prep(E["charge_product"])),
        ("germline conformity (E103)", prep((conf["a"] + conf["b"]) / 2)),
        ("lenfreq (E96)", prep(pd.Series(lf))),
        ("packing (E104)", prep(E["packing"])),
        ("TwoNN LID (E100)", prep(-pd.Series(tn))),
        ("cross-epitope LID (E101)", prep(-xl)),
        ("mutual-kNN density (E99)", prep(pd.Series(mk))),
    ]

    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = rank(d_sub)
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))
    rest = dr - d3(0.0) / 4 - pae / 4
    d3c = {p: d3(p) for p in (0.0, 0.5, 1.0)}
    SW = [0.0, 0.15, 0.3]

    def score(p, sw, ws):
        s = rest + d3c[p] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        for w, (_, r) in zip(ws, TERMS):
            if w:
                s = (1 - w) * s + w * r
        return s

    print(f"D {base:.4f}\n")
    print("Each row adds one more candidate. Every weight, plus the bandwidth")
    print("exponent and the structure split, is re-selected by greedy forward")
    print("search inside each fold -- the inner search never sees the held-out")
    print("epitope.\n")
    print(f"  {'k':>2}  {'newest component':<28}{'nested':>9}{'vs D':>9}"
          f"{'oracle':>9}{'cost':>9}{'picks':>7}")
    curve = []
    for k in range(len(TERMS) + 1):
        active = list(range(k))

        def best_config(train):
            """Greedy: bandwidth and structure first, then each weight in turn."""
            p, sw = 0.0, 0.0
            ws = [0.0] * len(TERMS)

            def val(p_, sw_, ws_):
                per = hyp.per_epitope(labels, score(p_, sw_, ws_))
                return float(np.mean([per[t] for t in train]))
            p = max(d3c, key=lambda q: val(q, sw, ws))
            sw = max(SW, key=lambda q: val(p, q, ws))
            for _ in range(2):                       # two sweeps
                for j in active:
                    cur = ws[:]
                    best_w, best_v = ws[j], -1.0
                    for w in WGRID:
                        cur[j] = w
                        v = val(p, sw, cur)
                        if v > best_v:
                            best_v, best_w = v, w
                    ws[j] = best_w
            return p, sw, tuple(ws)

        held, picks = [], []
        for e in pub:
            cfg = best_config([x for x in pub if x != e])
            picks.append(cfg)
            held.append(hyp.per_epitope(labels, score(*cfg[:2], list(cfg[2])))[e])
        nested = float(np.mean(held))
        ocfg = best_config(pub)
        operr = hyp.per_epitope(labels, score(*ocfg[:2], list(ocfg[2])))
        oracle = float(np.mean([operr[t] for t in pub]))
        newest = TERMS[k - 1][0] if k else "(D plus bandwidth/structure only)"
        print(f"  {k:>2}  {newest:<28}{nested:9.4f}{nested-base:+9.4f}"
              f"{oracle:9.4f}{nested-oracle:+9.4f}{len(set(picks)):7d}")
        curve.append([k, newest, nested, nested - base, oracle,
                      nested - oracle, len(set(picks))])

    # ---------- robustness: per-fold detail and a permutation null ----------
    print("\n=== all nine: per-fold detail ===")
    active = list(range(len(TERMS)))

    def best_config_on(train, lab):
        p_, sw_ = 0.0, 0.0
        ws = [0.0] * len(TERMS)

        def val(a, b, c_):
            per = hyp.per_epitope(lab, score(a, b, c_))
            return float(np.mean([per[t] for t in train if t in per]))
        p_ = max(d3c, key=lambda q: val(q, sw_, ws))
        sw_ = max(SW, key=lambda q: val(p_, q, ws))
        for _ in range(2):
            for j in active:
                cur = ws[:]
                bw, bv = ws[j], -1.0
                for w in WGRID:
                    cur[j] = w
                    v = val(p_, sw_, cur)
                    if v > bv:
                        bv, bw = v, w
                ws[j] = bw
        return p_, sw_, tuple(ws)

    print(f"  {'held out':<12}{'D':>9}{'stacked':>10}{'diff':>9}   weights chosen")
    diffs = []
    for e in pub:
        cfg = best_config_on([x for x in pub if x != e], labels)
        v = hyp.per_epitope(labels, score(*cfg[:2], list(cfg[2])))[e]
        d0 = hyp.per_epitope(labels, dr)[e]
        diffs.append(v - d0)
        print(f"  {e:<12}{d0:9.4f}{v:10.4f}{v-d0:+9.4f}   p={cfg[0]:g} sw={cfg[1]:g} "
              + " ".join(f"{w:g}" for w in cfg[2]))
    print(f"  {'MEAN':<12}{'':9}{'':10}{np.mean(diffs):+9.4f}"
          f"   sd across folds {np.std(diffs):.4f}")
    print(f"  folds improved: {sum(d > 0 for d in diffs)}/5")

    print("\n=== permutation null: the same procedure on shuffled labels ===")
    print("  Labels permuted within each epitope, so base rates are preserved")
    print("  and any apparent gain is the procedure's own optimism.")
    rng2 = np.random.default_rng(7)
    nulls = []
    for rep in range(48):
        lp = labels.copy()
        for e in pub:
            m = lp["epitope_aa"] == e
            lp.loc[m, "pos"] = rng2.permutation(lp.loc[m, "pos"].to_numpy())
        hd, bd = [], []
        for e in pub:
            cfg = best_config_on([x for x in pub if x != e], lp)
            per = hyp.per_epitope(lp, score(*cfg[:2], list(cfg[2])))
            hd.append(per[e])
            bd.append(hyp.per_epitope(lp, dr)[e])
        nulls.append(float(np.mean(hd)) - float(np.mean(bd)))
        if (rep + 1) % 8 == 0:
            print(f"    {rep+1:2d} done, running mean {np.mean(nulls):+.4f} "
                  f"sd {np.std(nulls):.4f} max {max(nulls):+.4f}")
    print(f"  null mean {np.mean(nulls):+.4f}  sd {np.std(nulls):.4f}  "
          f"max {max(nulls):+.4f}")
    print(f"  observed {curve[-1][3]:+.4f} -> "
          f"{(sum(n >= curve[-1][3] for n in nulls)+1)/(len(nulls)+1):.3f} "
          f"permutation p")

    bestk = max(curve, key=lambda r: r[2])
    print(f"\n  curve peaks at k={bestk[0]} ({bestk[1]}) "
          f"at {bestk[2]:.4f} ({bestk[3]:+.4f})")
    print(f"  all nine stacked: {curve[-1][2]:.4f} ({curve[-1][3]:+.4f})")

    lines = ["# E105 — Stacking every candidate that ever showed a gain", "",
             "Generated by `scripts/95_full_stack_curve.py`. Components are added",
             "in descending order of their residual-on-D (E103's diagnostic).",
             "Every weight, the bandwidth exponent and the structure split are",
             "re-selected by greedy forward search inside each fold.", "",
             f"D alone **{base:.4f}**.", "",
             "| k | newest component | nested | vs D | oracle | selection cost | distinct picks |",
             "|---:|---|---:|---:|---:|---:|---:|"]
    for r in curve:
        lines.append(f"| {r[0]} | {r[1]} | **{r[2]:.4f}** | {r[3]:+.4f} "
                     f"| {r[4]:.4f} | {r[5]:+.4f} | {r[6]} |")
    lines += ["", "## Robustness", "",
              "| held out | D | stacked | diff |", "|---|---:|---:|---:|"]
    for e, d_ in zip(pub, diffs):
        lines.append(f"| {e} | {hyp.per_epitope(labels, dr)[e]:.4f} "
                     f"| {hyp.per_epitope(labels, dr)[e]+d_:.4f} | {d_:+.4f} |")
    lines += ["", f"Mean {np.mean(diffs):+.4f}, sd across folds {np.std(diffs):.4f}, "
              f"{sum(d > 0 for d in diffs)}/5 folds improved.", "",
              "Permutation null, labels shuffled within epitope, the whole nested",
              f"procedure re-run {len(nulls)} times: mean **{np.mean(nulls):+.4f}**, "
              f"sd {np.std(nulls):.4f}, max {max(nulls):+.4f}. Observed "
              f"{curve[-1][3]:+.4f}, permutation p "
              f"**{(sum(n >= curve[-1][3] for n in nulls)+1)/(len(nulls)+1):.3f}**.", "",
              f"Peak at **k={bestk[0]}**, {bestk[2]:.4f} ({bestk[3]:+.4f}); "
              f"all nine stacked {curve[-1][2]:.4f} ({curve[-1][3]:+.4f}).", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
