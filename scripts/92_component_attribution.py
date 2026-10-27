#!/usr/bin/env python
"""Where does the +0.0035 live, and will it transfer?

E98 through E101 assembled five components on top of D, reaching +0.0035 fully
nested. The macro number hides the thing that decides whether any of it is real
on the private split: **the gain is not spread evenly**. Every experiment in the
sequence reported most of it on NLVPMVATV and TTDPSFLGRY. If that is a property
of those two cohorts rather than of the method, the private estimate has to come
down; if the contribution profile is consistent across epitopes, it transfers.

Three measurements, none of them a new candidate:

  1. **Exact Shapley attribution per epitope.** Five components give 32 subsets,
     so the Shapley value of each component on each epitope is computable
     exactly rather than approximated. Shapley is the right tool because the
     components interact — E98 measured 98% additivity, which is high but not
     100% — and a simple add-one-in or leave-one-out would mis-attribute the
     remainder.

  2. **Is the profile consistent?** If each component contributes in roughly
     the same proportion everywhere, the macro gain is a property of the
     method. If the profiles are uncorrelated across epitopes, it is five
     separate accidents.

  3. **How much does each component even move the private cohorts?** A
     component that leaves a private ranking untouched contributes exactly zero
     there whatever it does publicly — E97's bandwidth already has this problem
     on two of three. The rank perturbation is measurable on the private rows
     without labels, and calibrating public gain against public perturbation
     gives a projection with a stated basis rather than an assumption.

Run: .venv/bin/python scripts/92_component_attribution.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import combinations
from math import factorial
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/69_component_attribution.md"
EPS = 1e-9
FAMILY_RADIUS = 0.12
N_SUB, REPEATS, SEED = 110, 25, 0
FIXED = {"default": 0.1, "small": 0.1, "tcrdist": 20.0}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}
# the configuration E100/E101 selected
P, SW, LW, TW, XW = 0.5, 0.3, 0.05, 0.05, 0.05
NAMES = ["bandwidth", "structure", "lenfreq", "twonn", "xlid"]


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def twonn_rows(d):
    part = np.partition(d, 1, axis=1)[:, :2]
    r1 = np.maximum(part[:, 0], EPS)
    r2 = np.maximum(part[:, 1], EPS)
    return 1.0 / np.log(np.maximum(r2 / r1, 1.0 + EPS))


def main() -> int:
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    epitopes = sorted(test["epitope_aa"].unique())
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    priv = [e for e in epitopes if e not in pub]
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    ids_all = test["tcr_id"].to_numpy()
    ep_all = test["epitope_aa"].to_numpy()

    # ---------- rebuild every ingredient ----------
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
    nn1 = {s_: {e: float(np.median(np.sort(Dm[s_][e][0], axis=1)[:, 0]))
                for e in epitopes} for s_ in Dm}

    def d3(p):
        out = []
        for s_ in ("default", "small", "tcrdist"):
            gm = float(np.exp(np.mean(np.log(list(nn1[s_].values())))))
            sg = {e: float(GRID[s_][int(np.argmin(np.abs(
                np.log(GRID[s_]) - np.log(max(FIXED[s_] * (gm / nn1[s_][e]) ** p,
                                              1e-12)))))]) for e in epitopes}
            v = pd.Series(index=test["tcr_id"], dtype=float)
            for e in epitopes:
                Dd, idl = Dm[s_][e]
                v.loc[idl] = np.exp(-Dd / sg[e]).mean(axis=1)
            out.append(rank(v))
        return sum(out) / 3

    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    pae = rank(boltz["pae_tcr_pep_min"])
    fz = np.load(ROOT / "data/interim/contact_fingerprints.npz", allow_pickle=True)
    rf = fz["fingerprint"].astype(float)
    wh = {t: i for i, t in enumerate(fz["tcr_id"])}
    nz = rf / (np.linalg.norm(rf, axis=1, keepdims=True) + 1e-9)
    vv = {}
    for e, g in test.groupby("epitope_aa"):
        idl = [t for t in g["tcr_id"] if t in wh]
        B = nz[[wh[t] for t in idl]]
        dd = 1.0 - B @ B.T
        np.fill_diagonal(dd, np.inf)
        for i, t in enumerate(idl):
            vv[t] = -float(np.sort(dd[i])[:8].mean())
    fpd = rank(pd.Series(vv).reindex(test["tcr_id"]))
    lf = pd.Series(index=test["tcr_id"], dtype=float)
    tn = {}
    for e, g in test.groupby("epitope_aa"):
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        lf.loc[g["tcr_id"]] = (la.map(la.value_counts(normalize=True))
                               * lb.map(lb.value_counts(normalize=True))).to_numpy()
        idx = np.array([po[t] for t in g["tcr_id"]])
        Dd = 1.0 - W[idx] @ W[idx].T
        np.fill_diagonal(Dd, np.inf)
        tn.update(dict(zip(g["tcr_id"], twonn_rows(Dd))))
    lfr, tnr = rank(lf), rank(-pd.Series(tn))

    Dfull = 1.0 - W[[po[t] for t in ids_all]] @ W[[po[t] for t in ids_all]].T
    np.fill_diagonal(Dfull, 0.0)
    _, fam = connected_components(sp.csr_matrix(Dfull <= FAMILY_RADIUS),
                                  directed=False)
    members = {e: np.where(ep_all == e)[0] for e in epitopes}
    rng = np.random.default_rng(SEED)
    L = np.zeros((len(ids_all), len(epitopes)))
    for _ in range(REPEATS):
        pool = {e: rng.choice(members[e], size=min(N_SUB, len(members[e])),
                              replace=False) for e in epitopes}
        for j, e in enumerate(epitopes):
            ref = pool[e]
            d = Dfull[:, ref].copy()
            d[np.arange(len(ids_all))[:, None] == ref[None, :]] = np.inf
            d[fam[:, None] == fam[ref][None, :]] = np.inf
            L[:, j] += twonn_rows(d)
    L /= REPEATS
    oi = np.array([epitopes.index(e) for e in ep_all])
    own = L[np.arange(len(ids_all)), oi]
    oth = L.copy()
    oth[np.arange(len(ids_all)), oi] = np.nan
    xl = pd.Series(dict(zip(ids_all, own - np.nanmean(oth, axis=1))))
    xl = xl.replace([np.inf, -np.inf], np.nan)
    xlr = rank(-xl.fillna(xl.median()))

    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = rank(d_sub)
    rest = dr - d3(0.0) / 4 - pae / 4
    d3_0, d3_p = d3(0.0), d3(P)

    def build(on: frozenset) -> pd.Series:
        s = rest + (d3_p if "bandwidth" in on else d3_0) / 4 \
            + ((1 - SW) * pae + SW * fpd if "structure" in on else pae) / 4
        if "lenfreq" in on:
            s = (1 - LW) * s + LW * lfr
        if "twonn" in on:
            s = (1 - TW) * s + TW * tnr
        if "xlid" in on:
            s = (1 - XW) * s + XW * xlr
        return s

    subsets = [frozenset(cc) for k in range(6) for cc in combinations(NAMES, k)]
    val = {ss: hyp.per_epitope(labels, build(ss)) for ss in subsets}
    base = {e: val[frozenset()][e] for e in pub}
    full = {e: val[frozenset(NAMES)][e] for e in pub}

    # ---------- 1. exact Shapley, per epitope ----------
    n = len(NAMES)
    shap = {e: {} for e in pub}
    for e in pub:
        for comp in NAMES:
            tot = 0.0
            for ss in subsets:
                if comp in ss:
                    continue
                w = factorial(len(ss)) * factorial(n - len(ss) - 1) / factorial(n)
                tot += w * (val[ss | {comp}][e] - val[ss][e])
            shap[e][comp] = tot
    print("=== 1. exact Shapley attribution, per epitope (x1e4) ===")
    print(f"  {'epitope':<12}" + "".join(f"{c:>12}" for c in NAMES)
          + f"{'total':>10}{'check':>9}")
    for e in pub:
        tot = sum(shap[e].values())
        print(f"  {e:<12}" + "".join(f"{shap[e][c]*1e4:12.1f}" for c in NAMES)
              + f"{tot*1e4:10.1f}{(full[e]-base[e])*1e4:9.1f}")
        assert abs(tot - (full[e] - base[e])) < 1e-9
    macro_sh = {c: float(np.mean([shap[e][c] for e in pub])) for c in NAMES}
    print(f"  {'MACRO':<12}" + "".join(f"{macro_sh[c]*1e4:12.1f}" for c in NAMES)
          + f"{sum(macro_sh.values())*1e4:10.1f}")

    # ---------- 2. is the profile consistent across epitopes? ----------
    print("\n=== 2. is the contribution profile the same everywhere? ===")
    M = np.array([[shap[e][c] for c in NAMES] for e in pub])
    print(f"  {'':<12}" + "".join(f"{e[:4]:>8}" for e in pub))
    for i, e in enumerate(pub):
        print(f"  {e[:11]:<12}" + "".join(
            f"{spearmanr(M[i], M[j]).statistic:+8.2f}" for j in range(len(pub))))
    off = [spearmanr(M[i], M[j]).statistic
           for i in range(len(pub)) for j in range(i + 1, len(pub))]
    print(f"  mean off-diagonal profile correlation: {np.mean(off):+.3f}")
    print(f"  share of the macro gain on NLV+TTD: "
          f"{sum(full[e]-base[e] for e in ('NLVPMVATV','TTDPSFLGRY'))/sum(full[e]-base[e] for e in pub):.0%}")

    # ---------- 3. does each component move the private cohorts at all? ----------
    print("\n=== 3. rank perturbation each component causes, per cohort ===")
    print("  1 - Spearman(D, D+component) within each cohort; 0 means no effect")
    print(f"  {'epitope':<12}{'split':<9}" + "".join(f"{c:>12}" for c in NAMES))
    sp_of = test.drop_duplicates("epitope_aa").set_index("epitope_aa")["split"]
    pert = {}
    for e in epitopes:
        g = test[test["epitope_aa"] == e]["tcr_id"]
        row = []
        for comp in NAMES:
            s = build(frozenset([comp]))
            row.append(1 - spearmanr(dr.reindex(g), s.reindex(g)).statistic)
        pert[e] = row
        print(f"  {e:<12}{sp_of[e]:<9}" + "".join(f"{v:12.4f}" for v in row))

    # ---------- 4. projection ----------
    print("\n=== 4. projecting to the private split ===")
    print("  gain per unit perturbation, measured on the public five, applied")
    print("  to each private cohort's own measured perturbation.")
    proj = {}
    print(f"  {'component':<12}{'public gain':>12}{'public pert':>13}"
          f"{'gain/pert':>11}" + "".join(f"{e[:4]:>10}" for e in priv))
    for ci, comp in enumerate(NAMES):
        pg = macro_sh[comp]
        pp_ = float(np.mean([pert[e][ci] for e in pub]))
        ratio = pg / pp_ if pp_ > 1e-9 else 0.0
        est = [ratio * pert[e][ci] for e in priv]
        proj[comp] = est
        print(f"  {comp:<12}{pg*1e4:12.1f}{pp_:13.4f}{ratio*1e4:11.1f}"
              + "".join(f"{v*1e4:10.1f}" for v in est))
    tot_priv = float(np.mean([sum(proj[c][i] for c in NAMES)
                              for i in range(len(priv))]))
    print(f"  {'TOTAL':<12}{sum(macro_sh.values())*1e4:12.1f}"
          + " " * 24 + "".join(f"{sum(proj[c][i] for c in NAMES)*1e4:10.1f}"
                               for i in range(len(priv))))
    print(f"\n  projected private macro gain: {tot_priv:+.4f}"
          f"   against the public {sum(macro_sh.values()):+.4f}")

    lines = ["# E102 — Where the +0.0035 lives, and whether it transfers", "",
             "Generated by `scripts/92_component_attribution.py`. Not a new",
             "candidate: an attribution of the five components already assembled.",
             "", "## 1. Exact Shapley attribution, per epitope (×10⁴)", "",
             "Five components give 32 subsets, so the Shapley value is exact",
             "rather than sampled; the per-epitope totals are asserted against the",
             "measured full-minus-empty difference.", "",
             "| epitope | " + " | ".join(NAMES) + " | total |",
             "|---|" + "---:|" * (len(NAMES) + 1)]
    for e in pub:
        lines.append(f"| {e} | " + " | ".join(f"{shap[e][c]*1e4:.1f}"
                                              for c in NAMES)
                     + f" | **{sum(shap[e].values())*1e4:.1f}** |")
    lines.append("| **MACRO** | " + " | ".join(f"**{macro_sh[c]*1e4:.1f}**"
                                               for c in NAMES)
                 + f" | **{sum(macro_sh.values())*1e4:.1f}** |")
    lines += ["", "## 2. Consistency of the profile", "",
              f"Mean off-diagonal Spearman between per-epitope contribution "
              f"profiles: **{np.mean(off):+.3f}**.",
              f"Share of the macro gain on NLVPMVATV + TTDPSFLGRY: "
              f"**{sum(full[e]-base[e] for e in ('NLVPMVATV','TTDPSFLGRY'))/sum(full[e]-base[e] for e in pub):.0%}**.",
              "", "## 3. Rank perturbation per cohort", "",
              "`1 − Spearman(D, D+component)` inside each cohort. Zero means the",
              "component cannot contribute there whatever it does elsewhere.", "",
              "| epitope | split | " + " | ".join(NAMES) + " |",
              "|---|---|" + "---:|" * len(NAMES)]
    for e in epitopes:
        lines.append(f"| {e} | {sp_of[e]} | "
                     + " | ".join(f"{v:.4f}" for v in pert[e]) + " |")
    lines += ["", "## 4. Projection to the private split", "",
              "Gain per unit perturbation measured on the public five, applied to",
              "each private cohort's own measured perturbation.", "",
              "| component | public gain ×10⁴ | public perturbation | ratio | "
              + " | ".join(priv) + " |",
              "|---|---:|---:|---:|" + "---:|" * len(priv)]
    for ci, comp in enumerate(NAMES):
        pp_ = float(np.mean([pert[e][ci] for e in pub]))
        lines.append(f"| {comp} | {macro_sh[comp]*1e4:.1f} | {pp_:.4f} "
                     f"| {(macro_sh[comp]/pp_ if pp_>1e-9 else 0)*1e4:.1f} | "
                     + " | ".join(f"{v*1e4:.1f}" for v in proj[comp]) + " |")
    lines += ["",
              f"**Projected private macro gain {tot_priv:+.4f}** against the "
              f"public {sum(macro_sh.values()):+.4f}.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
