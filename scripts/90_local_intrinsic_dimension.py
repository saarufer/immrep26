#!/usr/bin/env python
"""Local intrinsic dimension: the exponent, not the scale.

E99 closed hubness by measuring that this space is not hub-afflicted — density
and in-degree agree at ρ +0.66–0.85, so the in-degree view is not a separate
quantity. That leaves one more functional of the same distance matrix that
genuinely is separate.

Density asks **how far** a TCR's neighbours are. Local intrinsic dimension asks
**how fast the neighbour count grows with radius** — the exponent rather than
the scale — and it is computed from *ratios* of log-distances, so it is
invariant to the overall scale that density measures. Two TCRs can sit at
identical density with very different LID.

The biology it would correspond to is concrete. A receptor inside a clonal
lineage has neighbours strung along a mutational path: few effective degrees of
freedom, **low LID**. A receptor in a generic region of sequence space has
neighbours scattered in every direction: **high LID**. E29 measured that family
recurrence carries part of the ranking (masking each TCR's own family costs
−0.0282), and LID is a continuous, label-free way of asking whether a TCR sits
in lineage-like structure without needing a family radius at all.

Two standard estimators, both within-cohort:

  **Levina–Bickel MLE**   LID_k(x) = −[ (1/(k−1)) Σ_{i<k} log(r_i/r_k) ]⁻¹
  **TwoNN (Facco)**       from the ratio μ = r₂/r₁ of the two nearest distances

The standing gates apply, plus the +0.02 bar, nested LOEO against E15's
−0.0048, and — because E98 already combined three things — the test that
actually matters: does it add **on top of the all-three combination**, with its
weight nested? E99's mutual-kNN density passed every earlier screen and failed
that one at −0.0007.

Run: .venv/bin/python scripts/90_local_intrinsic_dimension.py
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

REPORT = ROOT / "reports/67_local_intrinsic_dimension.md"
EPS = 1e-9
KS = (5, 10, 20, 40)
FIXED = {"default": 0.1, "small": 0.1, "tcrdist": 20.0}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def lid_mle(sorted_d: np.ndarray, k: int) -> np.ndarray:
    """Levina-Bickel maximum-likelihood LID from each row's sorted distances.

    Distances are floored at EPS: one GLCTLVAML pair is exactly coincident and
    log(0) would poison that row. Flooring biases that single row towards a
    high LID, which is the conservative direction for a feature meant to flag
    lineage structure.
    """
    r = np.maximum(sorted_d[:, :k], EPS)
    rk = r[:, [k - 1]]
    lg = np.log(r[:, :k - 1] / rk)
    mean = -lg.mean(axis=1)
    return np.where(mean > EPS, 1.0 / mean, np.nan)


def lid_twonn(sorted_d: np.ndarray) -> np.ndarray:
    """Facco's TwoNN, per point: the local exponent from mu = r2/r1."""
    r1 = np.maximum(sorted_d[:, 0], EPS)
    r2 = np.maximum(sorted_d[:, 1], EPS)
    mu = np.maximum(r2 / r1, 1.0 + EPS)
    return 1.0 / np.log(mu)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    epitopes = sorted(test["epitope_aa"].unique())
    pub = sorted(scored["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}

    feats: dict[str, dict] = {}
    summary = []
    for e, g in test.groupby("epitope_aa"):
        ids = g["tcr_id"].to_numpy()
        X = W[[po[t] for t in ids]]
        D = 1.0 - X @ X.T
        np.fill_diagonal(D, np.inf)
        sd = np.sort(D, axis=1)
        dens = -sd[:, :8].mean(axis=1)
        vals = {"density (reference)": dens,
                "LID TwoNN": lid_twonn(sd)}
        for k in KS:
            vals[f"LID MLE k={k}"] = lid_mle(sd, k)
        for nm, arr in vals.items():
            feats.setdefault(nm, {}).update(dict(zip(ids, arr)))
        l20 = vals["LID MLE k=20"]
        ok = np.isfinite(l20)
        summary.append([e, g["split"].mode().iloc[0], len(ids),
                        float(np.nanmedian(l20)),
                        float(np.nanpercentile(l20, 10)),
                        float(np.nanpercentile(l20, 90)),
                        float(spearmanr(dens[ok], l20[ok]).statistic)])

    print("=== 1. the local dimension of each cohort, and is it density again? ===")
    print(f"  {'epitope':<12}{'split':<9}{'n':>5}{'LID median':>12}"
          f"{'p10':>8}{'p90':>8}{'rho(dens,LID)':>15}")
    for s in summary:
        print(f"  {s[0]:<12}{s[1]:<9}{s[2]:5d}{s[3]:12.2f}{s[4]:8.2f}"
              f"{s[5]:8.2f}{s[6]:+15.3f}")

    series = {k: pd.Series(v) for k, v in feats.items()}
    print("\n=== 2. each estimator on task A ===")
    print(f"  {'feature':<22}{'macro':>8}{'inverted':>10}{'rho w/ typ':>12}"
          f"{'typ':>6}{'priv':>7}{'rho w/ dens':>13}")
    ref = series["density (reference)"]
    rows = []
    for nm, s in series.items():
        s = s.replace([np.inf, -np.inf], np.nan)
        s = s.fillna(s.median())
        series[nm] = s
        g = typicality_gate(test, scored, s)
        pv, _ = private_variance_gate(test, s)
        m = hyp.mean_auc(labels, s)
        rho = g["rho_with_typicality"]
        rd = float(np.mean([spearmanr(s.reindex(gg["tcr_id"]),
                                      ref.reindex(gg["tcr_id"])).statistic
                            for _, gg in test.groupby("epitope_aa")]))
        print(f"  {nm:<22}{m:8.4f}{1-m:10.4f}{rho:+12.3f}"
              f"{'pass' if abs(rho) < 0.6378 else 'FAIL':>6}"
              f"{'pass' if pv else 'FAIL':>7}{rd:+13.3f}")
        rows.append([nm, m, rho, pv, rd])

    # orient every LID so that higher = better, then blend
    print("\n=== 3. blended into D (each LID oriented by its own sign) ===")
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = rank(d_sub)
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))
    print(f"  D alone {base:.4f}")
    WS = [0.05, 0.10, 0.15, 0.20]
    blend = []
    for nm, s in series.items():
        if nm == "density (reference)":
            continue
        signed = s if hyp.mean_auc(labels, s) >= 0.5 else -s
        sr = rank(signed)
        per_w = {0.0: hyp.per_epitope(labels, dr)}
        for w in WS:
            per_w[w] = hyp.per_epitope(labels, (1 - w) * dr + w * sr)
        vals = {w: float(np.mean([per_w[w][e] for e in pub])) for w in per_w}
        bw = max(WS, key=lambda w: vals[w])
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            pick = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
            held.append(per_w[pick][e])
        nested = float(np.mean(held))
        print(f"  {nm:<22} oracle w={bw:.2f} {vals[bw]:.4f} ({vals[bw]-base:+.4f})"
              f"   nested {nested:.4f} ({nested-base:+.4f})")
        blend.append([nm, bw, vals[bw], vals[bw] - base, nested, nested - base])

    # ---------- the test that matters: on top of E98's all-three ----------
    print("\n=== 4. on top of E98's all-three combination ===")
    Dm = {}
    for sp in ("default", "small"):
        cc = np.load(ROOT / f"data/interim/sceptr/{sp}.npz", allow_pickle=True)
        WW = cc["vectors"] / np.linalg.norm(cc["vectors"], axis=1, keepdims=True)
        pp = {t: i for i, t in enumerate(cc["tcr_id"])}
        Dm[sp] = {}
        for e, g in test.groupby("epitope_aa"):
            Xx = WW[[pp[t] for t in g["tcr_id"]]]
            Dd = 1.0 - Xx @ Xx.T
            np.fill_diagonal(Dd, np.inf)
            Dm[sp][e] = (Dd, g["tcr_id"].to_numpy())
    Dm["tcrdist"] = {}
    for e, g in test.groupby("epitope_aa"):
        cc = np.load(ROOT / f"data/interim/tcrdist/{e}.npz", allow_pickle=True)
        pp = {t: i for i, t in enumerate(cc["tcr_id"])}
        r = [pp[t] for t in g["tcr_id"]]
        Dd = (cc["alpha"] + cc["beta"]).astype(float)[np.ix_(r, r)]
        np.fill_diagonal(Dd, np.inf)
        Dm["tcrdist"][e] = (Dd, g["tcr_id"].to_numpy())
    nn1 = {sp: {e: float(np.median(np.sort(Dm[sp][e][0], axis=1)[:, 0]))
                for e in epitopes} for sp in Dm}

    def d3(p):
        out = []
        for sp in ("default", "small", "tcrdist"):
            gm = float(np.exp(np.mean(np.log(list(nn1[sp].values())))))
            sg = {e: float(GRID[sp][int(np.argmin(np.abs(
                np.log(GRID[sp]) - np.log(max(FIXED[sp] * (gm / nn1[sp][e]) ** p,
                                              1e-12)))))]) for e in epitopes}
            v = pd.Series(index=test["tcr_id"], dtype=float)
            for e in epitopes:
                Dd, ids = Dm[sp][e]
                v.loc[ids] = np.exp(-Dd / sg[e]).mean(axis=1)
            out.append(rank(v))
        return sum(out) / 3

    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    pae = rank(boltz["pae_tcr_pep_min"])
    fpz = np.load(ROOT / "data/interim/contact_fingerprints.npz", allow_pickle=True)
    rawf = fpz["fingerprint"].astype(float)
    wh = {t: i for i, t in enumerate(fpz["tcr_id"])}
    nmz = rawf / (np.linalg.norm(rawf, axis=1, keepdims=True) + 1e-9)
    vv = {}
    for e, g in test.groupby("epitope_aa"):
        idl = [t for t in g["tcr_id"] if t in wh]
        B = nmz[[wh[t] for t in idl]]
        dd = 1.0 - B @ B.T
        np.fill_diagonal(dd, np.inf)
        for i, t in enumerate(idl):
            vv[t] = -float(np.sort(dd[i])[:8].mean())
    fpd = rank(pd.Series(vv).reindex(test["tcr_id"]))
    lf = pd.Series(index=test["tcr_id"], dtype=float)
    for e, g in test.groupby("epitope_aa"):
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        lf.loc[g["tcr_id"]] = (la.map(la.value_counts(normalize=True))
                               * lb.map(lb.value_counts(normalize=True))).to_numpy()
    all3 = (dr - d3(0.0) / 4 - pae / 4) + d3(0.5) / 4 + (0.7 * pae + 0.3 * fpd) / 4
    all3 = 0.95 * all3 + 0.05 * rank(lf)
    a3 = float(np.mean([hyp.per_epitope(labels, all3)[e] for e in pub]))
    print(f"  all-three {a3:.4f} ({a3-base:+.4f} vs D)")
    top = []
    for nm, s in series.items():
        if nm == "density (reference)":
            continue
        signed = s if hyp.mean_auc(labels, s) >= 0.5 else -s
        sr = rank(signed)
        per_w = {0.0: hyp.per_epitope(labels, all3)}
        for w in WS:
            per_w[w] = hyp.per_epitope(labels, (1 - w) * all3 + w * sr)
        vals = {w: float(np.mean([per_w[w][e] for e in pub])) for w in per_w}
        bw = max(WS, key=lambda w: vals[w])
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            pick = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
            held.append(per_w[pick][e])
        nested = float(np.mean(held))
        print(f"  + {nm:<20} oracle w={bw:.2f} {vals[bw]:.4f} "
              f"({vals[bw]-a3:+.4f})   nested {nested:.4f} ({nested-a3:+.4f})")
        top.append([nm, bw, vals[bw], vals[bw] - a3, nested, nested - a3])

    # ---------- 5. fully nested over all four parameters ----------
    # Section 4 holds the other three at values chosen on these same epitopes.
    # E98 showed that costs -0.0011 for three parameters; this re-selects all
    # four inside every fold, which is the only number worth quoting.
    print("\n=== 5. fully nested: all four parameters re-chosen inside each fold ===")
    from itertools import product as iproduct
    P_GRID, S_GRID, L_GRID = [0.0, 0.5, 1.0], [0.0, 0.15, 0.3], [0.0, 0.05]
    T_GRID = [0.0, 0.05, 0.10]
    twonn = series["LID TwoNN"]
    twonn_r = rank(twonn if hyp.mean_auc(labels, twonn) >= 0.5 else -twonn)
    d3c = {p_: d3(p_) for p_ in P_GRID}
    rest4 = dr - d3(0.0) / 4 - pae / 4
    gp = {}
    for p_, sw, lw, tw in iproduct(P_GRID, S_GRID, L_GRID, T_GRID):
        sc = rest4 + d3c[p_] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        if lw:
            sc = (1 - lw) * sc + lw * rank(lf)
        if tw:
            sc = (1 - tw) * sc + tw * twonn_r
        gp[(p_, sw, lw, tw)] = hyp.per_epitope(labels, sc)
    print(f"  {len(gp)} combinations "
          f"({len(P_GRID)}x{len(S_GRID)}x{len(L_GRID)}x{len(T_GRID)})")
    held4, picks4 = [], []
    for e in pub:
        tr = [x for x in pub if x != e]
        bk = max(gp, key=lambda k: np.mean([gp[k][t] for t in tr]))
        picks4.append((e, bk))
        held4.append(gp[bk][e])
        print(f"    {e:<12} p={bk[0]:g} struct={bk[1]:g} len={bk[2]:g} "
              f"twonn={bk[3]:g}  ->  {gp[bk][e]:.4f}")
    n4 = float(np.mean(held4))
    ok = max(gp, key=lambda k: np.mean([gp[k][t] for t in pub]))
    ov = float(np.mean([gp[ok][t] for t in pub]))
    print(f"  fully nested   {n4:.4f}  ({n4-base:+.4f} vs D)")
    print(f"  grid oracle    {ov:.4f}  ({ov-base:+.4f} vs D)  at {ok}")
    print(f"  selection cost of four parameters: {n4-ov:+.4f}")
    # The grid here is coarser than E98's (3x3x2 vs 5x4x3), and a coarser grid
    # has lower selection variance on its own. So the three-parameter version is
    # re-run on THIS grid, otherwise the comparison credits TwoNN with a
    # stability that the grid supplied.
    gp3 = {k: v for k, v in gp.items() if k[3] == 0.0}
    held3 = []
    picks3 = []
    for e in pub:
        tr = [x for x in pub if x != e]
        bk = max(gp3, key=lambda k: np.mean([gp3[k][t] for t in tr]))
        picks3.append((e, bk))
        held3.append(gp3[bk][e])
    n3 = float(np.mean(held3))
    ok3 = max(gp3, key=lambda k: np.mean([gp3[k][t] for t in pub]))
    ov3 = float(np.mean([gp3[ok3][t] for t in pub]))
    print(f"\n  --- the fair comparison, three parameters on THIS grid ---")
    print(f"  three-parameter fully nested  {n3:.4f}  ({n3-base:+.4f} vs D)"
          f"   selection cost {n3-ov3:+.4f}")
    print(f"  picks: " + "  ".join(f"{e[:3]}=({b[0]:g},{b[1]:g},{b[2]:g})"
                                   for e, b in picks3))
    print(f"  so TwoNN's own contribution is {n4-n3:+.4f}")
    print(f"  (E98's finer 5x4x3 grid gave +0.0009 with cost -0.0011)")

    lines = ["# E100 — Local intrinsic dimension: the exponent, not the scale", "",
             "Generated by `scripts/90_local_intrinsic_dimension.py`.", "",
             "## 1. Cohort local dimension, and is it density again?", "",
             "| epitope | split | n | LID median | p10 | p90 | ρ(density, LID) |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for s in summary:
        lines.append(f"| {s[0]} | {s[1]} | {s[2]} | {s[3]:.2f} | {s[4]:.2f} "
                     f"| {s[5]:.2f} | **{s[6]:+.3f}** |")
    lines += ["", "## 2. Each estimator on task A", "",
              "| feature | macro | inverted | ρ with typicality | typ gate | priv gate | ρ with density |",
              "|---|---:|---:|---:|---|---|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {1-r[1]:.4f} | {r[2]:+.3f} "
                     f"| {'pass' if abs(r[2])<0.6378 else 'FAIL'} "
                     f"| {'pass' if r[3] else 'FAIL'} | {r[4]:+.3f} |")
    lines += ["", "## 3. Blended into D", "", f"D alone **{base:.4f}**.", "",
              "| feature | oracle w | oracle | vs D | nested | vs D |",
              "|---|---:|---:|---:|---:|---:|"]
    for b in blend:
        lines.append(f"| {b[0]} | {b[1]:.2f} | {b[2]:.4f} | {b[3]:+.4f} "
                     f"| {b[4]:.4f} | **{b[5]:+.4f}** |")
    lines += ["", "## 4. On top of E98's all-three combination", "",
              f"all-three **{a3:.4f}** ({a3-base:+.4f} vs D).", "",
              "| feature | oracle w | oracle | vs all-three | nested | vs all-three |",
              "|---|---:|---:|---:|---:|---:|"]
    for b in top:
        lines.append(f"| {b[0]} | {b[1]:.2f} | {b[2]:.4f} | {b[3]:+.4f} "
                     f"| {b[4]:.4f} | **{b[5]:+.4f}** |")
    lines += ["", "## 5. Fully nested over all four parameters", "",
              f"{len(gp)}-point grid, every parameter re-selected inside each fold.",
              "", "| held out | p | struct_fp | lenfreq | twonn | AUC |",
              "|---|---:|---:|---:|---:|---:|"]
    for e, bk in picks4:
        lines.append(f"| {e} | {bk[0]:g} | {bk[1]:g} | {bk[2]:g} | {bk[3]:g} "
                     f"| {gp[bk][e]:.4f} |")
    lines += ["",
              f"- fully nested: **{n4:.4f}** ({n4-base:+.4f} vs D)",
              f"- grid oracle: {ov:.4f} ({ov-base:+.4f} vs D)",
              f"- selection cost of four parameters: **{n4-ov:+.4f}**",
              "", "Because a coarser grid has lower selection variance on its own,",
              "the three-parameter version is re-run on **this** grid rather than",
              "compared against E98's finer one:", "",
              f"- three parameters, fully nested on this grid: **{n3:.4f}** "
              f"({n3-base:+.4f} vs D), selection cost {n3-ov3:+.4f}",
              f"- **TwoNN's own contribution: {n4-n3:+.4f}**",
              "- E98's finer 5x4x3 grid gave +0.0009 at a cost of -0.0011", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
