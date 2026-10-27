#!/usr/bin/env python
"""Cross-epitope local dimension: is a TCR's lineage structure specific to its
own cohort?

E100 found TwoNN local intrinsic dimension is the first feature here with both
decorrelation from density (ρ −0.337) and margin over chance (0.6461 inverted),
and that adding it drove the four-parameter selection cost to exactly zero.
What it measures within a cohort is whether a receptor sits in lineage-like
structure — neighbours strung along a mutational path rather than scattered.

The obvious next question is whether that structure is **specific**. A genuine
member of an epitope-specific response should look lineage-like in *its own*
cohort and generic everywhere else. A spurious annotation should look much the
same in all of them. E88 built the machinery for exactly this comparison, and
E89's contrastive *density* version of the question failed the typicality gate
on all four variants — but TwoNN passes that gate at ρ −0.084, so the
contrastive form is worth asking again on this functional.

Two symmetry problems have to be handled or the answer is meaningless:

  **Cohort size.** r₁ and r₂ shrink as a cohort grows, so every cohort is
  subsampled to a common size before any LID is computed.

  **Relatives.** A TCR has clonal relatives in its own cohort and none in the
  others, so an unmasked contrast partly measures "does this TCR have relatives
  at all", which is density again. Both versions are computed — family-masked
  and not — and the masked one is the honest measure of specificity while the
  unmasked one is reported to show how much of the effect is relatives.

Everything is label-free and needs only cohort membership, so it is computable
on the private rows. The test that decides it is the same as E100's: does it
add on top of the four-component combination with **all five** parameters
re-selected inside every fold?

Run: .venv/bin/python scripts/91_cross_epitope_lid.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/68_cross_epitope_lid.md"
EPS = 1e-9
FAMILY_RADIUS = 0.12
N_SUB = 110
REPEATS = 25
SEED = 0
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


def twonn_from_rows(d: np.ndarray) -> np.ndarray:
    """TwoNN LID from each row's distances to a reference set (row-wise)."""
    part = np.partition(d, 1, axis=1)[:, :2]
    r1 = np.maximum(part[:, 0], EPS)
    r2 = np.maximum(part[:, 1], EPS)
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

    ids_all = test["tcr_id"].to_numpy()
    ep_all = test["epitope_aa"].to_numpy()
    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}
    X = W[[po[t] for t in ids_all]]
    Dfull = 1.0 - X @ X.T
    np.fill_diagonal(Dfull, 0.0)
    _, fam = connected_components(sp.csr_matrix(Dfull <= FAMILY_RADIUS),
                                  directed=False)
    members = {e: np.where(ep_all == e)[0] for e in epitopes}
    own_idx = np.array([epitopes.index(e) for e in ep_all])

    rng = np.random.default_rng(SEED)
    print("=== computing cross-epitope TwoNN LID "
          f"({REPEATS} subsamples of {N_SUB} per cohort) ===")
    acc = {False: np.zeros((len(ids_all), len(epitopes))),
           True: np.zeros((len(ids_all), len(epitopes)))}
    for _ in range(REPEATS):
        pool = {e: rng.choice(members[e], size=min(N_SUB, len(members[e])),
                              replace=False) for e in epitopes}
        for j, e in enumerate(epitopes):
            ref = pool[e]
            d = Dfull[:, ref].copy()
            d[np.arange(len(ids_all))[:, None] == ref[None, :]] = np.inf
            acc[False][:, j] += twonn_from_rows(d)
            dm = d.copy()
            dm[fam[:, None] == fam[ref][None, :]] = np.inf
            acc[True][:, j] += twonn_from_rows(dm)
    for k in acc:
        acc[k] /= REPEATS

    feats, summary = {}, []
    for masked in (False, True):
        L = acc[masked]
        tag = "family-masked" if masked else "unmasked"
        own = L[np.arange(len(ids_all)), own_idx]
        other = L.copy()
        other[np.arange(len(ids_all)), own_idx] = np.nan
        omean = np.nanmean(other, axis=1)
        omin = np.nanmin(other, axis=1)
        rank_own = (np.nansum(L < own[:, None], axis=1)).astype(float)
        for nm, arr in ((f"LID own ({tag})", own),
                        (f"LID own - mean other ({tag})", own - omean),
                        (f"LID own - min other ({tag})", own - omin),
                        (f"rank of own among 8 ({tag})", rank_own)):
            feats[nm] = pd.Series(dict(zip(ids_all, arr)))
        # r1 == r2 to numerical precision makes TwoNN diverge, so the summary
        # uses medians over finite values; the per-feature AUCs are unaffected
        # because inf is replaced before scoring.
        fin = lambda v: v[np.isfinite(v)]  # noqa: E731
        summary.append([tag, float(np.median(fin(own))),
                        float(np.median(fin(omean))),
                        float(np.median(fin(own - omean))),
                        int((~np.isfinite(own)).sum())])

    print(f"\n  {'version':<16}{'median own':>12}{'median other':>14}"
          f"{'median contrast':>17}{'non-finite':>12}")
    for s in summary:
        print(f"  {s[0]:<16}{s[1]:12.3f}{s[2]:14.3f}{s[3]:+17.3f}{s[4]:12d}")

    # reference quantities
    dens = {}
    for e, g in test.groupby("epitope_aa"):
        idx = np.array([po[t] for t in g["tcr_id"]])
        D = 1.0 - W[idx] @ W[idx].T
        np.fill_diagonal(D, np.inf)
        dens.update(dict(zip(g["tcr_id"], -np.sort(D, axis=1)[:, :8].mean(axis=1))))
    dens = pd.Series(dens)
    twonn_own_plain = {}
    for e, g in test.groupby("epitope_aa"):
        idx = np.array([po[t] for t in g["tcr_id"]])
        D = 1.0 - W[idx] @ W[idx].T
        np.fill_diagonal(D, np.inf)
        twonn_own_plain.update(dict(zip(g["tcr_id"], twonn_from_rows(D))))
    twonn_plain = pd.Series(twonn_own_plain)

    print("\n=== each contrast on task A ===")
    print(f"  {'feature':<38}{'macro':>8}{'inv':>8}{'typ ρ':>8}{'gates':>8}"
          f"{'ρ dens':>9}{'ρ TwoNN':>10}")
    rows = []
    for nm, s in feats.items():
        s = s.replace([np.inf, -np.inf], np.nan)
        s = s.fillna(s.median())
        feats[nm] = s
        g = typicality_gate(test, scored, s)
        pv, _ = private_variance_gate(test, s)
        m = hyp.mean_auc(labels, s)
        rho = g["rho_with_typicality"]
        rd = float(np.mean([spearmanr(s.reindex(gg["tcr_id"]),
                                      dens.reindex(gg["tcr_id"])).statistic
                            for _, gg in test.groupby("epitope_aa")]))
        rt = float(np.mean([spearmanr(s.reindex(gg["tcr_id"]),
                                      twonn_plain.reindex(gg["tcr_id"])).statistic
                            for _, gg in test.groupby("epitope_aa")]))
        ok = "pass" if (abs(rho) < 0.6378 and pv) else "FAIL"
        print(f"  {nm:<38}{m:8.4f}{1-m:8.4f}{rho:+8.3f}{ok:>8}{rd:+9.3f}{rt:+10.3f}")
        rows.append([nm, m, rho, ok, rd, rt])

    # ---------- the decisive test: on top of E100's four, all five nested ----------
    print("\n=== on top of E100's four components, all five parameters nested ===")
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
    for e, g in test.groupby("epitope_aa"):
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        lf.loc[g["tcr_id"]] = (la.map(la.value_counts(normalize=True))
                               * lb.map(lb.value_counts(normalize=True))).to_numpy()
    lfr = rank(lf)
    tnr = rank(-twonn_plain)
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = rank(d_sub)
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))
    rest = dr - d3(0.0) / 4 - pae / 4
    d3c = {p_: d3(p_) for p_ in (0.0, 0.5, 1.0)}

    # The question is whether the CONTRAST adds, so the fifth component is the
    # contrast itself. Picking whichever variant scores best would select
    # "LID own", which is plain TwoNN recomputed on subsamples and is already
    # the fourth component.
    best_name = "LID own - mean other (family-masked)"
    cand = feats[best_name]
    cand_r = rank(cand if hyp.mean_auc(labels, cand) >= 0.5 else -cand)
    print(f"  fifth component: {best_name}")

    def build(p_, sw, lw, tw, xw):
        s = rest + d3c[p_] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        if lw:
            s = (1 - lw) * s + lw * lfr
        if tw:
            s = (1 - tw) * s + tw * tnr
        if xw:
            s = (1 - xw) * s + xw * cand_r
        return s

    P, S, L, T, Xw = [0.0, 0.5, 1.0], [0.0, 0.15, 0.3], [0.0, 0.05], \
        [0.0, 0.05, 0.10], [0.0, 0.05, 0.10]
    gp = {k: hyp.per_epitope(labels, build(*k))
          for k in product(P, S, L, T, Xw)}
    print(f"  {len(gp)} combinations")
    out = {}
    for tag, grid in (("four (E100, xw=0)", {k: v for k, v in gp.items()
                                             if k[4] == 0.0}),
                      ("five", gp)):
        held, picks = [], []
        for e in pub:
            tr = [x for x in pub if x != e]
            bk = max(grid, key=lambda k: np.mean([grid[k][t] for t in tr]))
            picks.append((e, bk))
            held.append(grid[bk][e])
        n = float(np.mean(held))
        ok_ = max(grid, key=lambda k: np.mean([grid[k][t] for t in pub]))
        ov = float(np.mean([grid[ok_][t] for t in pub]))
        uniq = len(set(b for _, b in picks))
        print(f"  {tag:<20} nested {n:.4f} ({n-base:+.4f} vs D)   "
              f"oracle {ov:.4f}   cost {n-ov:+.4f}   distinct picks {uniq}")
        out[tag] = (n, ov, picks, uniq)

    lines = ["# E101 — Cross-epitope local dimension", "",
             "Generated by `scripts/91_cross_epitope_lid.py`.", "",
             "## Own versus other cohorts", "",
             "| version | mean LID own | mean LID other | mean contrast |",
             "|---|---:|---:|---:|"]
    for s in summary:
        lines.append(f"| {s[0]} | {s[1]:.3f} | {s[2]:.3f} | **{s[3]:+.3f}** |")
    lines += ["", "## Each contrast on task A", "",
              "| feature | macro | inverted | ρ typicality | gates | ρ density | ρ TwoNN-own |",
              "|---|---:|---:|---:|---|---:|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {1-r[1]:.4f} | {r[2]:+.3f} "
                     f"| {r[3]} | {r[4]:+.3f} | {r[5]:+.3f} |")
    lines += ["", "## On top of E100's four, all parameters nested", "",
              f"Fifth component: `{best_name}`. D alone {base:.4f}.", "",
              "| configuration | nested | vs D | oracle | selection cost | distinct picks |",
              "|---|---:|---:|---:|---:|---:|"]
    for tag, (n, ov, picks, uniq) in out.items():
        lines.append(f"| {tag} | **{n:.4f}** | {n-base:+.4f} | {ov:.4f} "
                     f"| {n-ov:+.4f} | {uniq} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
