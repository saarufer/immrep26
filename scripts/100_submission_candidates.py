#!/usr/bin/env python
"""Candidate entries for the private split, each betting on a different thing.

The rules permit several teams, so the question is not "which single
configuration is best" but "which small set of configurations covers the ways
the private split could differ from the public one". Five candidates are built,
each making one bet that the others do not, plus the three already submitted
for reference.

What the project established that bears on the choice:

  **E84/E27** the pool term contributes almost nothing privately — CINGVCWTV
  has 11 same-epitope pool rows, LTDEMIAQY 5, LLWNGPMAV 191 — and `D − F`
  is measured at +0.0044 publicly with a dose-response that vanishes below
  ~60 pool rows. D should still beat F privately, but by about a third of its
  public margin.

  **E102** the five-component stack transfers nearly in full: projected private
  +0.0039 against public +0.0042, because the one component that does not reach
  CINGVCWTV and LLWNGPMAV (the bandwidth) is also the smallest.

  **E86** CINGVCWTV and LLWNGPMAV are structural near-twins of TTDPSFLGRY at
  distance 1.40 and 1.47; LTDEMIAQY has no close public analogue.

  **E106** none of the stacked gains is significant — permutation p 0.200 —
  so these are hedges, not improvements, and the honest expectation for every
  candidate below is that it lands within noise of D.

  **E85** per-epitope sampling error on the private cohorts is 0.025–0.060,
  so differences of 0.005 between candidates are far inside the private
  evaluation's own noise, whatever happens publicly.

The frozen entry is never touched; its sha256 is verified before anything is
written.

Run: .venv/bin/python scripts/100_submission_candidates.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/77_submission_candidates.md"
OUT = ROOT / "submissions"
FROZEN_SHA = "523a4692df02467b"
EPS = 1e-9
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}
PRIVATE_BANDS = {"CINGVCWTV": (0.80, 0.83), "LLWNGPMAV": (0.80, 0.83),
                 "LTDEMIAQY": (0.66, 0.68)}
PUBLISHED = {"GILGFVFTL": 0.871, "YLQPRTFLL": 0.928, "NLVPMVATV": 0.809,
             "TTDPSFLGRY": 0.826, "GLCTLVAML": 0.690}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    sha = hashlib.sha256((OUT / "density3_struct03.csv").read_bytes()).hexdigest()
    assert sha.startswith(FROZEN_SHA), f"frozen entry changed: {sha[:16]}"
    print(f"frozen entry intact ({sha[:16]})")

    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    eps_all = sorted(test["epitope_aa"].unique())
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    # ---------- ingredients ----------
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
        sd = np.sort(D, axis=1)
        r1, r2 = np.maximum(sd[:, 0], EPS), np.maximum(sd[:, 1], EPS)
        tn.update(dict(zip(g["tcr_id"], 1.0 / np.log(np.maximum(r2 / r1, 1 + EPS)))))
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
    lfr, tnr, mkr = rank(pd.Series(lf)), rank(-pd.Series(tn)), rank(pd.Series(mk))
    Eg = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    pf = Eg["polar_frac"].reindex(test["tcr_id"])
    pfr = rank(-pf.fillna(pf.median()))
    cq = Eg["charge_product"].reindex(test["tcr_id"])
    cqr = rank(cq.fillna(cq.median()))
    pk = Eg["packing"].reindex(test["tcr_id"])
    pkr = rank(pk.fillna(pk.median()))
    PD = pd.read_csv(ROOT / "data/interim/polar_frac_parts.csv").set_index("tcr_id")
    pg = PD["pf_germline"].reindex(test["tcr_id"])
    pgr = rank(-pg.fillna(pg.median()))

    subs = {n: rank(pd.read_csv(OUT / f).set_index("tcr_id")["prediction"])
            for n, f in (("frozen", "density3_struct03.csv"),
                         ("D", "d_density3_m30pool_m50_struct.csv"),
                         ("F", "f_density3_m30_m50_struct_nopool.csv"))}

    def restructure(base, p, sw):
        """Replace the bandwidth-0 density3 and bare-pae structure slot."""
        return base - d3(0.0) / 4 - pae / 4 + d3(p) / 4 \
            + ((1 - sw) * pae + sw * fpd) / 4

    def mix(s, *pairs):
        for w, r in pairs:
            if w:
                s = (1 - w) * s + w * r
        return s

    # ---------- the candidates ----------
    cands = {}
    # S1: one new component, the most transferable and least selected
    cands["S1_D_plus_polar"] = mix(subs["D"], (0.15, pfr))
    # S2: E104's five-component configuration
    cands["S2_D_stack5"] = mix(restructure(subs["D"], 0.5, 0.3),
                               (0.05, lfr), (0.05, tnr), (0.15, pfr))
    # S3: the same five on the pool-free base
    cands["S3_F_stack5"] = mix(restructure(subs["F"], 0.5, 0.3),
                               (0.05, lfr), (0.05, tnr), (0.15, pfr))
    # S4: the nine-component maximum-public stack
    cands["S4_D_stack9"] = mix(restructure(subs["D"], 0.5, 0.3),
                               (0.15, pfr), (0.15, cqr), (0.10, pgr),
                               (0.05, lfr), (0.05, pkr), (0.05, tnr),
                               (0.10, mkr))
    # S5: weights chosen on TTDPSFLGRY alone, the private cohorts' twin
    best, bv = None, -1.0
    for wp in (0.05, 0.10, 0.15, 0.20):
        for wt in (0.0, 0.05, 0.10):
            for wl in (0.0, 0.05):
                s = mix(restructure(subs["D"], 0.5, 0.3),
                        (wl, lfr), (wt, tnr), (wp, pfr))
                v = hyp.per_epitope(labels, s)["TTDPSFLGRY"]
                if v > bv:
                    bv, best = v, (wl, wt, wp)
    cands["S5_D_ttd_tuned"] = mix(restructure(subs["D"], 0.5, 0.3),
                                  (best[0], lfr), (best[1], tnr), (best[2], pfr))
    print(f"S5 weights chosen on TTDPSFLGRY alone: lenfreq={best[0]} "
          f"twonn={best[1]} polar={best[2]} (TTD {bv:.4f})")

    # ---------- score and project ----------
    print(f"\n{'candidate':<20}{'public':>9}" + "".join(f"{e[:3]:>8}" for e in pub)
          + f"{'E74 private':>13}")
    rows = []
    for nm, s in list(subs.items()) + list(cands.items()):
        per = hyp.per_epitope(labels, s)
        m = float(np.mean([per[e] for e in pub]))
        fr = np.mean([(per[e] - PUBLISHED[e]) / (1 - PUBLISHED[e]) for e in pub])
        prv = float(np.mean([(sum(b) / 2) + fr * (1 - sum(b) / 2)
                             for b in PRIVATE_BANDS.values()]))
        print(f"{nm:<20}{m:9.4f}" + "".join(f"{per[e]:8.3f}" for e in pub)
              + f"{prv:13.4f}")
        rows.append([nm, m, [per[e] for e in pub], prv])

    OUT.mkdir(exist_ok=True)
    print()
    for nm, s in cands.items():
        p = OUT / f"{nm}.csv"
        pd.DataFrame({"tcr_id": test["tcr_id"],
                      "prediction": s.reindex(test["tcr_id"]).to_numpy()}
                     ).to_csv(p, index=False)
        h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        print(f"  wrote {p.name:<24} sha256 {h}")

    assert hashlib.sha256((OUT / "density3_struct03.csv").read_bytes()
                          ).hexdigest().startswith(FROZEN_SHA)
    print("\nfrozen entry still intact")

    lines = ["# E110 — Five candidate entries, each betting on a different thing",
             "", "Generated by `scripts/100_submission_candidates.py`.", "",
             "| candidate | public | " + " | ".join(e[:3] for e in pub)
             + " | E74 private projection |",
             "|---|---:|" + "---:|" * (len(pub) + 1)]
    for r in rows:
        lines.append(f"| {r[0]} | **{r[1]:.4f}** | "
                     + " | ".join(f"{v:.3f}" for v in r[2])
                     + f" | {r[3]:.4f} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
