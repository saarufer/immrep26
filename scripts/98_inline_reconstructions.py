#!/usr/bin/env python
"""The five remaining inline analyses, reconstructed as code.

E106 found three analyses behind E103 and E104 that had run as shell heredocs
and left no retained script, and `scripts/97` recovered two of them. A second
pass over the log found five such analyses in total, all producing numbers that
are quoted in E95, E99, E103 and E104:

  A  E95   zero-parameter rank-average fusion of the three submissions, and
           nested LOEO over the six candidates
  B  E99   mutual-kNN density added on top of the E98 all-three combination
  C  E103  germline conformity as a *replacement* for D's structure slot
           rather than an addition, over eight slot forms
  D  E104  polar_frac on top of the four-component stack, with and without
           charge_product, every parameter nested
  E  E104  raw polar_frac versus the sequence-residualised version inside the
           same stack

Each section prints the figure the log quotes next to the figure recomputed
here, and asserts they agree, so a future reader can see the reproduction
rather than take it on trust. Sections D and E carry the sign-orientation
defect E106 identified; they are reproduced **as published** so the logged
numbers are checkable, with the defect noted inline rather than silently fixed
— E106 already measured that it costs +0.0000 in the point estimate.

Run: .venv/bin/python scripts/98_inline_reconstructions.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/75_inline_reconstructions.md"
EPS = 1e-9
FIXED = {"default": .1, "small": .1, "tcrdist": 20.}
GRID = {"default": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "small": np.array([.02, .03, .05, .07, .1, .15, .2, .3, .5]),
        "tcrdist": np.array([4., 6., 10., 14., 20., 30., 45., 70.])}
# what the log says, for the assertions
EXPECT = {"A_D+F": 0.8582, "A_D+F+frozen": 0.8563, "A_nested": 0.8600,
          "B_all3": 0.8621, "B_nested": 0.8614,
          "C_pae_fp": 0.8610, "C_conf_only": 0.8534, "C_nested": 0.8586,
          "D_four": 0.8628, "D_five": 0.8659, "E_seqresid": 0.8628}
TOL = 0.0002


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def check(tag, got, lines):
    exp = EXPECT[tag]
    ok = abs(got - exp) <= TOL
    msg = f"{'OK ' if ok else 'DIFF'}  {tag:<16} log {exp:.4f}  recomputed {got:.4f}"
    print("   " + msg)
    lines.append(f"| {tag} | {exp:.4f} | {got:.4f} | {'reproduces' if ok else '**differs**'} |")
    return ok


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    eps_all = sorted(test["epitope_aa"].unique())
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    macro = lambda per: float(np.mean([per[e] for e in pub]))   # noqa: E731
    checks = []

    # ================= A: zero-parameter fusion of the three submissions =====
    print("=== A (E95): zero-parameter fusion of the three submissions ===")
    S = {n: pd.read_csv(ROOT / f"submissions/{f}").set_index("tcr_id")["prediction"]
         for n, f in (("frozen", "density3_struct03.csv"),
                      ("D", "d_density3_m30pool_m50_struct.csv"),
                      ("F", "f_density3_m30_m50_struct_nopool.csv"))}
    R = {n: rank(s) for n, s in S.items()}
    candA = {"frozen": R["frozen"], "D": R["D"], "F": R["F"],
             "D+F": (R["D"] + R["F"]) / 2,
             "D+frozen": (R["D"] + R["frozen"]) / 2,
             "D+F+frozen": (R["D"] + R["F"] + R["frozen"]) / 3}
    perA = {n: hyp.per_epitope(labels, s) for n, s in candA.items()}
    for n in candA:
        print(f"  {n:<14}{macro(perA[n]):9.4f}")
    heldA = []
    for e in pub:
        tr = [x for x in pub if x != e]
        b = max(candA, key=lambda k: np.mean([perA[k][t] for t in tr]))
        heldA.append(perA[b][e])
    checks.append(check("A_D+F", macro(perA["D+F"]), checks))
    checks.append(check("A_D+F+frozen", macro(perA["D+F+frozen"]), checks))
    checks.append(check("A_nested", float(np.mean(heldA)), checks))

    # ================= shared machinery for B-E ==============================
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
    dr = rank(S["D"])
    base = macro(hyp.per_epitope(labels, dr))
    rest = dr - d3(0.0) / 4 - pae / 4
    d3c = {p: d3(p) for p in (0.0, 0.5, 1.0)}

    # ================= B: mutual-kNN on top of the all-three =================
    print("\n=== B (E99): mutual-kNN on top of E98's all-three ===")
    all3 = rest + d3c[0.5] / 4 + (0.7 * pae + 0.3 * fpd) / 4
    all3 = 0.95 * all3 + 0.05 * lfr
    checks.append(check("B_all3", macro(hyp.per_epitope(labels, all3)), checks))
    perB = {w: hyp.per_epitope(labels, (1 - w) * all3 + w * mkr)
            for w in (0.0, 0.05, 0.10, 0.15, 0.20)}
    heldB = []
    for e in pub:
        tr = [x for x in pub if x != e]
        b = max(perB, key=lambda w: np.mean([perB[w][t] for t in tr]))
        heldB.append(perB[b][e])
    checks.append(check("B_nested", float(np.mean(heldB)), checks))

    # ================= C: conformity as a structure-slot replacement =========
    print("\n=== C (E103): germline conformity replacing D's structure slot ===")
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
    cf = ((conf["a"] + conf["b"]) / 2)
    cfr = rank(cf.fillna(cf.median()))
    restC = dr - pae / 4
    slots = {"pae (current)": pae, "conformity": cfr, "fingerprint (E83)": fpd,
             "pae .7 + fp .3 (E100)": 0.7 * pae + 0.3 * fpd,
             "pae .7 + conf .3": 0.7 * pae + 0.3 * cfr,
             "pae .5 + conf .5": 0.5 * pae + 0.5 * cfr,
             "pae .5 + fp .25 + conf .25": 0.5 * pae + 0.25 * fpd + 0.25 * cfr,
             "pae .6 + fp .2 + conf .2": 0.6 * pae + 0.2 * fpd + 0.2 * cfr}
    perC = {k: hyp.per_epitope(labels, restC + v / 4) for k, v in slots.items()}
    for k in slots:
        print(f"  {k:<28}{macro(perC[k]):9.4f}{macro(perC[k])-base:+9.4f}")
    heldC = []
    for e in pub:
        tr = [x for x in pub if x != e]
        b = max(perC, key=lambda k: np.mean([perC[k][t] for t in tr]))
        heldC.append(perC[b][e])
    checks.append(check("C_pae_fp", macro(perC["pae .7 + fp .3 (E100)"]), checks))
    checks.append(check("C_conf_only", macro(perC["conformity"]), checks))
    checks.append(check("C_nested", float(np.mean(heldC)), checks))

    # ================= D/E: polar_frac in the stack ==========================
    print("\n=== D/E (E104): polar_frac in the four-component stack ===")
    E_ = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    pf = E_["polar_frac"].reindex(test["tcr_id"])
    pf = pf.fillna(pf.median())
    # NOTE: published orientation uses all five epitopes (the E106 defect).
    # Reproduced as published; E106 measured the cost at +0.0000.
    pfr = rank(-pf)
    NO = {'A': 2, 'R': 5, 'N': 4, 'D': 4, 'C': 2, 'Q': 4, 'E': 4, 'G': 2, 'H': 4,
          'I': 2, 'L': 2, 'K': 3, 'M': 2, 'F': 2, 'P': 2, 'S': 3, 'T': 3,
          'W': 3, 'Y': 3, 'V': 2}
    HV = {'A': 5, 'R': 11, 'N': 8, 'D': 8, 'C': 6, 'Q': 9, 'E': 9, 'G': 4,
          'H': 10, 'I': 8, 'L': 8, 'K': 9, 'M': 8, 'F': 11, 'P': 7, 'S': 6,
          'T': 7, 'W': 14, 'Y': 12, 'V': 7}
    nof = lambda s: (sum(NO.get(ch, 2) for ch in str(s))        # noqa: E731
                     / max(sum(HV.get(ch, 8) for ch in str(s)), 1))
    seq = test.set_index("tcr_id").apply(
        lambda r: .5 * nof(r.cdr3_alpha_aa) + .5 * nof(r.cdr3_beta_aa), axis=1)
    pfres = pd.Series(index=pf.index, dtype=float)
    for e, g in test.groupby("epitope_aa"):
        ids = [t for t in g["tcr_id"] if t in pf.index]
        x = seq.reindex(ids).rank(pct=True).to_numpy()
        y = pf.reindex(ids).rank(pct=True).to_numpy()
        pfres.loc[ids] = y - np.polyval(np.polyfit(x, y, 1), x)

    def build(p, sw, lw, tw, pw, V):
        s = rest + d3c[p] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        for w, r in ((lw, lfr), (tw, tnr), (pw, V)):
            if w:
                s = (1 - w) * s + w * r
        return s

    P, Sg, L, T = [0., .5, 1.], [0., .15, .3], [0., .05], [0., .05, .10]
    for tag, V, PW in (("D_four", pfr, [0.]),
                       ("D_five", pfr, [0., .05, .10, .15, .20]),
                       ("E_seqresid", rank(pfres), [0., .05, .10, .15, .20])):
        gp = {k: hyp.per_epitope(labels, build(*k, V))
              for k in product(P, Sg, L, T, PW)}
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            b = max(gp, key=lambda k: np.mean([gp[k][t] for t in tr]))
            held.append(gp[b][e])
        checks.append(check(tag, float(np.mean(held)), checks))

    agree = [c for c in checks if isinstance(c, bool)]
    print(f"\n{sum(agree)}/{len(agree)} figures reproduce within {TOL}")

    lines = ["# E108 — The five remaining inline analyses, reconstructed", "",
             "Generated by `scripts/98_inline_reconstructions.py`. E106 found",
             "analyses behind E95, E99, E103 and E104 that ran as shell heredocs",
             "and left no code; `scripts/97` recovered two, this recovers the rest.",
             "Sections D and E reproduce the published sign orientation, which",
             "E106 identified as a leak and measured at +0.0000.", "",
             "| figure | log | recomputed | |", "|---|---:|---:|---|"]
    lines += [c for c in checks if isinstance(c, str)]
    lines += ["", f"**{sum(agree)}/{len(agree)} reproduce within {TOL}.**", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT.relative_to(ROOT)}")
    return 0 if all(agree) else 1


if __name__ == "__main__":
    raise SystemExit(main())
