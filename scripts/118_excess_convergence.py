#!/usr/bin/env python
"""E135 -- Excess convergence: not how many neighbours, but how unusual they are.

A collaborator's proposal, and the sharpest one yet, because it changes the
measured object rather than adding a descriptor. It also follows directly from
E127, which produced a genuinely odd pair of facts:

  for historical epitope ATTRIBUTION, the positionally-shuffled library C beat
  the real library A (+0.5459 against +0.5284, A-C = -0.0175 over five seeds)
  for functional REVALIDATION, D leans heavily on real-sequence local clustering

So two statistics may be tangled inside one distance density:

  low-order compositional typicality   does t match this epitope's V/J, length
                                       and per-position amino-acid distribution?
  high-order convergence evidence      does t sit in a real sequence family that
                                       the low-order composition cannot explain?

The new quantity, using E127's validated epitope-conditioned null:

    U(t,e) = log( (D_real(t,e) + eps) / (E[D_shuffle(t,e)] + eps) )

**Why this is not E78.** E78 subtracted a V/J + length + Pgen *generation*
background and failed. The null here is the epitope's own observed positional
composition, held exactly (E127 asserts the sorted-multiset equality), with V
and J preserved per member. The question is therefore different: once the
observed per-position composition is fixed, does the joint structure among real
CDR3 residues still predict functional revalidation? It remains a member of the
background-correction family and E78's negative result is the thing to beat.

**Discipline, as proposed.** Nothing is trained and nothing is fused on the
first pass. Three quantities are measured -- D_real, D_shuffle, U -- then put
through the standing gates and a five-epitope LOEO. Only if U carries stable
independent information does it earn the sixteen-fold test.

**Q2 checked before running:** the shuffle permutes each position among cohort
members of equal length, so a singleton length stratum cannot be shuffled.
Coverage is 99.1-100% in every cohort including all three private ones
(CINGVCWTV 99.6%, LLWNGPMAV 99.1%, LTDEMIAQY 99.2%), so the null is
constructible on the private split at essentially full precision.

Pre-registered residual-on-D prediction, as the standard requires: floor is
0.5098 (a term already inside D), the project record is 0.5609. If the
two-layer hypothesis holds, U should land **0.55-0.60**. If U is only D
rescaled, it lands near the floor.

Run: .venv/bin/python scripts/118_excess_convergence.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/96_excess_convergence.md"
SIGMA = 0.1
FAMILY_RADIUS = 0.12
DRAWS = 5
SEED = 0
EPS = 1e-9
PREDICTED_RESIDUAL = (0.55, 0.60)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def slash(g):
    return None if pd.isna(g) else str(g).split("*")[0] + "*01"


def to_frame(trav, ca, traj, trbv, cb, trbj):
    from libtcrlm.schema.tcr import TravGene, TrbvGene
    f = pd.DataFrame({"TRAV": trav, "CDR3A": ca, "TRAJ": traj,
                      "TRBV": trbv, "CDR3B": cb, "TRBJ": trbj})
    for col, vocab in (("TRAV", {g.name for g in TravGene}),
                       ("TRBV", {g.name for g in TrbvGene})):
        f.loc[~f[col].astype(str).str.split("*").str[0].isin(vocab), col] = None
    return f


def positional_shuffle(seqs, rng):
    out = list(seqs)
    by_len = defaultdict(list)
    for i, s in enumerate(seqs):
        by_len[len(s)].append(i)
    for L, idx in by_len.items():
        if len(idx) < 2:
            continue
        M = np.array([list(seqs[i]) for i in idx])
        for c in range(L):
            M[:, c] = rng.permutation(M[:, c])
        for j, i in enumerate(idx):
            out[i] = "".join(M[j])
    return out


def main() -> int:
    import sceptr.variant as variants

    test = load_test().reset_index(drop=True)
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)
    ep = test["epitope_aa"].to_numpy()
    epitopes = sorted(set(ep))
    rng = np.random.default_rng(SEED)

    a = test["cdr3_alpha_aa"].fillna("").tolist()
    b = test["cdr3_beta_aa"].fillna("").tolist()
    trav = test["TRAV_IMGT"].map(slash).tolist()
    traj = test["TRAJ_IMGT"].map(slash).tolist()
    trbv = test["TRBV_IMGT"].map(slash).tolist()
    trbj = test["TRBJ_IMGT"].map(slash).tolist()

    model = variants.default()
    print("=== embedding the real cohorts ===")
    VA = model.calc_vector_representations(to_frame(trav, a, traj, trbv, b, trbj))
    VA = VA / np.linalg.norm(VA, axis=1, keepdims=True)

    # families from the real space, used identically for both conditions
    fam = np.arange(len(test))
    for e in epitopes:
        idx = np.where(ep == e)[0]
        sub = (1.0 - VA[idx] @ VA[idx].T) <= FAMILY_RADIUS
        lab = np.arange(len(idx))
        for _ in range(len(idx)):
            u = np.array([lab[sub[i]].min() for i in range(len(idx))])
            if (u == lab).all():
                break
            lab = u
        fam[idx] = idx[0] + lab
    print(f"  families at radius {FAMILY_RADIUS}: {len(set(fam))}")

    def density(V):
        """Leave-self-and-family-out kernel density within each cohort."""
        out = np.zeros(len(test))
        for e in epitopes:
            idx = np.where(ep == e)[0]
            Dm = 1.0 - V[idx] @ V[idx].T
            K = np.exp(-Dm / SIGMA)
            keep = fam[idx][:, None] != fam[idx][None, :]
            n = keep.sum(axis=1)
            out[idx] = np.where(n > 0, (K * keep).sum(axis=1)
                                / np.maximum(n, 1), np.nan)
        return out

    d_real = density(VA)
    print("\n=== embedding the shuffled nulls ===")
    acc = np.zeros(len(test))
    for k in range(DRAWS):
        ca, cb = list(a), list(b)
        for e in epitopes:
            idx = np.where(ep == e)[0]
            for src, dst in ((a, ca), (b, cb)):
                sh = positional_shuffle([src[i] for i in idx], rng)
                for j, i in enumerate(idx):
                    dst[i] = sh[j]
        V = model.calc_vector_representations(
            to_frame(trav, ca, traj, trbv, cb, trbj))
        V = V / np.linalg.norm(V, axis=1, keepdims=True)
        acc += density(V)
        print(f"  draw {k+1}/{DRAWS} done")
    d_shuf = acc / DRAWS
    U = np.log((d_real + EPS) / (d_shuf + EPS))

    ids = test["tcr_id"]
    F = {"D_real": pd.Series(d_real, index=ids),
         "D_shuffle": pd.Series(d_shuf, index=ids),
         "U excess convergence": pd.Series(U, index=ids)}

    print("\n=== 1. do the real and shuffled densities even differ? ===")
    print(f"  {'epitope':<12}{'split':<9}{'D_real':>11}{'D_shuffle':>12}"
          f"{'mean U':>9}{'sd U':>8}")
    for e in epitopes:
        i = np.where(ep == e)[0]
        print(f"  {e:<12}{test['split'].iloc[i[0]]:<9}{np.nanmean(d_real[i]):11.5f}"
              f"{np.nanmean(d_shuf[i]):12.5f}{np.nanmean(U[i]):9.3f}"
              f"{np.nanstd(U[i]):8.3f}")

    print("\n=== 2. the standing gates, with the E112 leak check ===")
    print(f"  pre-registered residual for U: {PREDICTED_RESIDUAL[0]:.2f}"
          f"-{PREDICTED_RESIDUAL[1]:.2f}   (floor 0.5098, record 0.5609)")
    print(f"  {'quantity':<24}{'macro':>8}{'typ rho':>9}{'priv':>7}"
          f"{'resid':>8}{'leak':>7}{'valid':>7}")
    rows = []
    for nm, s in F.items():
        v = s.copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.mean()).to_numpy()
        tg = typicality_gate(test, scored, v)
        pv, _ = private_variance_gate(test, v)
        m = hyp.mean_auc(labels, v)
        res = pd.Series(index=v.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in v.index]
            x = dr.reindex(i).rank(pct=True).to_numpy()
            y = v.reindex(i).rank(pct=True).to_numpy()
            r = y - np.polyval(np.polyfit(x, y, 1), x)
            res.loc[i] = r
            leak.append(abs(spearmanr(r, x).statistic))
        rr = hyp.mean_auc(labels, res)
        resid, mx = max(rr, 1 - rr), max(leak)
        rho = tg["rho_with_typicality"]
        print(f"  {nm:<24}{max(m,1-m):8.4f}{rho:+9.3f}"
              f"{'pass' if pv else 'FAIL':>7}{resid:8.4f}{mx:7.3f}"
              f"{'yes' if mx < 0.1 else 'NO':>7}")
        rows.append([nm, max(m, 1 - m), rho, pv, resid, mx])

    print("\n=== 3. blended into D, five-epitope LOEO ===")
    pub = sorted(scored["epitope_aa"].unique())
    base = hyp.mean_auc(labels, dr)
    WS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
    print(f"  D alone {base:.4f}")
    blends = []
    for nm, s in F.items():
        v = s.copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.mean()).to_numpy()
        for sign in (1.0, -1.0):
            sr = hyp.rank_within(sign * v, test)
            per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
            vals = {w: float(np.mean([per[w][e] for e in pub])) for w in WS}
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pick][e])
            nested = float(np.mean(held))
            bw = max(WS, key=lambda w: vals[w])
            if sign > 0 or vals[bw] > 0:
                print(f"  {nm:<24}sign {int(sign):+d}  oracle w={bw:.2f} "
                      f"{vals[bw]:.4f} ({vals[bw]-base:+.4f})   "
                      f"nested {nested:.4f} ({nested-base:+.4f})")
            blends.append([nm, sign, bw, vals[bw]-base, nested-base])

    lines = ["# E135 — Excess convergence over an epitope-conditioned null", "",
             "Generated by `scripts/118_excess_convergence.py`.", "",
             f"`U(t,e) = log((D_real + eps)/(E[D_shuffle] + eps))`, "
             f"{DRAWS} shuffle draws, σ={SIGMA}, "
             f"leave-self-and-family-out at radius {FAMILY_RADIUS}.", "",
             "Shuffle coverage is 99.1–100% in every cohort including all three",
             "private ones, so the null is constructible on the private split.", "",
             "| quantity | macro | typicality ρ | private var | residual on D | leak | valid |",
             "|---|---:|---:|---|---:|---:|---|"]
    for nm, m, rho, pv, resid, mx in rows:
        lines.append(f"| `{nm}` | {m:.4f} | {rho:+.3f} "
                     f"| {'pass' if pv else 'FAIL'} | **{resid:.4f}** "
                     f"| {mx:.3f} | {'yes' if mx<0.1 else 'NO'} |")
    lines += ["", f"Pre-registered residual prediction for `U`: "
              f"{PREDICTED_RESIDUAL[0]:.2f}–{PREDICTED_RESIDUAL[1]:.2f}.", "",
              "| quantity | sign | oracle w | oracle Δ | nested Δ |",
              "|---|---:|---:|---:|---:|"]
    for nm, sg, bw, od, nd in blends:
        lines.append(f"| `{nm}` | {int(sg):+d} | {bw:.2f} | {od:+.4f} | **{nd:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
