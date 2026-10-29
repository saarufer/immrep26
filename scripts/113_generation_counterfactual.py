#!/usr/bin/env python
"""E127 -- V(D)J generation-background counterfactual for the E88 attribution.

E88 asked whether the density mechanism can pick the right *peptide* for a TCR
(task B) rather than only judge a given annotation (task A). This experiment
asks how much of whatever E88 measured is explained by recombination background
and positional composition, rather than by epitope-specific higher-order
structure.

E88's setup is frozen exactly: 8 epitopes, SIGMA 0.1, N_SUBSAMPLE 110,
REPEATS 30, SEED 0, family radius 0.12, per-column standardisation over the
QUERY set (the asymmetric version produced a degenerate predictor that sent 62%
of queries to NLVPMVATV).

Only the **reference libraries** change. Each is a 1:1 replacement for the real
cohort, so the exclusion mask, the subsample sizes and the surviving reference
counts are identical across conditions -- the comparison is controlled by
construction:

  A  the real epitope cohorts (E88)
  B  V and J genes kept from the real member; both CDR3s replaced by
     OLGA-generated sequences of the **same length**. V/J usage, CDR3 length
     and recombination background therefore match the real cohort, while
     epitope-specific CDR3 content is destroyed.
  C  V and J genes kept; each CDR3 position-wise shuffled across cohort
     members of the same length and chain. The cohort's **positional
     amino-acid composition is preserved exactly** and higher-order motifs
     (residue co-occurrence) are destroyed.

Reading, stated before the run:

  A ~ C   simple generation and composition features explain most of E88's
          attribution signal.
  A >> C  epitope-related higher-order structure exists -- though thymic
          selection, study provenance and sampling remain uncontrolled, so
          this would still not be a causal claim about antigen selection.

E89 already bears on the answer: an OLGA-trained SCEPTR separates the cohorts
at 2.67x against the real model's 2.60x, i.e. condition A does not beat a
generative background on the *separation* measure. This tests the attribution
measure instead, and adds the composition-matched rung that E89 lacked.

Run: .venv/bin/python scripts/113_generation_counterfactual.py
"""

from __future__ import annotations

import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

REPORT = ROOT / "reports/91_generation_counterfactual.md"
SIGMA = 0.1
N_SUBSAMPLE = 110
REPEATS = 30
SEED = 0
FAMILY_RADIUS = 0.12


def slash(g):
    return None if pd.isna(g) else str(g).split("*")[0] + "*01"


def to_frame(trav, cdr3a, traj, trbv, cdr3b, trbj) -> pd.DataFrame:
    from libtcrlm.schema.tcr import TravGene, TrbvGene
    f = pd.DataFrame({"TRAV": trav, "CDR3A": cdr3a, "TRAJ": traj,
                      "TRBV": trbv, "CDR3B": cdr3b, "TRBJ": trbj})
    for col, vocab in (("TRAV", {g.name for g in TravGene}),
                       ("TRBV", {g.name for g in TrbvGene})):
        bad = ~f[col].astype(str).str.split("*").str[0].isin(vocab)
        f.loc[bad, col] = None
    return f


def olga_pool(chain_folder: str, recomb: str, need: Counter, rng) -> dict:
    """OLGA CDR3s bucketed by length, enough of each length that is needed."""
    from tcrdist.pgen import OlgaModel
    m = OlgaModel(chain_folder=chain_folder, recomb_type=recomb)
    by_len: dict[int, list[str]] = defaultdict(list)
    tries = 0
    while tries < 60:
        for s in m.gen_cdr3s(n=50000):
            if s:
                by_len[len(s)].append(s)
        if all(len(by_len[L]) >= n for L, n in need.items()):
            break
        tries += 1
    return by_len


def positional_shuffle(seqs: list[str], rng) -> list[str]:
    """Shuffle each position independently among sequences of equal length."""
    out = list(seqs)
    by_len: dict[int, list[int]] = defaultdict(list)
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
    ep = test["epitope_aa"].to_numpy()
    epitopes = sorted(set(ep))
    rng = np.random.default_rng(SEED)
    print(f"=== {len(test)} TCRs, {len(epitopes)} epitopes ===")
    for e in epitopes:
        print(f"  {e:<12}{int((ep==e).sum()):5d}  {test.loc[ep==e,'split'].iloc[0]}")

    a_cdr3a = test["cdr3_alpha_aa"].fillna("").tolist()
    a_cdr3b = test["cdr3_beta_aa"].fillna("").tolist()
    trav = test["TRAV_IMGT"].map(slash).tolist()
    traj = test["TRAJ_IMGT"].map(slash).tolist()
    trbv = test["TRBV_IMGT"].map(slash).tolist()
    trbj = test["TRBJ_IMGT"].map(slash).tolist()

    print("\n=== building condition B: OLGA background at matched V/J and length ===")
    need_a = Counter(len(s) for s in a_cdr3a if s)
    need_b = Counter(len(s) for s in a_cdr3b if s)
    pa = olga_pool("human_T_alpha", "VJ", need_a, rng)
    pb = olga_pool("human_T_beta", "VDJ", need_b, rng)
    miss = 0

    def draw(pool, L):
        nonlocal miss
        c = pool.get(L)
        if not c:
            miss += 1
            return None
        return c[rng.integers(len(c))]

    b_cdr3a = [draw(pa, len(s)) or s if s else s for s in a_cdr3a]
    b_cdr3b = [draw(pb, len(s)) or s if s else s for s in a_cdr3b]
    print(f"  lengths OLGA could not supply: {miss} of {2*len(test)}")
    print(f"  identical to real by chance: alpha "
          f"{sum(x==y for x,y in zip(a_cdr3a,b_cdr3a))}, beta "
          f"{sum(x==y for x,y in zip(a_cdr3b,b_cdr3b))}")

    print("\n=== building condition C: positional shuffle within epitope x length ===")
    c_cdr3a, c_cdr3b = list(a_cdr3a), list(a_cdr3b)
    for e in epitopes:
        idx = np.where(ep == e)[0]
        for src, dst in ((a_cdr3a, c_cdr3a), (a_cdr3b, c_cdr3b)):
            sub = [src[i] for i in idx]
            sh = positional_shuffle(sub, rng)
            for j, i in enumerate(idx):
                dst[i] = sh[j]
    print(f"  unchanged after shuffle: alpha "
          f"{sum(x==y for x,y in zip(a_cdr3a,c_cdr3a))}, beta "
          f"{sum(x==y for x,y in zip(a_cdr3b,c_cdr3b))}")
    # verify C preserves positional composition exactly
    ok = True
    for e in epitopes[:3]:
        idx = np.where(ep == e)[0]
        for src, dst in ((a_cdr3a, c_cdr3a), (a_cdr3b, c_cdr3b)):
            for L in {len(src[i]) for i in idx if src[i]}:
                A = sorted("".join(src[i] for i in idx if len(src[i]) == L))
                B = sorted("".join(dst[i] for i in idx if len(dst[i]) == L))
                ok &= A == B
    print(f"  positional composition preserved in C: {ok}")

    print("\n=== embedding all three with sceptr.default ===")
    model = variants.default()
    E = {}
    for nm, (ca, cb) in (("A", (a_cdr3a, a_cdr3b)), ("B", (b_cdr3a, b_cdr3b)),
                         ("C", (c_cdr3a, c_cdr3b))):
        f = to_frame(trav, ca, traj, trbv, cb, trbj)
        V = model.calc_vector_representations(f)
        V = V / np.linalg.norm(V, axis=1, keepdims=True)
        E[nm] = V
        print(f"  {nm}: {V.shape}")

    # family labels and the exclusion mask come from A, used for all conditions
    DA = 1.0 - E["A"] @ E["A"].T
    fam = np.arange(len(test))
    for e in epitopes:
        idx = np.where(ep == e)[0]
        sub = DA[np.ix_(idx, idx)] <= FAMILY_RADIUS
        lab = np.arange(len(idx))
        for _ in range(len(idx)):
            upd = np.array([lab[sub[i]].min() for i in range(len(idx))])
            if (upd == lab).all():
                break
            lab = upd
        fam[idx] = idx[0] + lab
    print(f"\n  families at radius {FAMILY_RADIUS}: {len(set(fam))} "
          f"for {len(test)} TCRs")

    members = {e: np.where(ep == e)[0] for e in epitopes}
    query = np.arange(len(test))

    def attribute(V, permute=False):
        D = 1.0 - V @ V.T
        true = rng.permutation(ep[query]) if permute else ep[query]
        hits, ranks = [], []
        per = {e: [] for e in epitopes}
        for _ in range(REPEATS):
            pool = {e: rng.choice(members[e],
                                  size=min(N_SUBSAMPLE, len(members[e])),
                                  replace=False) for e in epitopes}
            Z = np.zeros((len(query), len(epitopes)))
            for j, e in enumerate(epitopes):
                ref = pool[e]
                d = D[np.ix_(query, ref)]
                keep = (query[:, None] != ref[None, :]) & \
                       (fam[query][:, None] != fam[ref][None, :])
                k = np.exp(-d / SIGMA)
                n = keep.sum(axis=1)
                Z[:, j] = np.where(n > 0, (k * keep).sum(axis=1)
                                   / np.maximum(n, 1), np.nan)
            Z = (Z - np.nanmean(Z, axis=0)) / (np.nanstd(Z, axis=0) + 1e-12)
            pred = np.array(epitopes)[np.nanargmax(Z, axis=1)]
            hits.append(pred == true)
            order = np.argsort(-Z, axis=1)
            ti = np.array([epitopes.index(t) for t in true])
            ranks.append(np.array([np.where(order[i] == ti[i])[0][0] + 1
                                   for i in range(len(query))]))
            for e in epitopes:
                m = true == e
                per[e].append(float((pred[m] == e).mean()) if m.any() else np.nan)
        H = np.mean([h.mean() for h in hits])
        R = np.mean([r.mean() for r in ranks])
        return H, R, {e: float(np.nanmean(per[e])) for e in epitopes}

    print("\n=== attribution: argmax accuracy and true-epitope rank ===")
    print(f"  (chance accuracy 1/{len(epitopes)} = {1/len(epitopes):.3f}, "
          f"chance rank {(len(epitopes)+1)/2:.1f})")
    out = {}
    for nm in ("A", "B", "C"):
        H, R, per = attribute(E[nm])
        out[nm] = (H, R, per)
        print(f"  {nm}  accuracy {H:.4f}   mean true rank {R:.3f}")
    Hp, Rp, _ = attribute(E["A"], permute=True)
    print(f"  null (A, labels permuted)  accuracy {Hp:.4f}   rank {Rp:.3f}")

    print("\n=== per-epitope accuracy ===")
    print(f"  {'epitope':<12}{'split':<9}{'A':>8}{'B':>8}{'C':>8}{'A-C':>8}")
    for e in epitopes:
        sp = test.loc[ep == e, "split"].iloc[0]
        print(f"  {e:<12}{sp:<9}{out['A'][2][e]:8.3f}{out['B'][2][e]:8.3f}"
              f"{out['C'][2][e]:8.3f}{out['A'][2][e]-out['C'][2][e]:+8.3f}")

    print("\n=== decomposition ===")
    hA, hB, hC = out["A"][0], out["B"][0], out["C"][0]
    span = hA - Hp
    print(f"  A {hA:.4f}   B {hB:.4f}   C {hC:.4f}   null {Hp:.4f}")
    print(f"  above null: A {hA-Hp:+.4f}  B {hB-Hp:+.4f}  C {hC-Hp:+.4f}")
    if span > 1e-9:
        print(f"  share of A's signal reached by V/J+length+background (B): "
              f"{(hB-Hp)/span*100:5.1f}%")
        print(f"  share reached by + positional composition (C):            "
              f"{(hC-Hp)/span*100:5.1f}%")
        print(f"  left for higher-order motif structure (A-C):              "
              f"{(hA-hC)/span*100:5.1f}%")

    lines = ["# E127 — V(D)J generation-background counterfactual", "",
             "Generated by `scripts/113_generation_counterfactual.py`.", "",
             "E88's setup frozen (σ 0.1, subsample 110, 30 repeats, family",
             "radius 0.12, per-column standardisation over the query set).",
             "Only the reference library changes; each is a 1:1 replacement so",
             "the exclusion mask and surviving reference counts match.", "",
             "| condition | library | accuracy | mean true rank |",
             "|---|---|---:|---:|",
             f"| **A** | real cohorts | **{hA:.4f}** | {out['A'][1]:.3f} |",
             f"| **B** | OLGA background, matched V/J + length | **{hB:.4f}** | {out['B'][1]:.3f} |",
             f"| **C** | positional shuffle, motifs destroyed | **{hC:.4f}** | {out['C'][1]:.3f} |",
             f"| null | A with permuted labels | {Hp:.4f} | {Rp:.3f} |", "",
             f"Chance accuracy 1/{len(epitopes)} = {1/len(epitopes):.3f}.", "",
             "## Per-epitope accuracy", "",
             "| epitope | split | A | B | C | A−C |", "|---|---|---:|---:|---:|---:|"]
    for e in epitopes:
        sp = test.loc[ep == e, "split"].iloc[0]
        lines.append(f"| {e} | {sp} | {out['A'][2][e]:.3f} | {out['B'][2][e]:.3f} "
                     f"| {out['C'][2][e]:.3f} | {out['A'][2][e]-out['C'][2][e]:+.3f} |")
    if span > 1e-9:
        lines += ["", "## Decomposition of A's signal above the null", "",
                  f"- V/J + length + recombination background (B): "
                  f"**{(hB-Hp)/span*100:.1f}%**",
                  f"- + positional amino-acid composition (C): "
                  f"**{(hC-Hp)/span*100:.1f}%**",
                  f"- higher-order motif structure (A−C): "
                  f"**{(hA-hC)/span*100:.1f}%**", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
