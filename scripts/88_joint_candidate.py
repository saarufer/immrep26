#!/usr/bin/env python
"""Do the three positive results add, or are they the same +0.001 three times?

Thirty-odd experiments produced three modifications with a *positive* nested
LOEO, which is three more than everything before them:

  E83  structure slot = 0.7 x pae + 0.3 x contact-fingerprint density  +0.0015
  E96  lenfreq (how common this CDR3 length pair is in the cohort)     +0.0005
  E97  per-epitope bandwidth, sigma_e = sigma_base x (gm nn1 / nn1_e)^0.5  +0.0007

Each is individually inside E15's -0.0048 selection cost and nowhere near the
+0.02 bar. The question is whether they are mechanistically independent — a
structural view, a sequence-length statistic and a kernel scale have no obvious
reason to move the same TCRs — in which case the sum is the thing to measure,
and measuring it **once jointly** matters: adding them one at a time pays a
selection cost per step, which is exactly how E15's -0.0048 is incurred.

Three things are reported and they answer different questions:

  1. **All eight combinations**, so interference is visible rather than assumed.
  2. **Nested LOEO with the combination re-chosen inside each fold**, which is
     what the selection would actually cost.
  3. **Nested LOEO with all three applied unconditionally**, no choice at all.
     This is the honest ceiling for a pre-registered change — though the three
     parameters were themselves tuned on these same five epitopes in E83/E96/E97,
     so even this is optimistic and is labelled as such.

And the part that decides whether any of it matters: **which modifications act
on the private rows at all.**

Run: .venv/bin/python scripts/88_joint_candidate.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/65_joint_candidate.md"
K = 8
FIXED = {"default": 0.1, "small": 0.1, "tcrdist": 20.0}
GRID = {"default": np.array([0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]),
        "small": np.array([0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]),
        "tcrdist": np.array([4.0, 6.0, 10.0, 14.0, 20.0, 30.0, 45.0, 70.0])}
P_BANDWIDTH = 0.5      # E97, chosen in all five folds
W_STRUCT_FP = 0.3      # E83
W_LENFREQ = 0.05       # E96


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    epitopes = sorted(test["epitope_aa"].unique())
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    # ---------- distance matrices ----------
    Dm = {}
    for space in ("default", "small"):
        c = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
        W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
        po = {t: i for i, t in enumerate(c["tcr_id"])}
        Dm[space] = {}
        for e, g in test.groupby("epitope_aa"):
            X = W[[po[t] for t in g["tcr_id"]]]
            D = 1.0 - X @ X.T
            np.fill_diagonal(D, np.inf)
            Dm[space][e] = (D, g["tcr_id"].to_numpy())
    Dm["tcrdist"] = {}
    for e, g in test.groupby("epitope_aa"):
        c = np.load(ROOT / f"data/interim/tcrdist/{e}.npz", allow_pickle=True)
        po = {t: i for i, t in enumerate(c["tcr_id"])}
        rows = [po[t] for t in g["tcr_id"]]
        D = (c["alpha"] + c["beta"]).astype(float)[np.ix_(rows, rows)]
        np.fill_diagonal(D, np.inf)
        Dm["tcrdist"][e] = (D, g["tcr_id"].to_numpy())

    nn1 = {sp: {e: float(np.median(np.sort(Dm[sp][e][0], axis=1)[:, 0]))
                for e in epitopes} for sp in Dm}

    def sigmas(space, p):
        sc = nn1[space]
        gm = float(np.exp(np.mean(np.log(list(sc.values())))))
        out = {}
        for e in epitopes:
            s = FIXED[space] * (gm / sc[e]) ** p
            out[e] = float(GRID[space][int(np.argmin(
                np.abs(np.log(GRID[space]) - np.log(max(s, 1e-12)))))])
        return out

    def density3(p) -> pd.Series:
        parts = []
        for sp in ("default", "small", "tcrdist"):
            sig = sigmas(sp, p)
            v = pd.Series(index=test["tcr_id"], dtype=float)
            for e in epitopes:
                D, ids = Dm[sp][e]
                v.loc[ids] = np.exp(-D / sig[e]).mean(axis=1)
            parts.append(rank(v))
        return sum(parts) / 3

    # ---------- the three modifications ----------
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    pae = rank(boltz["pae_tcr_pep_min"])
    fp = np.load(ROOT / "data/interim/contact_fingerprints.npz", allow_pickle=True)
    raw = fp["fingerprint"].astype(float)
    where = {t: i for i, t in enumerate(fp["tcr_id"])}
    normed = raw / (np.linalg.norm(raw, axis=1, keepdims=True) + 1e-9)
    vals = {}
    for e, g in test.groupby("epitope_aa"):
        ids = [t for t in g["tcr_id"] if t in where]
        B = normed[[where[t] for t in ids]]
        d = 1.0 - B @ B.T
        np.fill_diagonal(d, np.inf)
        for i, t in enumerate(ids):
            vals[t] = -float(np.sort(d[i])[:K].mean())
    fpd = rank(pd.Series(vals).reindex(test["tcr_id"]))
    struct_plain = pae
    struct_combo = (1 - W_STRUCT_FP) * pae + W_STRUCT_FP * fpd

    lf = pd.Series(index=test["tcr_id"], dtype=float)
    for e, g in test.groupby("epitope_aa"):
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        lf.loc[g["tcr_id"]] = (la.map(la.value_counts(normalize=True))
                               * lb.map(lb.value_counts(normalize=True))).to_numpy()
    lfr = rank(lf)

    # ---------- rebuild D exactly, then vary ----------
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = rank(d_sub)
    d3_base = density3(0.0)          # p=0 reproduces D's fixed bandwidths
    rest = dr - d3_base / 4 - struct_plain / 4     # m30_pooled + vdomain, /4 each

    def build(bw: bool, st: bool, lfq: bool) -> pd.Series:
        s = rest + (density3(P_BANDWIDTH) if bw else d3_base) / 4 \
            + (struct_combo if st else struct_plain) / 4
        return (1 - W_LENFREQ) * s + W_LENFREQ * lfr if lfq else s

    # sanity: the all-off build must reproduce D to numerical precision
    err = float(np.abs(build(False, False, False) - dr).max())
    print(f"reconstruction check: max |build(off,off,off) - D| = {err:.2e}")
    assert err < 1e-9, "D is not reproduced; the decomposition is wrong"

    names = {(False, False, False): "D (unchanged)",
             (True, False, False): "+ bandwidth (E97)",
             (False, True, False): "+ structure combo (E83)",
             (False, False, True): "+ lenfreq (E96)",
             (True, True, False): "+ bandwidth + structure",
             (True, False, True): "+ bandwidth + lenfreq",
             (False, True, True): "+ structure + lenfreq",
             (True, True, True): "ALL THREE"}
    cands = {names[k]: build(*k) for k in product([False, True], repeat=3)}
    per = {n: hyp.per_epitope(labels, s) for n, s in cands.items()}
    macro = {n: float(np.mean([per[n][e] for e in pub])) for n in cands}
    base = macro["D (unchanged)"]

    print("\n=== 1. all eight combinations ===")
    print(f"  {'combination':<28}{'macro':>9}{'vs D':>9}   per epitope")
    order = sorted(macro, key=lambda n: -macro[n])
    for n in order:
        print(f"  {n:<28}{macro[n]:9.4f}{macro[n]-base:+9.4f}   "
              + "  ".join(f"{e[:3]} {per[n][e]-per['D (unchanged)'][e]:+.4f}"
                          for e in pub))
    add = (macro["+ bandwidth (E97)"] - base) + \
          (macro["+ structure combo (E83)"] - base) + \
          (macro["+ lenfreq (E96)"] - base)
    print(f"\n  sum of the three individual gains : {add:+.4f}")
    print(f"  the joint gain actually measured   : {macro['ALL THREE']-base:+.4f}"
          f"   ({(macro['ALL THREE']-base)/add:.0%} of additive)")

    # ---------- 2. do they move the same TCRs? ----------
    print("\n=== 2. are the three mechanistically independent? ===")
    deltas = {"bandwidth": cands["+ bandwidth (E97)"] - dr,
              "structure": cands["+ structure combo (E83)"] - dr,
              "lenfreq": cands["+ lenfreq (E96)"] - dr}
    print(f"  within-epitope Spearman between the per-TCR score changes:")
    ks = list(deltas)
    pairs = []
    for i in range(3):
        for j in range(i + 1, 3):
            r = float(np.mean([
                spearmanr(deltas[ks[i]].reindex(g["tcr_id"]),
                          deltas[ks[j]].reindex(g["tcr_id"])).statistic
                for _, g in test.groupby("epitope_aa")]))
            print(f"    {ks[i]:<10} vs {ks[j]:<10} {r:+.3f}")
            pairs.append([ks[i], ks[j], r])

    # ---------- 3. nested LOEO, two ways ----------
    print("\n=== 3. nested LOEO ===")
    held_sel, chosen = [], []
    for e in pub:
        tr = [x for x in pub if x != e]
        best = max(cands, key=lambda n: np.mean([per[n][t] for t in tr]))
        chosen.append(best)
        held_sel.append(per[best][e])
    sel = float(np.mean(held_sel))
    fixed_all = float(np.mean([per["ALL THREE"][e] for e in pub]))
    print(f"  combination re-chosen inside each fold : {sel:.4f} "
          f"({sel-base:+.4f})")
    for e, c in zip(pub, chosen):
        print(f"      {e:<12} chose {c}")
    print(f"  all three applied unconditionally      : {fixed_all:.4f} "
          f"({fixed_all-base:+.4f})   <- optimistic, see header")

    # ---------- 3b. the test that actually decides it ----------
    # Section 3 holds the three parameters at the values E83/E96/E97 chose --
    # on these same five epitopes. The defensible number re-selects all three
    # inside every fold, so the held-out epitope has never seen any of them.
    print("\n=== 3b. fully nested: all three parameters re-chosen inside each fold ===")
    P_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
    S_GRID = [0.0, 0.15, 0.3, 0.5]
    L_GRID = [0.0, 0.05, 0.10]
    d3_cache = {p_: density3(p_) for p_ in P_GRID}
    grid_per = {}
    for p_, sw, lw in product(P_GRID, S_GRID, L_GRID):
        sc = rest + d3_cache[p_] / 4 + ((1 - sw) * pae + sw * fpd) / 4
        sc = (1 - lw) * sc + lw * lfr if lw else sc
        grid_per[(p_, sw, lw)] = hyp.per_epitope(labels, sc)
    print(f"  {len(grid_per)} parameter combinations "
          f"({len(P_GRID)}x{len(S_GRID)}x{len(L_GRID)})")
    held_full, picks_full = [], []
    for e in pub:
        tr = [x for x in pub if x != e]
        best = max(grid_per, key=lambda k: np.mean([grid_per[k][t] for t in tr]))
        picks_full.append((e, best))
        held_full.append(grid_per[best][e])
        print(f"    {e:<12} chose p={best[0]:g} struct_fp={best[1]:g} "
              f"lenfreq={best[2]:g}  ->  {grid_per[best][e]:.4f}")
    full_nested = float(np.mean(held_full))
    oracle_all = max(grid_per, key=lambda k: np.mean(
        [grid_per[k][t] for t in pub]))
    oracle_val = float(np.mean([grid_per[oracle_all][t] for t in pub]))
    print(f"  fully nested macro      {full_nested:.4f}  ({full_nested-base:+.4f})")
    print(f"  oracle over the grid    {oracle_val:.4f}  ({oracle_val-base:+.4f})"
          f"   at p={oracle_all[0]:g} struct_fp={oracle_all[1]:g} "
          f"lenfreq={oracle_all[2]:g}")
    print(f"  selection cost of the three parameters: "
          f"{full_nested-oracle_val:+.4f}")

    # ---------- 4. what survives to the private split ----------
    print("\n=== 4. which modifications act on the private rows at all? ===")
    print(f"  {'epitope':<12}{'split':<9}{'bandwidth sigma':>18}"
          f"{'fingerprint':>13}{'lenfreq levels':>16}")
    sig05 = sigmas("default", P_BANDWIDTH)
    priv_rows = []
    sp_of = test.drop_duplicates("epitope_aa").set_index("epitope_aa")["split"]
    for e in epitopes:
        g = test[test["epitope_aa"] == e]
        have_fp = sum(t in where for t in g["tcr_id"])
        lv = lf.reindex(g["tcr_id"]).nunique()
        changed = "changed" if abs(sig05[e] - FIXED["default"]) > 1e-9 else "UNCHANGED"
        print(f"  {e:<12}{sp_of[e]:<9}{sig05[e]:>10g} {changed:<8}"
              f"{have_fp:>8}/{len(g):<4}{lv:>16}")
        priv_rows.append([e, sp_of[e], sig05[e], changed, have_fp, len(g), lv])

    lines = ["# E98 — Do the three positive results add?", "",
             "Generated by `scripts/88_joint_candidate.py`. The all-off build",
             f"reproduces the submitted D to {err:.0e}, asserted in code.", "",
             "## 1. All eight combinations", "",
             "| combination | macro | vs D |", "|---|---:|---:|"]
    for n in order:
        lines.append(f"| {n} | {macro[n]:.4f} | **{macro[n]-base:+.4f}** |")
    lines += ["",
              f"Sum of the three individual gains: **{add:+.4f}**. "
              f"Joint gain measured: **{macro['ALL THREE']-base:+.4f}** "
              f"({(macro['ALL THREE']-base)/add:.0%} of additive).", "",
              "## 2. Are they independent?", "",
              "Within-epitope Spearman between the per-TCR score changes each",
              "modification induces:", "",
              "| a | b | ρ |", "|---|---|---:|"]
    for a, b, r in pairs:
        lines.append(f"| {a} | {b} | {r:+.3f} |")
    lines += ["", "## 3. Nested LOEO", "",
              f"- combination re-chosen inside each fold: **{sel:.4f}** "
              f"({sel-base:+.4f})",
              f"- all three applied unconditionally: **{fixed_all:.4f}** "
              f"({fixed_all-base:+.4f}) — optimistic, since the three parameters",
              "  were themselves tuned on these five epitopes", "",
              "| held out | combination chosen on the other four |", "|---|---|"]
    for e, c in zip(pub, chosen):
        lines.append(f"| {e} | {c} |")
    lines += ["", "### 3b. Fully nested — all three parameters re-chosen per fold", "",
              "Section 3 holds the parameters at the values E83/E96/E97 chose on",
              "these same five epitopes. This re-selects all three inside every",
              f"fold over a {len(grid_per)}-point grid.", "",
              "| held out | p | struct_fp | lenfreq | AUC |",
              "|---|---:|---:|---:|---:|"]
    for e, b in picks_full:
        lines.append(f"| {e} | {b[0]:g} | {b[1]:g} | {b[2]:g} "
                     f"| {grid_per[b][e]:.4f} |")
    lines += ["",
              f"- fully nested macro: **{full_nested:.4f}** ({full_nested-base:+.4f})",
              f"- oracle over the same grid: {oracle_val:.4f} "
              f"({oracle_val-base:+.4f}) at p={oracle_all[0]:g}, "
              f"struct_fp={oracle_all[1]:g}, lenfreq={oracle_all[2]:g}",
              f"- selection cost of the three parameters: "
              f"**{full_nested-oracle_val:+.4f}**", "",
              "## 4. What reaches the private split", "",
              "| epitope | split | bandwidth σ | | fingerprint coverage | lenfreq levels |",
              "|---|---|---:|---|---:|---:|"]
    for r in priv_rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:g} | {r[3]} | {r[4]}/{r[5]} "
                     f"| {r[6]} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
