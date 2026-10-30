#!/usr/bin/env python
"""E146 -- Sparse-reference stress test, and whether the pool term should shrink.

An external review's highest-priority experiment, and it names a hole this
project's own audit missed:

> **If a component's information is already weak, percentile ranking still
> stretches it into a full 0-1 distribution.**

E123 verified that all four of D's terms carry *variance* in all three private
cohorts and treated that as reassurance. Variance only rules out a constant. It
does not show the term deserves its weight. D gives `m30pool` **a flat quarter**
of the rank average, and the private pools hold **11** receptors (CINGVCWTV),
**191** (LLWNGPMAV) and **5** (LTDEMIAQY), against 3,005 and 9,881 on the
public cohorts the weight was validated on.

So the question is not "does the pool term vary on private" but "does a pool of
5 to 11 receptors support a quarter of the model".

**Design.** The public epitopes are forced into the private regime by
subsampling their reference pools, and the same queries are re-scored:

  pool sizes   5, 11, 25, 50, 100, all     (5 and 11 match LTD and CIN exactly)
  repeats      8 independent draws per size
  variants     D as shipped; F (pool deleted); D with a pre-fixed shrinkage

**The shrinkage, fixed before any outcome is seen.** E84 measured
independently, long before this experiment, that the pool term is dead below
~60 same-epitope rows. That constant is reused rather than fitted:

    lambda(n) = min(1, n / 60)
    r_pool* = 0.5 + lambda(n) * (r_pool - 0.5)

and, per the review's warning, **the shrunk component is NOT re-ranked** --
re-percentiling would undo the shrinkage exactly.

**Pre-registered.** If the pool term is actively harmful at private pool sizes,
shrinkage should beat D at n=5 and n=11 and converge to D at large n. Predicted
gain at n<=11: **+0.002 to +0.010**. If instead shrinkage never beats D, the
flat quarter is defensible and the review's concern, though correct in
principle, does not bind here.

Run: .venv/bin/python scripts/124_sparse_reference_stress.py
"""

from __future__ import annotations

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

REPORT = ROOT / "reports/102_sparse_reference_stress.md"
SIZES = [5, 11, 25, 50, 100, None]
REPEATS = 8
TOPK = 8
SEED = 0
N0 = 60.0          # E84's measured threshold, reused not fitted
PREDICTED = (0.002, 0.010)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    import sceptr
    sceptr.disable_hardware_acceleration()
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    pm = load_module("pool", "scripts/33_pool_negative_control.py")

    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    lab = labels.dropna(subset=["pos"]).copy()
    lab["pos"] = lab["pos"].astype(int)
    pub = sorted(lab["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s.reindex(test["tcr_id"]), test)  # noqa: E731

    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv"
                        ).set_index("tcr_id")
    structure = rank(boltz["pae_tcr_pep_min"])
    vdomain = rank(rival.vdomain_topk(test))
    frozen = decomp.submitted_score(test)
    density3 = rank((frozen - 0.3 * structure) / 0.7)
    sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                      ).set_index("tcr_id")["prediction"]
    # D's shipped m30 term, recovered the same way E123 verified to 1e-16
    m30_full = rank(4 * sub.reindex(test["tcr_id"]) - density3 - vdomain - structure)

    print("=== embedding each cohort's queries and its reference pool ===")
    pool = pm.load_pool()
    Q, EX = {}, {}
    for e, g in test.groupby("epitope_aa"):
        g = g.reset_index(drop=True)
        fr = pd.DataFrame({"cdr3a": g["cdr3_alpha_aa"], "cdr3b": g["cdr3_beta_aa"],
                           "trav": g["TRAV_IMGT"], "traj": g["TRAJ_IMGT"],
                           "trbv": g["TRBV_IMGT"], "trbj": g["TRBJ_IMGT"], "ep": e})
        Q[e] = (g, pm.embed(fr))
        keys = set(zip(g["cdr3_alpha_aa"], g["cdr3_beta_aa"]))
        same = pool[pool["ep"] == e]
        ex = same[[(a, b) not in keys for a, b in
                   zip(same["cdr3a"], same["cdr3b"])]].drop_duplicates(
                       ["cdr3a", "cdr3b"]).reset_index(drop=True)
        EX[e] = (ex, pm.embed(ex) if len(ex) else np.zeros((0, Q[e][1].shape[1])))
        print(f"  {e:<12}{len(g):4d} queries, {len(ex):5d} pool receptors")

    def m30_with(e, idx):
        """m30 for cohort e using only the pool rows in idx (raw, not ranked)."""
        g, vec = Q[e]
        ex, exv = EX[e]
        ref = vec if len(idx) == 0 else np.vstack([vec, exv[idx]])
        d = 1.0 - vec @ ref.T
        d[np.arange(len(g)), np.arange(len(g))] = np.inf
        k = min(TOPK, d.shape[1] - 1)
        return pd.Series(-np.sort(d, axis=1)[:, :k].mean(axis=1),
                         index=g["tcr_id"])

    rng = np.random.default_rng(SEED)
    base_per = hyp.per_epitope(labels, rank(
        (density3 + m30_full + vdomain + structure) / 4))
    print(f"\n=== D as shipped, public macro "
          f"{np.mean([base_per[e] for e in pub]):.4f} ===")
    print(f"  pre-registered shrinkage gain at n<=11: "
          f"{PREDICTED[0]:+.3f} to {PREDICTED[1]:+.3f}")

    print("\n=== forcing each public cohort into the private regime ===")
    print(f"  lambda(n) = min(1, n/60), E84's threshold, reused not fitted")
    print(f"  {'pool n':>8}{'D (pool kept)':>15}{'F (pool deleted)':>18}"
          f"{'D + shrinkage':>15}{'shrink - D':>12}")
    out = []
    for n in SIZES:
        accD, accF, accS = [], [], []
        reps = 1 if n is None else REPEATS
        for _ in range(reps):
            m30 = pd.Series(index=test["tcr_id"], dtype=float)
            m30s = pd.Series(index=test["tcr_id"], dtype=float)
            for e in pub + [x for x in Q if x not in pub]:
                ex, _ = EX[e]
                navail = len(ex)
                take = navail if n is None else min(n, navail)
                idx = (np.arange(navail) if take == navail
                       else rng.choice(navail, take, replace=False))
                raw = m30_with(e, idx)
                r = raw.rank(pct=True)
                m30.loc[r.index] = r.to_numpy()
                lam = min(1.0, take / N0)
                m30s.loc[r.index] = (0.5 + lam * (r - 0.5)).to_numpy()
            dD = rank((density3 + m30 + vdomain + structure) / 4)
            dS = rank((density3 + m30s + vdomain + structure) / 4)
            pD = hyp.per_epitope(labels, dD)
            pS = hyp.per_epitope(labels, dS)
            accD.append(float(np.mean([pD[e] for e in pub])))
            accS.append(float(np.mean([pS[e] for e in pub])))
        # F: pool deleted entirely, independent of n
        dF = rank((density3 + rank(pd.concat(
            [m30_with(e, np.array([], dtype=int)) for e in Q])) + vdomain
            + structure) / 4)
        pF = hyp.per_epitope(labels, dF)
        accF = [float(np.mean([pF[e] for e in pub]))]
        lbl = "all" if n is None else str(n)
        print(f"  {lbl:>8}{np.mean(accD):15.4f}{np.mean(accF):18.4f}"
              f"{np.mean(accS):15.4f}{np.mean(accS)-np.mean(accD):+12.4f}")
        out.append([lbl, np.mean(accD), np.std(accD), np.mean(accF),
                    np.mean(accS), np.std(accS)])

    print("\n=== the shipped configuration, with shrinkage at the TRUE pool sizes ===")
    print("  (each cohort keeps its own real pool; only lambda changes)")
    m30s = pd.Series(index=test["tcr_id"], dtype=float)
    for e in Q:
        ex, _ = EX[e]
        raw = m30_with(e, np.arange(len(ex)))
        r = raw.rank(pct=True)
        lam = min(1.0, len(ex) / N0)
        m30s.loc[r.index] = (0.5 + lam * (r - 0.5)).to_numpy()
    dS = rank((density3 + m30s + vdomain + structure) / 4)
    pS = hyp.per_epitope(labels, dS)
    print(f"  {'epitope':<12}{'pool':>7}{'lambda':>8}{'D':>9}{'shrunk':>9}{'diff':>9}")
    for e in pub:
        ex, _ = EX[e]
        print(f"  {e:<12}{len(ex):7d}{min(1.0,len(ex)/N0):8.2f}"
              f"{base_per[e]:9.4f}{pS[e]:9.4f}{pS[e]-base_per[e]:+9.4f}")
    mb = float(np.mean([base_per[e] for e in pub]))
    ms = float(np.mean([pS[e] for e in pub]))
    print(f"  {'MACRO':<12}{'':7}{'':8}{mb:9.4f}{ms:9.4f}{ms-mb:+9.4f}")
    print(f"\n  lambda on the PRIVATE cohorts: CINGVCWTV "
          f"{min(1.0,11/N0):.2f}, LLWNGPMAV {min(1.0,191/N0):.2f}, "
          f"LTDEMIAQY {min(1.0,5/N0):.2f}")

    lines = ["# E146 — Sparse-reference stress test and pool shrinkage", "",
             "Generated by `scripts/124_sparse_reference_stress.py`.", "",
             "The review's point: percentile ranking stretches a weak component",
             "into a full 0–1 distribution, so E123's finding that all four of",
             "D's terms carry *variance* on the private cohorts rules out a",
             "constant but does not show the pool term deserves its flat quarter",
             "when the private pools hold 11, 191 and 5 receptors.", "",
             f"Shrinkage `λ(n) = min(1, n/60)` with E84's threshold reused, not",
             "fitted; the shrunk component is **not** re-ranked.", "",
             f"Pre-registered gain at n≤11: {PREDICTED[0]:+.3f} to {PREDICTED[1]:+.3f}.",
             "", "| pool n | D (pool kept) | F (pool deleted) | D + shrinkage | shrink − D |",
             "|---|---:|---:|---:|---:|"]
    for lbl, dm, ds, fm, sm, ss in out:
        lines.append(f"| {lbl} | {dm:.4f} ± {ds:.4f} | {fm:.4f} "
                     f"| {sm:.4f} ± {ss:.4f} | **{sm-dm:+.4f}** |")
    lines += ["", "## Shrinkage at the true pool sizes", "",
              "| epitope | pool | λ | D | shrunk | diff |",
              "|---|---:|---:|---:|---:|---:|"]
    for e in pub:
        ex, _ = EX[e]
        lines.append(f"| {e} | {len(ex)} | {min(1.0,len(ex)/N0):.2f} "
                     f"| {base_per[e]:.4f} | {pS[e]:.4f} "
                     f"| {pS[e]-base_per[e]:+.4f} |")
    lines.append(f"| **MACRO** | | | **{mb:.4f}** | **{ms:.4f}** "
                 f"| **{ms-mb:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
