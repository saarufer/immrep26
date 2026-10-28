#!/usr/bin/env python
"""An independent structural replicate, used as a noise model rather than a feature.

Two complete Boltz runs of all 2,636 complexes exist on disk, `out` and `out5`.
Every structural result in this project — E36, E75, E83, E103, E104, E109 —
used only `out5`, and E19's consistency work resampled the *diffusion head*
within one run, which shares the trunk. Comparing the two runs is different:
on 25 sampled complexes the TCR C-alpha RMSD after superposing on MHC and
peptide has median **1.19 A** and reaches 5.76 A, with none identical. These
are genuinely independent predictions.

That is worth more as a **noise model** than as another feature. `polar_frac`
is computed from one structure; if that structure is uncertain, the value is
noisy, and two things follow that are both testable:

  1. **Averaging the two replicates should sharpen it** — if the signal is
     real. Two independent noisy estimates of the same quantity average to a
     cleaner one. If instead the AUC is unchanged or worse, `polar_frac` is
     measuring something specific to a single predicted pose rather than a
     property of the receptor.

  2. **Disagreement should predict unreliability.** A TCR whose two folds
     disagree has a less trustworthy `polar_frac`, and shrinking those values
     toward the cohort mean should help more than shrinking everything.

This is the first time the project has had a replicate to test a structural
feature against, and E104's `polar_frac` is the only structural feature that
ever carried anything, so it is the one to test.

Run: .venv/bin/python scripts/104_replicate_structures.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/81_replicate_structures.md"
CACHE = ROOT / "data/interim/polar_frac_run_out.csv"
RUN_A = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
NEAR = 9.0


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def polar_frac(path: Path):
    """Exactly E104's definition, on whichever run's PDB it is given."""
    xyz, el, ch = [], [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            e = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if e == "H":
                continue
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            el.append(e)
            ch.append(line[21])
    if not xyz:
        return np.nan
    xyz = np.asarray(xyz, float)
    el, ch = np.asarray(el), np.asarray(ch)
    tcr = np.isin(ch, ["A", "B"])
    pmhc = np.isin(ch, ["M", "P"])
    if not (tcr.any() and pmhc.any()):
        return np.nan
    t = cKDTree(xyz[pmhc])
    ti = np.where(tcr)[0]
    near = np.array([len(x) > 0 for x in t.query_ball_point(xyz[ti], r=NEAR)])
    if not near.any():
        return np.nan
    return float(np.isin(el, ["N", "O"])[ti][near].mean())


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))

    if CACHE.exists():
        print(f"loading {CACHE.name}")
        A = pd.read_csv(CACHE).set_index("tcr_id")["polar_frac_out"]
    else:
        vals, ids = [], test["tcr_id"].tolist()
        for i, t in enumerate(ids):
            p = RUN_A / t / f"{t}_model_0.pdb"
            vals.append(polar_frac(p) if p.exists() else np.nan)
            if (i + 1) % 400 == 0:
                print(f"  {i+1}/{len(ids)}")
        A = pd.Series(vals, index=ids, name="polar_frac_out")
        A.rename_axis("tcr_id").to_frame().to_csv(CACHE)
        print(f"wrote {CACHE.name}")

    B = pd.read_csv(ROOT / "data/interim/interface_energetics.csv"
                    ).set_index("tcr_id")["polar_frac"]
    A = A.reindex(test["tcr_id"])
    B = B.reindex(test["tcr_id"])
    both = A.notna() & B.notna()
    print(f"\n{int(both.sum())} complexes with both runs")

    def residual_auc(s):
        out = pd.Series(index=s.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            ids = [t for t in g["tcr_id"] if t in s.index and pd.notna(s.get(t))]
            if len(ids) < 10:
                continue
            x = dr.reindex(ids).rank(pct=True).to_numpy()
            y = s.reindex(ids).rank(pct=True).to_numpy()
            r = y - np.polyval(np.polyfit(x, y, 1), x)
            out.loc[ids] = r
            leak.append(abs(spearmanr(r, x).statistic))
        m = hyp.mean_auc(labels, out)
        return max(m, 1 - m), float(np.mean(leak))

    def blend(s):
        sv = s if hyp.mean_auc(labels, s) >= 0.5 else -s
        sr = rank(sv)
        per_w = {0.0: hyp.per_epitope(labels, dr)}
        for w in (0.05, 0.10, 0.15, 0.20):
            per_w[w] = hyp.per_epitope(labels, (1 - w) * dr + w * sr)
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            pk = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
            held.append(per_w[pk][e])
        return float(np.mean(held))

    # ---------- 1. how reproducible is polar_frac itself? ----------
    print("\n=== 1. reproducibility of polar_frac across the two runs ===")
    r_all = float(spearmanr(A[both], B[both]).statistic)
    print(f"  pooled Spearman(run out, run out5) = {r_all:+.3f}")
    print(f"  {'epitope':<12}{'split':<9}{'rho':>8}{'sd(out)':>10}{'sd(out5)':>10}"
          f"{'sd(diff)':>10}")
    rep = []
    for e, g in test.groupby("epitope_aa"):
        ids = [t for t in g["tcr_id"] if both.get(t, False)]
        r = float(spearmanr(A.reindex(ids), B.reindex(ids)).statistic)
        d = (A.reindex(ids) - B.reindex(ids))
        print(f"  {e:<12}{g['split'].mode().iloc[0]:<9}{r:+8.3f}"
              f"{A.reindex(ids).std():10.4f}{B.reindex(ids).std():10.4f}"
              f"{d.std():10.4f}")
        rep.append([e, g["split"].mode().iloc[0], r, float(A.reindex(ids).std()),
                    float(B.reindex(ids).std()), float(d.std())])

    # ---------- 2. does averaging the replicates sharpen it? ----------
    print("\n=== 2. one run, the other run, and their average ===")
    print(f"  {'variant':<28}{'macro':>8}{'resid on D':>12}{'leak':>8}"
          f"{'nested in D':>13}")
    rows = []
    for nm, s in (("run out5 (E104's)", B), ("run out", A),
                  ("mean of the two", (A + B) / 2)):
        s2 = s.copy()
        s2 = s2.fillna(s2.median())
        ra, lk = residual_auc(s2)
        m = hyp.mean_auc(labels, s2)
        nb = blend(s2)
        print(f"  {nm:<28}{max(m,1-m):8.4f}{ra:12.4f}{lk:8.3f}"
              f"{nb:9.4f} ({nb-base:+.4f})")
        rows.append([nm, max(m, 1 - m), ra, lk, nb, nb - base])

    # ---------- 3. does disagreement mark unreliable values? ----------
    print("\n=== 3. disagreement as a reliability weight ===")
    dis = (A - B).abs()
    dis = dis.fillna(dis.median())
    mean_pf = ((A + B) / 2).fillna(((A + B) / 2).median())
    print(f"  disagreement |out - out5|: median {dis.median():.4f}, "
          f"p90 {dis.quantile(0.9):.4f}")
    ra_d, _ = residual_auc(dis)
    print(f"  disagreement as a feature by itself: residual-on-D {ra_d:.4f}")
    shr = []
    for lam in (0.0, 0.5, 1.0, 2.0, 4.0):
        w = 1.0 / (1.0 + lam * dis / (dis.median() + 1e-12))
        s = pd.Series(index=mean_pf.index, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            ids = list(g["tcr_id"])
            v = mean_pf.reindex(ids)
            ww = w.reindex(ids)
            mu = float(np.average(v, weights=ww))
            s.loc[ids] = mu + ww * (v - mu)      # shrink the unreliable to the mean
        ra, _ = residual_auc(s)
        nb = blend(s)
        print(f"  shrinkage lambda={lam:<4} residual {ra:.4f}   "
              f"nested {nb:.4f} ({nb-base:+.4f})")
        shr.append([lam, ra, nb, nb - base])

    lines = ["# E115 — An independent structural replicate as a noise model", "",
             "Generated by `scripts/104_replicate_structures.py`. Two complete",
             "Boltz runs of all 2,636 complexes exist; every structural result so",
             "far used only `out5`. On 25 sampled complexes the TCR C-alpha RMSD",
             "between runs, after superposing on MHC and peptide, has median",
             "**1.19 A** and reaches 5.76 A with none identical.", "",
             "## 1. How reproducible is `polar_frac`?", "",
             f"Pooled Spearman between runs: **{r_all:+.3f}**.", "",
             "| epitope | split | ρ | sd(out) | sd(out5) | sd(difference) |",
             "|---|---|---:|---:|---:|---:|"]
    for r in rep:
        lines.append(f"| {r[0]} | {r[1]} | **{r[2]:+.3f}** | {r[3]:.4f} "
                     f"| {r[4]:.4f} | {r[5]:.4f} |")
    lines += ["", "## 2. One run, the other, and their average", "",
              f"D alone {base:.4f}; floor for residual-on-D is 0.5098.", "",
              "| variant | macro | residual on D | leak | nested in D |",
              "|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.4f} | **{r[2]:.4f}** | {r[3]:.3f} "
                     f"| {r[4]:.4f} ({r[5]:+.4f}) |")
    lines += ["", "## 3. Disagreement as a reliability weight", "",
              f"`|out − out5|` as a feature by itself: residual-on-D "
              f"**{ra_d:.4f}**.", "",
              "| shrinkage λ | residual on D | nested in D |",
              "|---:|---:|---:|"]
    for r in shr:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {r[2]:.4f} ({r[3]:+.4f}) |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
