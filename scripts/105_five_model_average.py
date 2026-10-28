#!/usr/bin/env python
"""Averaging all five Boltz models, to spend the headroom E115 measured.

E115 decomposed `polar_frac` with an independent run and found single-structure
reliability **0.812**: 81% reproducible signal, 19% prediction noise. Spearman–
Brown then says averaging *m* structures gives reliability
`m·r / (1 + (m−1)·r)`, so five give **0.956** and the disattenuated residual-on-D
ceiling of 0.5650 should be almost entirely reachable.

Every structural result in this project used `model_0` alone. `out5` holds
**five** models per complex and `out` holds a sixth independent one, all
already on disk. This computes `polar_frac` on all of them and tests the
prediction quantitatively: the measured reliability of the m-model mean should
track Spearman–Brown, and the residual-on-D should climb toward 0.5650 and stop.

If it climbs and stops where predicted, the noise model is right and the
feature is as good as the structures allow. If it keeps climbing, something
other than noise averaging is happening and the model is wrong.

Run: .venv/bin/python scripts/105_five_model_average.py
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

REPORT = ROOT / "reports/82_five_model_average.md"
CACHE = ROOT / "data/interim/polar_frac_all_models.csv"
OUT5 = ROOT / "data/interim/boltz/out5/boltz_results_inputs/predictions"
OUT1 = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
NEAR = 9.0


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def polar_frac(path: Path):
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
    tcr, pmhc = np.isin(ch, ["A", "B"]), np.isin(ch, ["M", "P"])
    if not (tcr.any() and pmhc.any()):
        return np.nan
    ti = np.where(tcr)[0]
    near = np.array([len(x) > 0 for x in
                     cKDTree(xyz[pmhc]).query_ball_point(xyz[ti], r=NEAR)])
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
        F = pd.read_csv(CACHE).set_index("tcr_id")
    else:
        rows, ids = [], test["tcr_id"].tolist()
        for i, t in enumerate(ids):
            rec = {"tcr_id": t}
            for k in range(5):
                p = OUT5 / t / f"{t}_model_{k}.pdb"
                rec[f"m{k}"] = polar_frac(p) if p.exists() else np.nan
            p = OUT1 / t / f"{t}_model_0.pdb"
            rec["runA"] = polar_frac(p) if p.exists() else np.nan
            rows.append(rec)
            if (i + 1) % 300 == 0:
                print(f"  {i+1}/{len(ids)}")
        F = pd.DataFrame(rows).set_index("tcr_id")
        F.to_csv(CACHE)
        print(f"wrote {CACHE.name}")

    cols = [f"m{k}" for k in range(5)] + ["runA"]
    F = F.reindex(test["tcr_id"])
    print(f"\n{len(F)} complexes, {len(cols)} structures each, "
          f"{int(F[cols].notna().all(axis=1).sum())} complete")

    def residual_auc(s):
        out = pd.Series(index=s.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            ids = [t for t in g["tcr_id"] if t in s.index]
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

    # ---------- measured reliability of the m-model mean ----------
    print("\n=== reliability as structures are averaged ===")
    print("  measured by split-half over the five out5 models, Spearman-Brown")
    print("  corrected; the single-structure value is E115's 0.812.")
    print(f"  {'m':>2}{'measured reliability':>22}{'Spearman-Brown predicts':>26}"
          f"{'residual on D':>15}{'nested in D':>14}")
    r1 = 0.812
    rows_out = []
    for m in (1, 2, 3, 5, 6):
        use = cols[:m] if m <= 5 else cols
        s = F[use].mean(axis=1)
        s = s.fillna(s.median())
        if m == 1:
            rel_meas = r1
        else:
            h = max(1, m // 2)
            a = F[use[:h]].mean(axis=1)
            b = F[use[h:]].mean(axis=1)
            ok = a.notna() & b.notna()
            rho = float(spearmanr(a[ok], b[ok]).statistic)
            rel_meas = 2 * rho / (1 + rho)          # Spearman-Brown to full length
        pred = m * r1 / (1 + (m - 1) * r1)
        ra, _ = residual_auc(s)
        nb = blend(s)
        print(f"  {m:>2}{rel_meas:>22.3f}{pred:>26.3f}{ra:>15.4f}"
              f"{nb:>9.4f} ({nb-base:+.4f})")
        rows_out.append([m, rel_meas, pred, ra, nb, nb - base])

    ceiling = (0.5586 - 0.5) / np.sqrt(r1) + 0.5
    print(f"\n  E115's disattenuated ceiling: {ceiling:.4f}")
    best = max(rows_out, key=lambda r: r[3])
    print(f"  best achieved: {best[3]:.4f} at m={best[0]} "
          f"({best[3]/ceiling*100:.0f}% of the way to the ceiling from 0.5)")

    lines = ["# E116 — Averaging all five Boltz models", "",
             "Generated by `scripts/105_five_model_average.py`. E115 measured",
             "single-structure reliability at **0.812**; Spearman–Brown predicts",
             "`m·r / (1 + (m−1)·r)` for the m-structure mean, and the",
             "disattenuated residual-on-D ceiling is **0.5650**.", "",
             f"D alone {base:.4f}; floor for residual-on-D is 0.5098.", "",
             "| structures averaged | measured reliability | Spearman–Brown predicts | residual on D | nested in D |",
             "|---:|---:|---:|---:|---:|"]
    for r in rows_out:
        lines.append(f"| {r[0]} | {r[1]:.3f} | {r[2]:.3f} | **{r[3]:.4f}** "
                     f"| {r[4]:.4f} ({r[5]:+.4f}) |")
    lines += ["", f"E115's ceiling {ceiling:.4f}; best achieved **{best[3]:.4f}** "
              f"at m={best[0]}.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
