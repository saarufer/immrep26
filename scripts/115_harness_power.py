#!/usr/bin/env python
"""E131 -- Minimum detectable effect of this project's own evaluation harness.

Prompted by the right question: are hypotheses being rejected because the data
cannot resolve them, rather than because they are false?

The first thing to separate is *which* data. The private split is never an
evaluation set -- there are no labels there -- so no rejection in this project
was ever decided by private scarcity. Of sixteen closed candidates, thirteen
were decided by nested LOEO or the residual gate on the **five public
epitopes**, one by circularity (E121), one by an external panel's n (E128/E129),
and exactly one -- E122's LTDEMIAQY pool bet -- was genuinely unmeasurable
because the private cohort could not be checked.

So the question becomes: **how small a true effect can nested LOEO over five
epitopes even see?** That is a property of the harness, not of any candidate,
and it is measurable by injecting a signal of known strength.

Method. Build a synthetic feature with a *genuinely transferable* effect:

    f = y + Normal(0, sigma)   per epitope, then ranked within epitope

Because it is built from the labels, it transfers perfectly across epitopes --
every cohort gets exactly the same quality of signal. That makes the resulting
MDE a **lower bound**: any real feature, whose usefulness varies by cohort,
is strictly harder to detect than this.

The true (population) effect of a given sigma is estimated by averaging the
per-epitope gain over many draws before selecting the weight, so selection
noise does not inflate it. Then the actual gate is run on each draw -- the same
weight grid, the same leave-one-epitope-out weight selection as every candidate
faced -- and the detection rate recorded.

Reported: the true effect at which nested LOEO exceeds zero 80% of the time,
and the true effect needed to clear the pre-registered +0.02 bar.

Run: .venv/bin/python scripts/115_harness_power.py
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

REPORT = ROOT / "reports/93_harness_power.md"
WS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
SIGMAS = [8.0, 5.0, 3.5, 2.5, 2.0, 1.6, 1.3, 1.1, 0.9, 0.7, 0.5]
DRAWS = 300
SEED = 0
E118_NULL_SD = 0.0045      # measured null sd of the 16-fold estimator
BAR = 0.02                 # pre-registered inclusion bar


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)
    pub = sorted(scored["epitope_aa"].unique())
    y = scored.set_index("tcr_id")["pos"]
    base_per = hyp.per_epitope(labels, dr)
    base = float(np.mean([base_per[e] for e in pub]))
    print(f"=== D on the five public epitopes: {base:.4f} ===")
    for e in pub:
        print(f"  {e:<12}{base_per[e]:.4f}  n={int((scored.epitope_aa==e).sum())}")
    print(f"\n  between-epitope sd of D's own AUC: "
          f"{np.std([base_per[e] for e in pub]):.4f}  <- this is what the")
    print("  5-fold mean has to see through, and it does not shrink with rows")

    rng = np.random.default_rng(SEED)
    ids = test["tcr_id"]
    ep_of = test.set_index("tcr_id")["epitope_aa"]

    def make_feature(sigma, rg):
        """Label plus noise, ranked within epitope. Perfectly transferable."""
        v = pd.Series(rg.normal(0, sigma, len(ids)), index=ids)
        v = v + y.reindex(ids).fillna(y.mean())
        out = pd.Series(index=ids, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            s = v.reindex(g["tcr_id"])
            out.loc[g["tcr_id"]] = s.rank(pct=True).to_numpy()
        return out

    print(f"\n=== power curve, {DRAWS} draws per sigma ===")
    print(f"  {'sigma':>6}{'TRUE effect':>13}{'nested mean':>13}{'nested sd':>11}"
          f"{'P(>0)':>8}{'P(>1sd)':>9}{'P(>+0.02)':>11}")
    rows = []
    for sg in SIGMAS:
        per_w_acc = {w: {e: [] for e in pub} for w in WS}
        nested = []
        for _ in range(DRAWS):
            f = make_feature(sg, rng)
            sr = hyp.rank_within(f, test)
            per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
            for w in WS:
                for e in pub:
                    per_w_acc[w][e].append(per[w][e])
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pick][e])
            nested.append(float(np.mean(held)) - base)
        # TRUE effect: average per-epitope gain across draws FIRST, then pick w
        exp_gain = {w: float(np.mean([np.mean(per_w_acc[w][e]) for e in pub])) - base
                    for w in WS}
        true_eff = max(exp_gain.values())
        nested = np.array(nested)
        rows.append([sg, true_eff, nested.mean(), nested.std(),
                     float((nested > 0).mean()),
                     float((nested > E118_NULL_SD).mean()),
                     float((nested > BAR).mean())])
        print(f"  {sg:6.1f}{true_eff:+13.4f}{nested.mean():+13.4f}"
              f"{nested.std():11.4f}{rows[-1][4]:8.2f}{rows[-1][5]:9.2f}"
              f"{rows[-1][6]:11.2f}")

    R = np.array(rows)
    print("\n=== minimum detectable effect ===")
    for col, name, tgt in ((4, "nested > 0", 0.80), (5, "nested > 1 null sd", 0.80),
                           (6, "nested > the +0.02 bar", 0.80)):
        ok = R[R[:, col] >= tgt]
        if len(ok):
            mde = ok[:, 1].min()
            print(f"  {name:<24}80% power at a TRUE effect of {mde:+.4f}")
        else:
            print(f"  {name:<24}never reaches 80% power in this sweep "
                  f"(max {R[:,col].max():.2f})")
    # what the closed candidates actually measured
    print("\n=== where the closed candidates sat on this curve ===")
    cands = [("E104 polar_frac (accepted)", 0.0016), ("E100 TwoNN", 0.0019),
             ("E97 bandwidth", 0.0007), ("E96 lenfreq", 0.0005),
             ("E126 Minervina margin", -0.0003), ("E118 the whole stack", -0.0017),
             ("E116 6-structure average", 0.0047)]
    for nm, v in cands:
        # probability the harness would report >0 if the TRUE effect were v
        j = int(np.argmin(np.abs(R[:, 1] - abs(v))))
        print(f"  {nm:<28}nested {v:+.4f}   if TRUE effect were "
              f"{abs(v):+.4f}, P(nested>0) ~ {R[j,4]:.2f}")

    lines = ["# E131 — Minimum detectable effect of the evaluation harness", "",
             "Generated by `scripts/115_harness_power.py`.", "",
             "The private split is never an evaluation set, so no rejection was",
             "decided by private scarcity. Thirteen of sixteen closed candidates",
             "were decided by nested LOEO or the residual gate on the five",
             "public epitopes. This measures what that procedure can see.", "",
             f"D on the five public epitopes: **{base:.4f}**; between-epitope sd",
             f"of D's own per-epitope AUC: **{np.std([base_per[e] for e in pub]):.4f}**.",
             "", "The injected feature is built from the labels, so it transfers",
             "perfectly — the MDE below is a **lower bound** on the真 difficulty.",
             "", "| σ | true effect | nested mean | nested sd | P(>0) | P(>1 null sd) | P(>+0.02) |",
             "|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r[0]:.1f} | {r[1]:+.4f} | {r[2]:+.4f} | {r[3]:.4f} "
                     f"| {r[4]:.2f} | {r[5]:.2f} | {r[6]:.2f} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
