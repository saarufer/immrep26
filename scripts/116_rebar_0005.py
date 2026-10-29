#!/usr/bin/env python
"""E132 -- Re-evaluate every closed candidate against a +0.005 bar.

The +0.02 bar was pre-registered, and E131 showed what it really demands: a
measured +0.02 needs a TRUE effect of +0.0309 to be cleared 80% of the time,
which is ~19x the largest effect any channel here has shown. Lowering the bar
to +0.005 is therefore a reasonable request, but a bar is only meaningful
together with its two error rates, so both are measured first:

  power  the TRUE effect needed to clear a measured +0.005 at 80%
  size   how often a feature with NO signal clears +0.005

The size calculation has to include the search each experiment actually ran.
E104 tested nine interface features and reported the best; E99 tested six
hubness variants; E112 tested eleven pMHC-side features. A bar applied to the
best of k is not the same bar applied once, so both are reported.

Then the re-evaluation itself is two-layered, because lowering an effect-size
bar must not resurrect a candidate that failed a **validity** gate:

  layer 1  typicality gate, private variance gate, residual leak < 0.1,
           and Q3 independence from the labelling process
  layer 2  nested LOEO vs +0.005

A candidate that fails layer 1 is not readmitted at any bar.

Run: .venv/bin/python scripts/116_rebar_0005.py
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

REPORT = ROOT / "reports/94_rebar_0005.md"
WS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
NEW_BAR = 0.005
OLD_BAR = 0.02
DRAWS = 400
KS = (1, 6, 9, 11)
SIGMAS = [None, 8.0, 5.0, 3.5, 2.8, 2.2, 1.8, 1.5]   # None = pure noise
SEED = 0


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
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)
    pub = sorted(scored["epitope_aa"].unique())
    y = scored.set_index("tcr_id")["pos"]
    base_per = hyp.per_epitope(labels, dr)
    base = float(np.mean([base_per[e] for e in pub]))
    ids = test["tcr_id"]
    rng = np.random.default_rng(SEED)

    def feature(sigma, rg):
        v = pd.Series(rg.normal(0, 1.0, len(ids)), index=ids)
        if sigma is not None:
            v = v * sigma + y.reindex(ids).fillna(y.mean())
        out = pd.Series(index=ids, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            out.loc[g["tcr_id"]] = v.reindex(g["tcr_id"]).rank(pct=True).to_numpy()
        return out

    def nested_of(f):
        sr = hyp.rank_within(f, test)
        per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
        held = []
        for e in pub:
            trn = [x for x in pub if x != e]
            pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
            held.append(per[pick][e])
        return float(np.mean(held)) - base, per

    print("=== 1. the +0.005 bar, characterised before it is used ===")
    print(f"  {'sigma':>7}{'true effect':>13}{'nested mean':>13}"
          f"{'P(>+0.005)':>12}{'P(>+0.02)':>11}")
    curve = []
    for sg in SIGMAS:
        acc = {w: {e: [] for e in pub} for w in WS}
        ns = []
        for _ in range(DRAWS):
            n, per = nested_of(feature(sg, rng))
            ns.append(n)
            for w in WS:
                for e in pub:
                    acc[w][e].append(per[w][e])
        te = max(float(np.mean([np.mean(acc[w][e]) for e in pub])) - base
                 for w in WS)
        ns = np.array(ns)
        curve.append([sg, te, ns.mean(), float((ns > NEW_BAR).mean()),
                      float((ns > OLD_BAR).mean()), ns])
        lab = "noise" if sg is None else f"{sg:.1f}"
        print(f"  {lab:>7}{te:+13.4f}{ns.mean():+13.4f}"
              f"{curve[-1][3]:12.2f}{curve[-1][4]:11.2f}")

    null = curve[0][5]
    print(f"\n  NULL (no signal): nested mean {null.mean():+.4f} "
          f"sd {null.std():.4f}  max {null.max():+.4f}")
    print(f"  one random feature clears +0.005: {float((null>NEW_BAR).mean())*100:.1f}%"
          f"   clears +0.02: {float((null>OLD_BAR).mean())*100:.1f}%")
    print("\n  and the search each experiment actually ran -- best of k features:")
    for k in KS:
        p = 1 - (1 - float((null > NEW_BAR).mean())) ** k
        p2 = 1 - (1 - float((null > OLD_BAR).mean())) ** k
        print(f"    best of {k:2d}:  P(clears +0.005) {p*100:5.1f}%"
              f"    P(clears +0.02) {p2*100:5.1f}%")

    ok = [c for c in curve if c[3] >= 0.80]
    mde = min(c[1] for c in ok) if ok else np.nan
    print(f"\n  TRUE effect needed to clear +0.005 at 80% power: {mde:+.4f}")
    print(f"  (for comparison, E131: +0.0309 for the +0.02 bar)")

    # ---------------- 2. the re-evaluation ----------------
    # nested value as logged, plus layer-1 validity status
    C = [
        # name, nested, layer-1 verdict, note
        ("E105 full stack (signed)",       +0.0103, "FAIL",
         "E106 sign-selection leak; honest null p 0.200; E118 16-fold -0.0017"),
        ("E116 6-structure average",       +0.0047, "pass",
         "residual flat -0.0005; gain is 1.04 null sd; no new information"),
        ("E98 joint candidate (best)",     +0.0020, "pass",
         "best of 3 orthogonal components; 98% additive"),
        ("E112 pMHC-side (best of 11)",    +0.0019, "FAIL",
         "pep_frac_covered residual leak 0.806; all 11 re-audited"),
        ("E100 TwoNN LID (best)",          +0.0019, "pass", "same distance matrix"),
        ("E104 polar_frac",                +0.0016, "pass",
         "ACCEPTED already; the one live channel"),
        ("E113 V-gene prior (best)",       +0.0016, "FAIL",
         "Q3: aggregates the competition's own labels; 54% pair coverage"),
        ("E99 hubness (best of 6)",        +0.0009, "pass", "same distance matrix"),
        ("E101 cross-epitope LID",         +0.0008, "pass", "same distance matrix"),
        ("E97 per-epitope bandwidth",      +0.0007, "pass", "reparameterisation"),
        ("E96 lenfreq",                    +0.0005, "pass", "length already in SCEPTR"),
        ("E103 germline docking",          +0.0000, "pass", "reducible to m50"),
        ("E122 Atlas pool",                +0.0000, "pass", "CIN gets 0 new rows"),
        ("E126 Minervina margin",          -0.0003, "pass",
         "residual 0.5609 ties the record; still negative nested"),
        ("E118 the stack, 16 folds",       -0.0017, "pass", "p 0.683"),
        ("E121 Hi-TpH",                    np.nan,  "FAIL",
         "Q3: 100% training overlap, 872 inverted labels"),
    ]
    print(f"\n=== 2. re-evaluation at the +0.005 bar ===")
    print(f"  {'candidate':<30}{'nested':>9}{'layer 1':>9}{'+0.02':>7}{'+0.005':>8}")
    newpass = []
    for nm, v, l1, note in C:
        p_old = "pass" if (not np.isnan(v) and v > OLD_BAR) else "no"
        if l1 == "FAIL":
            p_new = "blocked"
        elif np.isnan(v):
            p_new = "n/a"
        else:
            p_new = "PASS" if v > NEW_BAR else "no"
        if p_new == "PASS":
            newpass.append((nm, v, note))
        vs = "   n/a" if np.isnan(v) else f"{v:+.4f}"
        print(f"  {nm:<30}{vs:>9}{l1:>9}{p_old:>7}{p_new:>8}")

    print(f"\n  candidates that clear +0.005 and pass layer 1: {len(newpass)}")
    for nm, v, note in newpass:
        print(f"    {nm}  {v:+.4f}")
        print(f"      -> {note}")
    print(f"\n  expected FALSE passes if 16 candidates are each the best of ~9:")
    fp1 = float((null > NEW_BAR).mean())
    print(f"    per-candidate FPR (best of 9) "
          f"{(1-(1-fp1)**9)*100:.1f}%  x 16 = "
          f"{16*(1-(1-fp1)**9):.1f} expected false passes")

    lines = ["# E132 — Re-evaluating every closed candidate at a +0.005 bar", "",
             "Generated by `scripts/116_rebar_0005.py`.", "",
             "## The bar, characterised first", "",
             "| σ | true effect | nested mean | P(>+0.005) | P(>+0.02) |",
             "|---|---:|---:|---:|---:|"]
    for sg, te, nm_, p5, p20, _ in curve:
        lab = "**noise**" if sg is None else f"{sg:.1f}"
        lines.append(f"| {lab} | {te:+.4f} | {nm_:+.4f} | {p5:.2f} | {p20:.2f} |")
    lines += ["", f"- null: mean {null.mean():+.4f}, sd {null.std():.4f}, "
              f"max {null.max():+.4f}",
              f"- one random feature clears +0.005 **{fp1*100:.1f}%** of the time",
              f"- best of 9 random features clears it "
              f"**{(1-(1-fp1)**9)*100:.1f}%** of the time",
              f"- TRUE effect needed to clear +0.005 at 80% power: **{mde:+.4f}** "
              f"(vs +0.0309 for the +0.02 bar)", "",
              "## Re-evaluation", "",
              "| candidate | nested | layer 1 (validity) | +0.02 | +0.005 |",
              "|---|---:|---|---|---|"]
    for nm, v, l1, note in C:
        vs = "n/a" if np.isnan(v) else f"{v:+.4f}"
        p_old = "pass" if (not np.isnan(v) and v > OLD_BAR) else "no"
        p_new = ("blocked" if l1 == "FAIL" else
                 "n/a" if np.isnan(v) else
                 "**PASS**" if v > NEW_BAR else "no")
        lines.append(f"| {nm} | {vs} | {l1} | {p_old} | {p_new} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
