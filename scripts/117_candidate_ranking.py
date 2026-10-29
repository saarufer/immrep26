#!/usr/bin/env python
"""E134 -- All submission candidates ranked by expected private score.

Two things have to be kept apart, because only one of them is estimable.

**The private LEVEL is not identifiable.** The class-composition projection
that produced the often-quoted 0.799 fails its own leave-one-epitope-out test:
MAE 0.0833, worst error 0.2129 (GLCTLVAML predicted 0.954 against an actual
0.741), and its sensitivity surface spans 0.829-0.857 across an unobservable
parameter. No per-candidate private level is defensible.

**The DIFFERENCES between candidates are estimable**, and since the level is
shared by the whole family, the differences alone determine the ranking.

For every candidate the public macro is exact -- E47 established harness ==
leaderboard. The expected private difference is then built as follows:

  unfitted candidates   architecture ablations (F, density4, density3, the
                        single-component baseline). Their difference from D is
                        a real architectural difference, not a fitted gain, so
                        it is carried across unshrunk.
  fitted candidates     the gain was chosen on the five public epitopes, so the
                        public difference is in-sample. Each one's weights are
                        REFIT leave-one-epitope-out here and the held-out gain
                        reported -- the same procedure every feature faced.
  already measured      S4 (E118, 16 pseudo-epitope folds, -0.0017) and
                        S7 (E122, +0.0000) keep their measured values.

Run: .venv/bin/python scripts/117_candidate_ranking.py
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

REPORT = ROOT / "reports/95_candidate_ranking.md"
GRID_W = [0.0, 0.05, 0.10, 0.15, 0.20]


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    pub = sorted(scored["epitope_aa"].unique())

    def rank(s):
        return hyp.rank_within(s.reindex(test["tcr_id"]), test)

    sub = {n: rank(pd.read_csv(ROOT / f"submissions/{f}"
                               ).set_index("tcr_id")["prediction"])
           for n, f in (("D", "d_density3_m30pool_m50_struct.csv"),
                        ("F", "f_density3_m30_m50_struct_nopool.csv"))}
    per = {n: hyp.per_epitope(labels, s) for n, s in sub.items()}

    # the three components S1/S2/S3/S5 are built from
    Eg = pd.read_csv(ROOT / "data/interim/interface_energetics.csv"
                     ).set_index("tcr_id")
    pf = Eg["polar_frac"].reindex(test["tcr_id"])
    pfr = rank(-pf.fillna(pf.median()))
    lf = {}
    for e, g in test.groupby("epitope_aa"):
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        for t, v in zip(g["tcr_id"],
                        (la.map(la.value_counts(normalize=True))
                         * lb.map(lb.value_counts(normalize=True))).to_numpy()):
            lf[t] = v
    lfr = rank(pd.Series(lf))
    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}
    tn, EPS = {}, 1e-12
    for e, g in test.groupby("epitope_aa"):
        idx = np.array([po[t] for t in g["tcr_id"]])
        Dd = 1.0 - W[idx] @ W[idx].T
        np.fill_diagonal(Dd, np.inf)
        sd = np.sort(Dd, axis=1)
        r1, r2 = np.maximum(sd[:, 0], EPS), np.maximum(sd[:, 1], EPS)
        tn.update(dict(zip(g["tcr_id"],
                           1.0 / np.log(np.maximum(r2 / r1, 1 + EPS)))))
    tnr = rank(-pd.Series(tn))
    COMP = {"polar_frac": pfr, "lenfreq": lfr, "TwoNN": tnr}

    def mix(base, weights):
        s = base
        for nm, w in weights:
            if w:
                s = (1 - w) * s + w * COMP[nm]
        return s

    def nested_refit(base_name, comps):
        """Refit the component weights leave-one-epitope-out."""
        base = sub[base_name]
        bper = per[base_name]
        grids = [GRID_W] * len(comps)
        combos = [()]
        for g in grids:
            combos = [c + (w,) for c in combos for w in g]
        cache = {}
        for cb in combos:
            s = mix(base, list(zip(comps, cb)))
            cache[cb] = hyp.per_epitope(labels, s)
        held, picks = [], []
        for e in pub:
            trn = [x for x in pub if x != e]
            best = max(combos, key=lambda cb: np.mean([cache[cb][t] for t in trn]))
            held.append(cache[best][e] - bper[e])
            picks.append((e, best))
        return float(np.mean(held)), held, picks, cache, combos

    print("=== refitting the fitted candidates leave-one-epitope-out ===")
    results = {}
    for nm, bn, comps in (("S1_D_plus_polar", "D", ["polar_frac"]),
                          ("S2_D_stack5", "D", ["lenfreq", "TwoNN", "polar_frac"]),
                          ("S3_F_stack5", "F", ["lenfreq", "TwoNN", "polar_frac"]),
                          ("S5_D_ttd_tuned", "D", ["lenfreq", "TwoNN", "polar_frac"])):
        oof, held, picks, cache, combos = nested_refit(bn, comps)
        # in-sample optimum, for the contrast
        ins = max(np.mean([cache[cb][t] for t in pub]) for cb in combos) \
            - np.mean([per[bn][t] for t in pub])
        results[nm] = (oof, ins, held, picks, bn, comps)
        print(f"  {nm:<20}base {bn}  components {comps}")
        print(f"    in-sample best  {ins:+.4f}")
        print(f"    OUT-OF-FOLD     {oof:+.4f}   per fold "
              f"{['%+.4f' % h for h in held]}")

    # S5 is the special case: its published weights came from ONE epitope
    print("\n  S5 note: the shipped S5 chose its weights on TTDPSFLGRY alone.")
    oof5, held5 = [], []
    base, bper = sub["D"], per["D"]
    comps = ["lenfreq", "TwoNN", "polar_frac"]
    combos = [()]
    for _ in comps:
        combos = [cb + (w,) for cb in combos for w in GRID_W]
    cache = {cb: hyp.per_epitope(labels, mix(base, list(zip(comps, cb))))
             for cb in combos}
    for e in pub:
        others = [x for x in pub if x != e]
        for anchor in others:
            best = max(combos, key=lambda cb: cache[cb][anchor])
            held5.append(cache[best][e] - bper[e])
    print(f"    weights from ONE other epitope, scored on a held-out one: "
          f"{np.mean(held5):+.4f} +- {np.std(held5):.4f} over {len(held5)} pairs")

    print("\n=== final ranking by expected private difference vs D ===")
    MEASURED = {"S4_D_stack9": (-0.0017, "E118, 16 pseudo-epitope folds"),
                "S7_D_atlas_pool": (0.0000, "E122, measured directly")}
    import glob, os
    allsub = {}
    for f in sorted(glob.glob(str(ROOT / "submissions/*.csv"))):
        n = os.path.basename(f)[:-4]
        allsub[n] = rank(pd.read_csv(f).set_index("tcr_id")["prediction"])
    AP = {n: hyp.per_epitope(labels, s) for n, s in allsub.items()}
    Dn = "d_density3_m30pool_m50_struct"
    METHOD = {
        "S4_D_stack9": "D restructured + 7 mixed components (only ~50% D)",
        "S2_D_stack5": "D restructured + lenfreq + TwoNN + polar_frac",
        "S5_D_ttd_tuned": "same 3 components, weights fitted on TTDPSFLGRY alone",
        "S1_D_plus_polar": "D + polar_frac at w=0.15, one component",
        "S3_F_stack5": "stack5 on the pool-free base F",
        "S7_D_atlas_pool": "D with the TCR Motif Atlas added to the pool",
        Dn: "density3 + m30pool + m50 + structure (the audited entry)",
        "f_density3_m30_m50_struct_nopool": "D with the pool term removed",
        "density4_struct03": "4-space density + structure",
        "density3_struct03": "3-space density + structure (the frozen entry)",
        "density_sceptr_default_sigma0.1": "single SCEPTR density, sigma 0.1",
    }
    rows = []
    for n in allsub:
        g = np.array([AP[n][e] - AP[Dn][e] for e in pub])
        pubm = float(np.mean([AP[n][e] for e in pub]))
        if n in MEASURED:
            exp, basis = MEASURED[n]
        elif n in results:
            exp, basis = results[n][0], "weights refit leave-one-epitope-out here"
        elif n == Dn:
            exp, basis = 0.0, "reference, 14-check audit (E123)"
        else:
            exp, basis = g.mean(), "architecture difference, not fitted"
        rows.append([n, pubm, g.mean(), g.min(), int((g < 0).sum()), exp, basis])
    rows.sort(key=lambda r: -r[5])
    print(f"  {'#':<3}{'candidate':<34}{'public':>8}{'pub vs D':>10}"
          f"{'worst ep':>10}{'neg':>4}{'EXP priv vs D':>15}")
    for i, r in enumerate(rows, 1):
        print(f"  {i:<3}{r[0]:<34}{r[1]:8.4f}{r[2]:+10.4f}{r[3]:+10.4f}"
              f"{r[4]:4d}{r[5]:+15.4f}")

    lines = ["# E134 — All candidates ranked by expected private score", "",
             "Generated by `scripts/117_candidate_ranking.py`.", "",
             "**The private level is not identifiable.** The class-composition",
             "projection behind the quoted 0.799 fails its own leave-one-epitope-out",
             "test (MAE 0.0833, worst 0.2129) and its sensitivity surface spans",
             "0.829–0.857. Only differences are estimable, and since the level is",
             "shared the differences determine the ranking.", "",
             "| # | candidate | public (exact) | public vs D | worst epitope | loses on | **expected private vs D** | main method | basis |",
             "|---:|---|---:|---:|---:|---:|---:|---|---|"]
    for i, r in enumerate(rows, 1):
        lines.append(f"| {i} | `{r[0]}` | {r[1]:.4f} | {r[2]:+.4f} | {r[3]:+.4f} "
                     f"| {r[4]}/5 | **{r[5]:+.4f}** | {METHOD.get(r[0],'')} "
                     f"| {r[6]} |")
    lines += ["", "## In-sample against out-of-fold, for the fitted candidates", "",
              "| candidate | in-sample best | out-of-fold | collapse |",
              "|---|---:|---:|---:|"]
    for nm, (oof, ins, held, picks, bn, comps) in results.items():
        lines.append(f"| `{nm}` | {ins:+.4f} | **{oof:+.4f}** | "
                     f"{(1-oof/ins)*100 if ins else 0:.0f}% |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
