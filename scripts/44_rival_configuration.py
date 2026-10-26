#!/usr/bin/env python
"""Why is the other configuration ahead on public, and does that hold on private?

A parallel effort reports public macro 0.8552 against our 0.8498, from two
components rank-averaged at equal weight:

    pool-expanded neighbourhood density, top-k = 8
    whole-V-domain alignment distance, top-k = 8

We have measured both independently and rejected both, so the question is not
whether they help the public score -- E27 and E25 say they do -- but whether
the public advantage survives the only transfer that matters.

Everything here is rebuilt from our own implementations: the pool machinery
from `scripts/33`, the V-domain distances cached by `scripts/22`, our metric
and our downsampling protocol. Nothing is ported.

Three questions, in order:

  1. Does our harness reproduce the reported public advantage? If not, the
     comparison is between different measurements and nothing follows.
  2. What does the advantage look like after the private-realisable
     downsampling E27 applied to our own configuration -- each public pool cut
     to the private epitopes' own distinct-clonotype multipliers?
  3. Under E47's published anchor, which configuration is expected to score
     higher on the private split?

Run: python scripts/44_rival_configuration.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_tcrdist_baseline, load_test  # noqa: E402

REPORT = ROOT / "reports/21_rival_configuration.md"
TOPK = 8
REPEATS = 20
SEED = 0

# E47: tcrdist3's published per-epitope AUROC on the private epitopes.
PRIVATE_BANDS = {"CINGVCWTV": (0.80, 0.83), "LLWNGPMAV": (0.80, 0.83),
                 "LTDEMIAQY": (0.66, 0.68)}
# E27: distinct-clonotype pool multipliers the private epitopes can actually get.
PRIVATE_MULTIPLIERS = {"LTDEMIAQY": 1.038, "CINGVCWTV": 1.045, "LLWNGPMAV": 1.834}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def rank_within(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    return frame.assign(
        r=frame.groupby("epitope_aa")["v"].rank(pct=True)).set_index("tcr_id")["r"]


def mean_auc(labels, series) -> float:
    return float(np.mean([roc_auc_score(b["pos"], b["tcr_id"].map(series))
                          for _, b in labels.groupby("epitope_aa")]))


def topk_density(distances: np.ndarray, k: int = TOPK) -> np.ndarray:
    d = distances.astype(float).copy()
    np.fill_diagonal(d, np.inf)
    k = min(k, d.shape[1] - 1)
    return -np.sort(d, axis=1)[:, :k].mean(axis=1)


def vdomain_topk(test: pd.DataFrame) -> pd.Series:
    """Their m50, from the distances scripts/22 already cached."""
    out = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/fullseq/{epitope}.npz", allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [position[t] for t in cohort["tcr_id"]]
        d = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
        out.loc[cohort["tcr_id"]] = topk_density(d)
    return out


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    pool_mod = load_module("pool", "scripts/33_pool_negative_control.py")
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")

    import sceptr
    sceptr.disable_hardware_acceleration()
    rng = np.random.default_rng(SEED)

    print("embedding cohorts and the same-epitope pool ...")
    pool = pool_mod.load_pool()
    query, real_extra = {}, {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        frame = pd.DataFrame({
            "cdr3a": cohort["cdr3_alpha_aa"], "cdr3b": cohort["cdr3_beta_aa"],
            "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
            "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"], "ep": epitope})
        query[epitope] = (cohort, pool_mod.embed(frame))
        keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
        same = pool[pool["ep"] == epitope]
        real_extra[epitope] = same[
            [(a, b) not in keys for a, b in zip(same["cdr3a"], same["cdr3b"])]
        ].drop_duplicates(["cdr3a", "cdr3b"])

    def pooled_density(extra_by_epitope) -> pd.Series:
        out = pd.Series(index=test["tcr_id"], dtype=float)
        for epitope, (cohort, vectors) in query.items():
            extra = extra_by_epitope.get(epitope)
            reference = vectors if extra is None or extra.empty else np.vstack(
                [vectors, pool_mod.embed(extra)])
            d = 1.0 - vectors @ reference.T
            d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
            out.loc[cohort["tcr_id"]] = -np.sort(d, axis=1)[:, :TOPK].mean(axis=1)
        return out

    print("\n=== 1. does our harness reproduce the reported advantage? ===")
    no_pool = pooled_density({e: None for e in query})
    full_pool = pooled_density(real_extra)
    vdom = vdomain_topk(test)

    rivals = {
        "m30 pool-expanded top-8": full_pool,
        "m50 V-domain top-8": vdom,
        "m30 + m50 (their S04)": None,
        "m30 only, no pool (floor)": no_pool,
    }
    combo = (rank_within(full_pool, test) + rank_within(vdom, test)) / 2
    rivals["m30 + m50 (their S04)"] = combo

    ours = rank_within(decomp.submitted_score(test), test)
    baseline = load_tcrdist_baseline()
    column = "prediction" if "prediction" in baseline else baseline.columns[1]
    tcr_series = baseline.set_index("tcr_id")[column]
    # mean_auc maps by tcr_id itself; passing a pre-mapped series maps twice.
    tcrdist_macro = mean_auc(labels, tcr_series)

    rows = []
    for name, series in rivals.items():
        rows.append({"configuration": name, "public": mean_auc(labels, series)})
        print(f"  {name:<28} {rows[-1]['public']:.4f}")
    ours_public = mean_auc(labels, ours)
    print(f"  {'our frozen entry':<28} {ours_public:.4f}")
    print(f"  {'tcrdist3 baseline':<28} {tcrdist_macro:.4f}")
    their_public = rows[2]["public"]
    print(f"\n  reported 0.8552, we rebuild {their_public:.4f} "
          f"(gap {their_public - 0.8552:+.4f})")
    print(f"  their advantage over us on public: {their_public - ours_public:+.4f}")

    print("\n=== 2. the advantage after private-realisable downsampling ===")
    print("  Each public pool is cut to the private epitopes' own")
    print("  distinct-clonotype multipliers, which is how E27 priced our own")
    print("  pool gain at +0.0004 rather than +0.0112.\n")
    discounted = []
    for private_epitope, multiplier in PRIVATE_MULTIPLIERS.items():
        scores = []
        for _ in range(REPEATS):
            extra = {}
            for epitope, (cohort, _) in query.items():
                available = real_extra[epitope]
                take = int(round((multiplier - 1.0) * len(cohort)))
                extra[epitope] = (available.sample(min(take, len(available)),
                                                   random_state=int(rng.integers(1 << 30)))
                                  if take > 0 and len(available) else available.iloc[:0])
            thin = pooled_density(extra)
            scores.append(mean_auc(labels, (rank_within(thin, test)
                                            + rank_within(vdom, test)) / 2))
        discounted.append({"private epitope": private_epitope,
                           "multiplier": multiplier,
                           "their config": float(np.mean(scores))})
        print(f"  {private_epitope:<11} x{multiplier:.3f}   "
              f"their config {np.mean(scores):.4f}")
    discounted = pd.DataFrame(discounted)
    their_realisable = float(discounted["their config"].mean())
    print(f"\n  their config at private pool density: {their_realisable:.4f}")
    print(f"  at full public pool:                  {their_public:.4f}")
    print(f"  pool advantage that cannot transfer:  "
          f"{their_public - their_realisable:+.4f}")
    print(f"  our frozen entry (no pool at all):    {ours_public:.4f}")
    print(f"  advantage that survives:              "
          f"{their_realisable - ours_public:+.4f}")

    print("\n=== 3. under E47's published anchor ===")
    lo = float(np.mean([v[0] for v in PRIVATE_BANDS.values()]))
    hi = float(np.mean([v[1] for v in PRIVATE_BANDS.values()]))
    print(f"  tcrdist3 private macro, published: {lo:.3f}-{hi:.3f}")
    table = []
    for name, public_score, transferable in (
        ("our frozen entry", ours_public, ours_public - tcrdist_macro),
        ("their S04, public margin", their_public, their_public - tcrdist_macro),
        ("their S04, realisable margin", their_public,
         their_realisable - tcrdist_macro),
    ):
        table.append({"configuration": name, "public": public_score,
                      "margin used": transferable,
                      "private lo": lo + transferable, "private hi": hi + transferable,
                      "private mid": (lo + hi) / 2 + transferable})
        print(f"  {name:<30} margin {transferable:+.4f}   "
              f"private {(lo + hi) / 2 + transferable:.3f}")
    anchor = pd.DataFrame(table)

    print("\n  The public gap is real. Most of it is pool, and E27 measured that")
    print("  the private epitopes cannot receive it: two of three can add about")
    print("  4% more distinct clonotypes.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Rebuilding the rival configuration in our harness\n\n"
        "## Public scores\n\n"
        + pd.DataFrame(rows + [{"configuration": "our frozen entry",
                                "public": ours_public},
                               {"configuration": "tcrdist3 baseline",
                                "public": tcrdist_macro}]).round(4).to_markdown(index=False)
        + "\n\n## After private-realisable downsampling\n\n"
        + discounted.round(4).to_markdown(index=False)
        + "\n\n## Under E47's published anchor\n\n"
        + anchor.round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
