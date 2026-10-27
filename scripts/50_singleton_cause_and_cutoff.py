#!/usr/bin/env python
"""Why is a TCR a singleton, and is there a principled radius at all?

Two questions, both the user's, and they turn out to be the same question.

**Why a singleton?** Three causes compete, and they make opposite predictions:

  sampling          the cohort is finite and the relatives were never collected.
                    Predicts singleton status is NOT predictable from anything
                    intrinsic to the TCR.
  hard to generate   convergent recombination cannot produce relatives for a
                    sequence it rarely makes. Predicts singletons have LOWER
                    paired Pgen.
  structurally odd   the fold or interface is unusual (the user's hypothesis).
                    Predicts singletons have worse or less reproducible Boltz
                    structures.

These are separable because the target -- is this TCR a singleton -- carries no
label, so it can be predicted from features directly and measured per epitope,
which holds cohort size fixed. The Pgen *direction* is what distinguishes the
second from the third, and E7 already found higher Pgen goes with higher density,
so the prediction is sharp.

**Is the radius principled?** E52 found `nn3_dist` scores 0.7793 *within*
singletons, above kernel density's 0.7702, which suggests 0.12 is a cut through
a continuum rather than a boundary between two kinds of TCR. Three ways to check:

  1. the shape of the nearest-neighbour distance distribution, against a
     background of distances to OTHER epitopes' TCRs -- the TCRNET/ALICE
     enrichment logic, used here to choose a radius rather than as a feature
     (that use failed in E8);
  2. validation rate by nn3 decile: a threshold should show a knee;
  3. the decisive one -- sweep the radius and ask whether the BINARY clustered
     indicator at its best radius beats the CONTINUOUS nn3_dist. If no radius
     does, the binary notion is strictly worse than the variable underneath it
     and should be dropped.

Run: .venv/bin/python scripts/50_singleton_cause_and_cutoff.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import typicality  # noqa: E402

REPORT = ROOT / "reports/26_singleton_cause_and_cutoff.md"
RADIUS = 0.12
SIGMA = 0.1
RADII = np.round(np.arange(0.02, 0.401, 0.02), 3)
SEED = 0

# A caution on reading these as independent causes: they are not independent.
# ms_iptm_sd correlates with ms_iptm_mean at Spearman +0.843 and
# consist_tcr_rmsd_mean with it at +0.676, so the two quantities that look like
# "reproducibility" measures co-vary with the confidence mean rather than
# opposing it. A first reading of the output treated them as contradicting the
# confidence features; they do not, they are attenuated copies of them.
CAUSE_FEATURES = {
    "log10_pgen_paired": "hard to generate",
    "ms_iptm_mean": "structurally odd",
    "ms_iptm_sd": "structurally odd",
    "consist_tcr_rmsd_mean": "structurally odd",
    "complex_plddt": "structurally odd",
    "pae_tcr_pep_min": "structurally odd",
    "iptm": "structurally odd",
    "typicality": "chain-feature atypicality",
    "cdr3_len_total": "chain-feature atypicality",
}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def space():
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors, {t: i for i, t in enumerate(cached["tcr_id"])}


def neighbour_table(test, vectors, position) -> pd.DataFrame:
    """Within-epitope nearest-neighbour distances, and the same against the
    other epitopes' TCRs as a background. Both label-free."""
    rows = []
    all_ids = test["tcr_id"].to_numpy()
    index_of = {t: i for i, t in enumerate(all_ids)}
    matrix = vectors[[position[t] for t in all_ids]]
    epitope_of = test.set_index("tcr_id")["epitope_aa"].to_dict()
    same = np.array([epitope_of[t] for t in all_ids])
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        rows_idx = [index_of[t] for t in ids]
        block = matrix[rows_idx]
        within = 1.0 - block @ block.T
        np.fill_diagonal(within, np.inf)
        other = 1.0 - block @ matrix[same != epitope].T
        w_sorted = np.sort(within, axis=1)
        o_sorted = np.sort(other, axis=1)
        _, fam = connected_components(
            sp.csr_matrix(within <= RADIUS), directed=False)
        sizes = pd.Series(fam).value_counts()
        clustered = np.isin(fam, sizes[sizes >= 2].index)
        for i, t in enumerate(ids):
            rows.append({
                "tcr_id": t, "epitope_aa": epitope, "cohort_n": len(ids),
                "nn1": w_sorted[i, 0], "nn3": w_sorted[i, :3].mean(),
                "bg_nn1": o_sorted[i, 0], "bg_nn3": o_sorted[i, :3].mean(),
                "clustered_012": bool(clustered[i]),
                "family_size": int(sizes[fam[i]]),
            })
    frame = pd.DataFrame(rows).set_index("tcr_id")
    frame["enrichment_nn3"] = frame["bg_nn3"] / frame["nn3"]
    return frame


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    vectors, position = space()
    nb = neighbour_table(test, vectors, position)

    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    cons = pd.read_csv(ROOT / "data/interim/boltz_consistency.csv").set_index("tcr_id")
    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    feats = pd.DataFrame(index=nb.index)
    feats["log10_pgen_paired"] = pgen["log10_pgen_paired"].reindex(nb.index)
    for column in ("ms_iptm_mean", "ms_iptm_sd", "consist_tcr_rmsd_mean"):
        feats[column] = cons[column].reindex(nb.index)
    for column in ("complex_plddt", "pae_tcr_pep_min", "iptm"):
        feats[column] = boltz[column].reindex(nb.index)
    feats["typicality"] = typicality(test).reindex(nb.index)
    ordered = test.set_index("tcr_id")
    feats["cdr3_len_total"] = (ordered["cdr3_alpha_aa"].astype(str).str.len()
                               + ordered["cdr3_beta_aa"].astype(str).str.len()
                               ).reindex(nb.index)

    # --------------------------------------------------------------
    print("=== Part 1a: is singleton status predictable from the TCR itself? ===")
    print("  AUC for predicting 'is a singleton', computed WITHIN each epitope")
    print("  so cohort size is held fixed. The target carries no label, so this")
    print("  is measurable on the private epitopes too.\n")
    target = (~nb["clustered_012"]).astype(int)
    rows = []
    for column, cause in CAUSE_FEATURES.items():
        per = {}
        for epitope, block in nb.groupby("epitope_aa"):
            y = target.loc[block.index].to_numpy()
            x = feats.loc[block.index, column].to_numpy(dtype=float)
            keep = np.isfinite(x)
            if len(np.unique(y[keep])) > 1 and keep.sum() > 20:
                per[epitope] = roc_auc_score(y[keep], x[keep])
        rows.append({"feature": column, "cause": cause,
                     "mean_auc": float(np.mean(list(per.values()))),
                     "min": float(np.min(list(per.values()))),
                     "max": float(np.max(list(per.values())))})
    cause_table = pd.DataFrame(rows).sort_values("mean_auc", ascending=False)
    print(cause_table.round(3).to_string(index=False))
    print("\n  0.5 means the feature says nothing about whether a TCR has")
    print("  relatives; above 0.5 means the feature is HIGH for singletons,")
    print("  below 0.5 means it is high for clustered TCRs.")

    combined = []
    for epitope, block in nb.groupby("epitope_aa"):
        x = feats.loc[block.index].to_numpy(dtype=float)
        x = np.nan_to_num(x, nan=np.nanmedian(x))
        x = (x - x.mean(axis=0)) / (x.std(axis=0) + 1e-9)
        y = target.loc[block.index].to_numpy()
        if len(np.unique(y)) < 2:
            continue
        model = LogisticRegression(max_iter=2000, C=0.5).fit(x, y)
        combined.append(roc_auc_score(y, model.decision_function(x)))
    print(f"\n  all nine features together, in-sample within epitope: "
          f"{np.mean(combined):.3f}")
    print("  (in-sample, so an upper bound on how intrinsic singleton status is)")

    print("\n=== Part 1b: the sampling explanation, at the epitope level ===")
    rate = nb.groupby("epitope_aa").agg(
        cohort_n=("cohort_n", "first"),
        singleton_rate=("clustered_012", lambda s: 1 - s.mean()))
    rho = spearmanr(rate["cohort_n"], rate["singleton_rate"])
    print(rate.round(3).to_string())
    print(f"\n  singleton rate against cohort size: rho {rho.statistic:+.3f} "
          f"(n={len(rate)})")

    # --------------------------------------------------------------
    print("\n=== Part 2a: the distance distribution against a background ===")
    print("  'bg' is the distance to the nearest TCRs of OTHER epitopes -- how")
    print("  close a TCR would be to an unrelated one anyway.\n")
    shape = nb.groupby("epitope_aa").agg(
        nn3_q10=("nn3", lambda s: s.quantile(0.10)),
        nn3_median=("nn3", "median"),
        nn3_q90=("nn3", lambda s: s.quantile(0.90)),
        bg_nn3_median=("bg_nn3", "median"),
        enrich_median=("enrichment_nn3", "median"))
    print(shape.round(3).to_string())

    print("\n  Where does within-epitope distance stop being special? The radius")
    print("  at which the observed count first matches the background count:\n")
    cross = []
    # The background median sits near 0.5, so a first version that stopped the
    # search at 0.40 returned NaN for every epitope.
    wide = np.round(np.arange(0.02, 0.901, 0.01), 3)
    for epitope, block in nb.groupby("epitope_aa"):
        point = np.nan
        for r in wide:
            observed = float((block["nn3"] <= r).mean())
            expected = float((block["bg_nn3"] <= r).mean())
            if observed > 0.01 and expected / observed >= 0.5:
                point = r
                break
        within = float((block["nn3"] <= RADIUS).mean())
        background = float((block["bg_nn3"] <= RADIUS).mean())
        cross.append({"epitope": epitope, "half_background_radius": point,
                      "frac_within_012": within,
                      "bg_frac_within_012": background,
                      "enrichment_at_012": (within / background
                                            if background > 0 else np.inf)})
    crossing = pd.DataFrame(cross).set_index("epitope")
    print(crossing.round(4).to_string())

    print("\n=== Part 2b: validation rate by nn3 decile -- knee or smooth? ===")
    joined = labels.join(nb, on="tcr_id", rsuffix="_nb")
    joined["decile"] = joined.groupby("epitope_aa")["nn3"].transform(
        lambda s: pd.qcut(s, 10, labels=False, duplicates="drop"))
    # Base rates differ fourfold across epitopes (NLV singletons 15.5% positive
    # against GLC's 51.3%), so pooling the rate directly would mix them. Take
    # each epitope's own decile rate and average the epitopes.
    per_epitope = (joined.groupby(["epitope_aa", "decile"])["pos"].mean()
                   .unstack("epitope_aa"))
    knee = pd.DataFrame({
        "n": joined.groupby("decile")["pos"].size(),
        "nn3_mean": joined.groupby("decile")["nn3"].mean(),
        "rate_pooled": joined.groupby("decile")["pos"].mean(),
        "rate_epitope_mean": per_epitope.mean(axis=1),
        "epitopes_declining": (per_epitope.diff() < 0).sum(axis=1),
    })
    print(knee.round(3).to_string())
    diffs = knee["rate_epitope_mean"].diff().dropna()
    print(f"\n  largest single-decile drop {diffs.min():+.3f} between deciles "
          f"{int(diffs.idxmin()) - 1} and {int(diffs.idxmin())}")
    print(f"  mean absolute step {diffs.abs().mean():.3f} -- a real threshold")
    print("  would put one step far above the rest.")

    print("\n=== Part 2b2: does enrichment order the PUBLISHED bands? n=8 ===")
    print("  E49 wanted this and could only get five points, because typicality's")
    print("  AUC needs labels. Enrichment is label-free, so all eight epitopes")
    print("  have an x value, and E47's published tcrdist3 AUROC gives the y.\n")
    published = {"GILGFVFTL": 0.90, "YLQPRTFLL": 0.90, "CINGVCWTV": 0.815,
                 "NLVPMVATV": 0.815, "LLWNGPMAV": 0.815, "TTDPSFLGRY": 0.815,
                 "GLCTLVAML": 0.67, "LTDEMIAQY": 0.67}
    band_rows = []
    from itertools import permutations as _perms
    for column in ("enrich_median", "frac_within_012", "nn3_median"):
        source = shape if column in shape else crossing
        x = [float(source.loc[e, column]) for e in published]
        y = [published[e] for e in published]
        rho = spearmanr(x, y)
        # n=8 and y holds only three distinct band values, so scipy's
        # asymptotic p is not trustworthy here for the same reason it was not
        # at n=5 in E49. 8! = 40,320 permutations is cheap, so take the exact
        # one. (Measured: 0.457 asymptotic against 0.510 exact -- immaterial
        # to a null result, but there is no reason to report the wrong one.)
        draws = np.array([spearmanr(x, list(o)).statistic for o in _perms(y)])
        exact = float((np.abs(draws) >= abs(rho.statistic) - 1e-12).mean())
        band_rows.append({"label_free_index": column,
                          "rho_vs_published_auroc": rho.statistic,
                          "p_exact": exact, "p_asymptotic": rho.pvalue,
                          "n": len(x)})
        print(f"  {column:<18} rho {rho.statistic:+.3f}  exact p {exact:.3f} "
              f"(asymptotic {rho.pvalue:.3f})")
    bands = pd.DataFrame(band_rows)
    print("\n  Per epitope, so the outlier is visible:")
    for epitope in sorted(published, key=lambda e: -float(shape.loc[e, 'enrich_median'])):
        print(f"    {epitope:<11} enrichment {float(shape.loc[epitope, 'enrich_median']):5.3f}"
              f"   published {published[epitope]:.3f}")

    print("\n=== Part 2c: the decisive test -- binary at its best radius ===")
    print("  against the continuous nn3_dist it is a threshold of.\n")
    sweep = []
    for r in RADII:
        indicator = pd.Series(index=test["tcr_id"], dtype=float)
        for epitope, cohort in test.groupby("epitope_aa"):
            ids = cohort["tcr_id"].to_numpy()
            block = vectors[[position[t] for t in ids]]
            d = 1.0 - block @ block.T
            np.fill_diagonal(d, np.inf)
            _, fam = connected_components(
                sp.csr_matrix(d <= r), directed=False)
            sizes = pd.Series(fam).value_counts()
            indicator.loc[ids] = np.isin(fam, sizes[sizes >= 2].index).astype(float)
        per = [roc_auc_score(b["pos"], b["tcr_id"].map(indicator))
               for _, b in labels.groupby("epitope_aa")
               if b["tcr_id"].map(indicator).nunique() > 1]
        sweep.append({"radius": r, "binary_auc": float(np.mean(per)),
                      "clustered_frac": float(indicator.mean())})
    sweep = pd.DataFrame(sweep)
    continuous = float(np.mean([
        roc_auc_score(b["pos"], -b["tcr_id"].map(nb["nn3"]))
        for _, b in labels.groupby("epitope_aa")]))
    best = sweep.loc[sweep["binary_auc"].idxmax()]
    print(sweep.round(4).to_string(index=False))
    print(f"\n  best binary indicator: radius {best['radius']:.2f} -> "
          f"{best['binary_auc']:.4f} ({best['clustered_frac']:.1%} clustered)")
    print(f"  the continuous nn3_dist:              {continuous:.4f}")
    print(f"  the radius currently used (0.12):     "
          f"{sweep.loc[sweep['radius'] == RADIUS, 'binary_auc'].iloc[0]:.4f}")
    print(f"\n  continuous beats the best binary by "
          f"{continuous - best['binary_auc']:+.4f}")
    if continuous > best["binary_auc"]:
        print("  No radius recovers what the continuous variable carries, so")
        print("  'singleton' is a lossy summary of nn3_dist, not a category.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Why a TCR is a singleton, and whether the radius is principled\n\n"
        "## Part 1a - predicting singleton status from the TCR itself\n\n"
        + cause_table.round(4).to_markdown(index=False)
        + f"\n\nAll nine together, in-sample within epitope: "
        f"**{np.mean(combined):.3f}** (an upper bound).\n\n"
        "## Part 1b - the sampling explanation\n\n"
        + rate.round(4).to_markdown()
        + f"\n\nSingleton rate against cohort size: rho "
        f"**{rho.statistic:+.3f}** (n={len(rate)}).\n\n"
        "## Part 2a - distance distribution and the background\n\n"
        + shape.round(4).to_markdown() + "\n\n" + crossing.round(4).to_markdown()
        + "\n\n## Part 2b - validation rate by nn3 decile\n\n"
        + knee.round(4).to_markdown()
        + "\n\n## Part 2b2 - label-free enrichment against the published bands\n\n"
        + bands.round(4).to_markdown(index=False)
        + "\n\n## Part 2c - binary indicator against the continuous variable\n\n"
        + sweep.round(4).to_markdown(index=False)
        + f"\n\n- best binary: radius {best['radius']:.2f}, AUC "
        f"**{best['binary_auc']:.4f}**\n"
        f"- continuous nn3_dist: **{continuous:.4f}**\n"
        f"- difference **{continuous - best['binary_auc']:+.4f}**\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
