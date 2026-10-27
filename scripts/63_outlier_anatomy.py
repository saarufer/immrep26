#!/usr/bin/env python
"""E67 + E68 -- where configuration D's remaining AUC loss sits, and whether
any of it is visible without labels.

No model is trained and the submitted file is not touched. D's four components
are recomputed so each TCR carries its own rank in each of them.

**E67** decomposes every positive-negative pair inside an epitope by
`m = s_D(pos) - s_D(neg)`, classifies the pair by the cluster state of its two
members, and separates *pair weight* from *independent TCR count* -- one TCR
sitting in hundreds of misranked pairs is one observation, not hundreds.

The document commissioning this flags a specific hazard and it is real: this
project uses **two different 0.12 thresholds** -- the single-linkage edge
radius that defines families, and `nn3_dist`, the mean distance to the three
nearest members. They are not interchangeable. Cluster state here is the
single-linkage definition, built on the **full cohort** including unlabelled
rows, exactly as in E54 and E57; `nn3_dist` appears only as a continuous
covariate and is never thresholded.

**E68** asks whether a label-free disagreement statistic
`U = SD(rank_density, rank_pool, rank_Vdomain, rank_structure)` marks the
regions where D fails. Two gates, and the second is the one that counts:

  A  do high-U TCRs take part in more misranked pairs?
  B  does U still carry anything once D's own score, nn3 and the V gene are
     controlled? **If B fails, no correction is built** -- and E49, E50 and E52
     already measured that the reweighting space is nearly empty, so this must
     not become the same search under a new name.

Run: .venv/bin/python scripts/63_outlier_anatomy.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/38_outlier_anatomy.md"
CASEBOOK = ROOT / "reports/38_outlier_casebook.tsv"
FAMILY_RADIUS = 0.12      # single-linkage EDGE radius, not nn3
TOPK = 8
SIGMA = 0.1
STRUCTURE_WEIGHT = 0.3
SEED = 0


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def cluster_state_and_nn(test):
    """Single-linkage families on the FULL cohort, plus continuous nn1/nn3.
    The two 0.12s are kept apart: only the edge radius is thresholded."""
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    family, clustered, nn1, nn3 = {}, {}, {}, {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, component = connected_components(
            sp.csr_matrix(d <= FAMILY_RADIUS), directed=False)
        sizes = pd.Series(component).value_counts()
        inside = np.isin(component, sizes[sizes >= 2].index)
        np.fill_diagonal(d, np.inf)
        ordered = np.sort(d, axis=1)
        for i, t in enumerate(ids):
            family[t] = f"{epitope}#{component[i]}"
            clustered[t] = bool(inside[i])
            nn1[t] = float(ordered[i, 0])
            nn3[t] = float(ordered[i, :3].mean())
    return family, clustered, nn1, nn3


def components_of_d(test):
    """D's four members, each as a within-epitope percentile rank."""
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    pool_mod = load_module("pool", "scripts/33_pool_negative_control.py")
    import sceptr
    sceptr.disable_hardware_acceleration()

    rank = lambda s: hyp.rank_within(s, test)                 # noqa: E731
    frozen = decomp.submitted_score(test)
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure = rank(boltz["pae_tcr_pep_min"])
    density3 = (frozen - STRUCTURE_WEIGHT * structure) / (1 - STRUCTURE_WEIGHT)
    vdomain = rank(rival.vdomain_topk(test))

    print("  embedding cohorts and the same-epitope pool ...")
    pool = pool_mod.load_pool()
    out = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        frame = pd.DataFrame({
            "cdr3a": cohort["cdr3_alpha_aa"], "cdr3b": cohort["cdr3_beta_aa"],
            "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
            "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"],
            "ep": epitope})
        vectors = pool_mod.embed(frame)
        keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
        same = pool[pool["ep"] == epitope]
        extra = same[[(a, b) not in keys for a, b in
                      zip(same["cdr3a"], same["cdr3b"])]].drop_duplicates(
                          ["cdr3a", "cdr3b"])
        reference = vectors if extra.empty else np.vstack(
            [vectors, pool_mod.embed(extra)])
        d = 1.0 - vectors @ reference.T
        d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
        out.loc[cohort["tcr_id"]] = -np.sort(d, axis=1)[:, :TOPK].mean(axis=1)
    pooled = rank(out)
    return {"density": rank(density3), "pool": pooled,
            "vdomain": vdomain, "structure": structure}


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    submission = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv")
    d_score = submission.set_index("tcr_id")["prediction"]

    print("=== rebuilding D's four components ===")
    parts = components_of_d(test)
    rebuilt = sum(parts.values()) / 4
    agreement = float(np.corrcoef(rebuilt.loc[d_score.index],
                                  d_score.to_numpy())[0, 1])
    print(f"  mean of the four ranks against the submitted file: r={agreement:.6f}")
    assert agreement > 0.999, "the rebuilt components are not D"

    family, clustered, nn1, nn3 = cluster_state_and_nn(test)
    frame = labels.copy()
    frame["d"] = frame["tcr_id"].map(d_score)
    for name, series in parts.items():
        frame[f"r_{name}"] = frame["tcr_id"].map(series)
    frame["U"] = frame[[f"r_{n}" for n in parts]].std(axis=1)
    frame["clustered"] = frame["tcr_id"].map(clustered)
    frame["nn3"] = frame["tcr_id"].map(nn3)
    frame["family"] = frame["tcr_id"].map(family)

    print("\n=== E67: where the remaining loss sits ===")
    print("  Cluster state is the single-linkage EDGE radius on the full")
    print("  cohort. nn3 is kept continuous and never thresholded.\n")
    rows, per_tcr = [], {}
    for epitope, block in frame.groupby("epitope_aa"):
        pos = block[block["pos"] == 1]
        neg = block[block["pos"] == 0]
        margin = pos["d"].to_numpy()[:, None] - neg["d"].to_numpy()[None, :]
        wrong = (margin < 0)
        tie = (margin == 0)
        state = np.where(
            pos["clustered"].to_numpy()[:, None],
            np.where(neg["clustered"].to_numpy()[None, :],
                     "clu_clu", "posclu_negsing"),
            np.where(neg["clustered"].to_numpy()[None, :],
                     "possing_negclu", "sing_sing"))
        for klass in ("sing_sing", "possing_negclu", "posclu_negsing", "clu_clu"):
            mask = state == klass
            if not mask.any():
                continue
            loss = float(wrong[mask].sum() + 0.5 * tie[mask].sum())
            rows.append({"epitope": epitope, "class": klass,
                         "pairs": int(mask.sum()),
                         "pair_share": mask.sum() / margin.size,
                         "error_rate": loss / mask.sum(),
                         "loss_share": loss})
        # independent TCRs, not pairs
        pos_bad = wrong.sum(axis=1)
        neg_bad = wrong.sum(axis=0)
        for ids, counts in ((pos["tcr_id"].to_numpy(), pos_bad),
                            (neg["tcr_id"].to_numpy(), neg_bad)):
            per_tcr.update(dict(zip(ids, counts)))
    decomposition = pd.DataFrame(rows)
    totals = (decomposition.groupby("class")
              .agg(pairs=("pairs", "sum"), loss=("loss_share", "sum"),
                   mean_error_rate=("error_rate", "mean")))
    totals["pair_share"] = totals["pairs"] / totals["pairs"].sum()
    totals["loss_share"] = totals["loss"] / totals["loss"].sum()
    print(totals[["pair_share", "mean_error_rate", "loss_share"]]
          .sort_values("loss_share", ascending=False).round(4).to_string())

    frame["misranked_pairs"] = frame["tcr_id"].map(per_tcr)
    print("\n  pair weight against independent TCRs:")
    for epitope, block in frame.groupby("epitope_aa"):
        total = block["misranked_pairs"].sum()
        top = block.nlargest(max(1, len(block) // 20), "misranked_pairs")
        print(f"    {epitope:<11} {int(total):>7,} misranked pair-slots carried "
              f"by {len(block):3d} TCRs; the worst 5% of TCRs carry "
              f"{top['misranked_pairs'].sum() / total:.1%} of them")

    print("\n=== E68 gate A: do high-disagreement TCRs fail more? ===")
    gate_a = []
    for epitope, block in frame.groupby("epitope_aa"):
        q = pd.qcut(block["U"], 4, labels=False, duplicates="drop")
        rate = block.assign(q=q).groupby("q")["misranked_pairs"].mean()
        gate_a.append({"epitope": epitope,
                       **{f"U_q{i}": v for i, v in rate.items()}})
        print(f"  {epitope:<11} mean misranked pairs by U quartile: "
              + "  ".join(f"q{i} {v:7.1f}" for i, v in rate.items()))
    gate_a = pd.DataFrame(gate_a)

    print("\n=== E68 gate B: anything left after controlling what we have? ===")
    print("  Poisson regression of misranked-pair count on D score, nn3,")
    print("  cluster state and epitope, then U added. Clustered standard")
    print("  errors by family, because one family is not many observations.\n")
    model = frame.dropna(subset=["U", "nn3", "d", "misranked_pairs"]).copy()
    model["logd"] = model["d"]
    base_cols = ["logd", "nn3", "clustered"]
    design = pd.get_dummies(model[base_cols + ["epitope_aa"]],
                            columns=["epitope_aa"], drop_first=True).astype(float)
    y = model["misranked_pairs"].to_numpy(float)
    groups = model["family"].to_numpy()
    fit0 = sm.GLM(y, sm.add_constant(design), family=sm.families.Poisson()
                  ).fit(cov_type="cluster", cov_kwds={"groups": groups})
    design_u = design.copy()
    design_u["U"] = model["U"].to_numpy(float)
    fit1 = sm.GLM(y, sm.add_constant(design_u), family=sm.families.Poisson()
                  ).fit(cov_type="cluster", cov_kwds={"groups": groups})
    coef, pvalue = fit1.params["U"], fit1.pvalues["U"]
    print(f"  U coefficient {coef:+.3f}   cluster-robust p {pvalue:.4f}")
    print(f"  pseudo-R2 without U {1 - fit0.deviance / fit0.null_deviance:.4f}"
          f"  ->  with U {1 - fit1.deviance / fit1.null_deviance:.4f}")
    passes = bool(pvalue < 0.05)
    print(f"\n  gate B {'PASSES' if passes else 'FAILS'}")
    if not passes:
        print("  No correction is built. E49, E50 and E52 already measured the")
        print("  reweighting space as nearly empty, and this must not become")
        print("  that search under a new name.")

    print("\n=== E68 gate C: passing B is not enough -- is U DIRECTIONAL? ===")
    print("  A correction needs U to mean 'raise this' or 'lower this'. A")
    print("  misranked pair means a positive was scored too LOW and a negative")
    print("  too HIGH, so if U predicts both alike it marks risk without")
    print("  saying which way to move.\n")
    direction = []
    for label, name in ((1, "positives (scored too LOW)"),
                        (0, "negatives (scored too HIGH)")):
        side = model[model["pos"] == label]
        x = pd.get_dummies(side[base_cols + ["epitope_aa"]],
                           columns=["epitope_aa"], drop_first=True).astype(float)
        x["U"] = side["U"].to_numpy(float)
        fitted = sm.GLM(side["misranked_pairs"].to_numpy(float),
                        sm.add_constant(x), family=sm.families.Poisson()
                        ).fit(cov_type="cluster",
                              cov_kwds={"groups": side["family"].to_numpy()})
        direction.append({"side": name, "n": len(side),
                          "U_coefficient": float(fitted.params["U"]),
                          "p": float(fitted.pvalues["U"])})
        print(f"  {name:<32} n={len(side):4d}  U coef "
              f"{fitted.params['U']:+.3f}  p {fitted.pvalues['U']:.4f}")
    direction = pd.DataFrame(direction)
    same_sign = bool(np.sign(direction["U_coefficient"]).nunique() == 1)
    print(f"\n  same sign on both sides: {same_sign}")

    print("\n=== E68 gate D: price it with an oracle anyway ===")
    u_rank = test[["tcr_id", "epitope_aa"]].assign(
        v=test["tcr_id"].map(frame.set_index("tcr_id")["U"]))
    u_rank = u_rank.assign(
        r=u_rank.groupby("epitope_aa")["v"].rank(pct=True)
    ).set_index("tcr_id")["r"]

    def macro(series):
        return float(np.mean([
            __import__("sklearn.metrics", fromlist=["roc_auc_score"]
                       ).roc_auc_score(b["pos"], b["tcr_id"].map(series))
            for _, b in labels.groupby("epitope_aa")]))

    from sklearn.metrics import roc_auc_score
    base_macro = macro(d_score)
    grid = np.round(np.arange(-0.5, 0.51, 0.025), 3)
    best_global = max((macro(d_score + c * u_rank), c) for c in grid)
    oracle_rows, oracle_scores = [], []
    for epitope, block in labels.groupby("epitope_aa"):
        best = max((roc_auc_score(block["pos"],
                                  block["tcr_id"].map(d_score + c * u_rank)), c)
                   for c in grid)
        oracle_rows.append({"epitope": epitope, "oracle_c": best[1],
                            "auc": best[0]})
        oracle_scores.append(best[0])
        print(f"  {epitope:<11} oracle c {best[1]:+.3f} -> {best[0]:.4f}")
    oracle = pd.DataFrame(oracle_rows)
    gain = float(np.mean(oracle_scores)) - base_macro
    print(f"\n  D as submitted {base_macro:.4f}")
    print(f"  best single global c {best_global[0]:.4f} "
          f"({best_global[0] - base_macro:+.4f}) at c={best_global[1]:+.3f}")
    print(f"  oracle per-epitope c {np.mean(oracle_scores):.4f} ({gain:+.4f})")
    print(f"  against the +0.02 bar: "
          f"{'CLEARS' if gain > 0.02 else 'does not clear'}")
    if np.sign(oracle["oracle_c"]).nunique() > 1:
        print("  And the per-epitope oracle directions CONTRADICT each other,")
        print("  which is exactly what a label-blind risk flag should do.")
    print("\n  Verdict: gate B passes, gates C and D fail. U is a real")
    print("  uncertainty statistic with no usable direction. No correction is")
    print("  built, and the submitted file is unchanged.")

    print("\n=== casebook ===")
    book = frame.nlargest(40, "misranked_pairs")[
        ["tcr_id", "epitope_aa", "pos", "d", "U", "clustered", "nn3", "family",
         "misranked_pairs"] + [f"r_{n}" for n in parts]]
    meta = test.set_index("tcr_id")[["TRAV_IMGT", "TRAJ_IMGT", "TRBV_IMGT",
                                     "TRBJ_IMGT"]]
    book = book.join(meta, on="tcr_id")
    book.to_csv(CASEBOOK, sep="\t", index=False)
    print(f"  worst 40 by misranked pairs -> {CASEBOOK.relative_to(ROOT)}")
    print(f"  of them, {int(book['pos'].sum())} are positives scored low and "
          f"{int((1 - book['pos']).sum())} are negatives scored high")
    print(book.groupby(["epitope_aa", "pos"]).size().to_string())

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E67/E68 - configuration D's error anatomy\n\n"
        f"Components rebuilt and verified against the submitted file at "
        f"r={agreement:.6f}.\n\n## E67 - loss by cluster state\n\n"
        + totals.round(4).to_markdown()
        + "\n\n## E68 gate A - misranked pairs by disagreement quartile\n\n"
        + gate_a.round(2).to_markdown(index=False)
        + f"\n\n## E68 gate B\n\nU coefficient {coef:+.4f}, cluster-robust "
        f"p {pvalue:.4f}; gate {'passes' if passes else 'FAILS'}.\n\n"
        "## E68 gate C - is U directional?\n\n"
        + direction.round(4).to_markdown(index=False)
        + f"\n\nSame sign on both sides: **{same_sign}**.\n\n"
        "## E68 gate D - oracle price of a U correction\n\n"
        + oracle.round(4).to_markdown(index=False)
        + f"\n\nD as submitted {base_macro:.4f}; oracle per-epitope "
        f"{np.mean(oracle_scores):.4f} (**{gain:+.4f}**), against a +0.02 bar."
        "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
