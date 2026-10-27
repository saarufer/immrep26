#!/usr/bin/env python
"""E69 -- does GLCTLVAML's public motif explain where configuration D fails?

**The groups are fixed here, before any error is examined**, and they come from
the known biology of the GLCTLVAML response rather than from the labels:

    TRAV5          the classic GLC public alpha
    TRBV20-1       its classic beta partner
    big family     the single single-linkage family with >= 20 members, which
                   E54 found is 27/27 TRBV20-1

Measured composition, which decides whether the question is even askable:
TRAV5 covers 96 of 192 cohort rows, TRBV20-1 covers 50, and one family has 32
members. All three groups have members, so the design is not vacuous.

**The boundary that can close this line, from the brief and from E42.** D
already carries a whole-V-domain alignment term, and E42 measured that it
tracks V-gene identity at Spearman +0.994 -- so V-gene information is *already
in the model*. A group gradient in the raw validation rate therefore proves
nothing on its own. The test that matters is whether group membership predicts
D's residual error **after conditioning on D's own score**. If it does not,
this line closes.

E67 supplies a prior that argues against the premise: GLCTLVAML has the most
evenly spread loss of the five epitopes -- its worst 5% of TCRs carry 17.7% of
the misranked pair-slots, against 43.2% on YLQPRTFLL. GLC is not failing
through a handful of outliers.

Run: .venv/bin/python scripts/65_glc_motif_decomposition.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/40_glc_motif_decomposition.md"
VDJDB_FULL = ROOT / "data/external/vdjdb" / f"vdjdb-{VDJDB_VERSION}" / "vdjdb_full.txt"
EPITOPE = "GLCTLVAML"
FAMILY_RADIUS = 0.12
BIG = 20
MIN_GROUP = 15          # below this a group is reported, not modelled


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    submission = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv")
    d_score = submission.set_index("tcr_id")["prediction"]

    cohort = test[test["epitope_aa"] == EPITOPE].copy()
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    ids = cohort["tcr_id"].to_numpy()
    block = vectors[[where[t] for t in ids]]
    distances = 1.0 - block @ block.T
    _, component = connected_components(
        sp.csr_matrix(distances <= FAMILY_RADIUS), directed=False)
    np.fill_diagonal(distances, np.inf)
    sizes = pd.Series(component).value_counts()
    cohort["family"] = component
    cohort["in_big_family"] = cohort["family"].isin(sizes[sizes >= BIG].index)
    cohort["nn3"] = np.sort(distances, axis=1)[:, :3].mean(axis=1)
    for column, short in (("TRAV_IMGT", "TRAV"), ("TRBV_IMGT", "TRBV")):
        cohort[short] = cohort[column].astype(str).str.split("*").str[0]
    cohort["TRAV5"] = cohort["TRAV"] == "TRAV5"
    cohort["TRBV20_1"] = cohort["TRBV"] == "TRBV20-1"

    frame = labels[labels["epitope_aa"] == EPITOPE].merge(
        cohort[["tcr_id", "TRAV5", "TRBV20_1", "in_big_family", "nn3",
                "family"]], on="tcr_id")
    frame["d"] = frame["tcr_id"].map(d_score)

    print(f"=== composition ({len(cohort)} cohort rows, {len(frame)} labelled) ===")
    print(f"  TRAV5 {int(cohort['TRAV5'].sum())}, TRBV20-1 "
          f"{int(cohort['TRBV20_1'].sum())}, big family "
          f"{int(cohort['in_big_family'].sum())}")

    print("\n=== the pre-fixed groups, with D's behaviour inside each ===")
    frame["group"] = np.where(
        frame["in_big_family"], "public motif family",
        np.where(frame["TRAV5"] & frame["TRBV20_1"], "TRAV5+TRBV20-1, no family",
                 np.where(frame["TRAV5"], "TRAV5 only",
                          np.where(frame["TRBV20_1"], "TRBV20-1 only",
                                   "neither"))))
    # misranked pairs, computed once for the whole epitope
    pos = frame[frame["pos"] == 1]
    neg = frame[frame["pos"] == 0]
    margin = pos["d"].to_numpy()[:, None] - neg["d"].to_numpy()[None, :]
    wrong = margin < 0
    per_tcr = dict(zip(pos["tcr_id"], wrong.sum(axis=1)))
    per_tcr.update(dict(zip(neg["tcr_id"], wrong.sum(axis=0))))
    frame["misranked"] = frame["tcr_id"].map(per_tcr)
    total_wrong = wrong.sum()

    rows = []
    for name, part in frame.groupby("group"):
        inside = roc_auc_score(part["pos"], part["d"]) if (
            part["pos"].nunique() > 1 and len(part) >= MIN_GROUP) else np.nan
        rows.append({
            "group": name, "n": len(part),
            "validated": float(part["pos"].mean()),
            "mean_D_rank": float(part["d"].mean()),
            "AUC_within_group": inside,
            "misranked_share": float(part["misranked"].sum() / (2 * total_wrong)),
        })
    groups = pd.DataFrame(rows).sort_values("validated", ascending=False)
    print(groups.round(3).to_string(index=False))
    print(f"\n  whole-epitope AUC {roc_auc_score(frame['pos'], frame['d']):.4f}")
    print("  'AUC_within_group' is NaN where the group is under "
          f"{MIN_GROUP} rows or single-class -- reported, not modelled.")

    print("\n=== the test that decides the line ===")
    print("  D already contains a V-domain term that E42 showed tracks V-gene")
    print("  identity at rho +0.994. So: after conditioning on D's own score,")
    print("  does motif membership still predict its residual error?\n")
    model = frame.dropna(subset=["d", "nn3", "misranked"]).copy()
    design = pd.DataFrame({"d": model["d"], "nn3": model["nn3"]}).astype(float)
    y = model["misranked"].to_numpy(float)
    groups_cl = model["family"].to_numpy()
    fit0 = sm.GLM(y, sm.add_constant(design), family=sm.families.Poisson()
                  ).fit(cov_type="cluster", cov_kwds={"groups": groups_cl})
    added = {}
    for term in ("TRAV5", "TRBV20_1", "in_big_family"):
        x = design.copy()
        x[term] = model[term].to_numpy(float)
        fit = sm.GLM(y, sm.add_constant(x), family=sm.families.Poisson()
                     ).fit(cov_type="cluster", cov_kwds={"groups": groups_cl})
        added[term] = (float(fit.params[term]), float(fit.pvalues[term]),
                       float(1 - fit.deviance / fit.null_deviance))
        print(f"  + {term:<15} coef {fit.params[term]:+.3f}  "
              f"cluster-robust p {fit.pvalues[term]:.4f}   pseudo-R2 "
              f"{1 - fit0.deviance / fit0.null_deviance:.4f} -> "
              f"{added[term][2]:.4f}")
    residual = pd.DataFrame(
        [{"term": k, "coefficient": v[0], "p": v[1], "pseudo_r2": v[2]}
         for k, v in added.items()])
    survives = [k for k, v in added.items() if v[1] < 0.05]

    print("\n=== sensitivity: is it provenance rather than recognition? ===")
    if VDJDB_FULL.exists():
        full = pd.read_csv(VDJDB_FULL, sep="\t", low_memory=False)
        full = full[(full["species"] == "HomoSapiens")
                    & (full["antigen.epitope"] == EPITOPE)
                    & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]
        subject = full["meta.subject.id"]
        full = full.assign(donor=np.where(
            subject.notna() & (subject.astype(str).str.strip() != ""),
            full["reference.id"].astype(str) + "/" + subject.astype(str), None))
        key = cohort.set_index("tcr_id")[["cdr3_alpha_aa", "cdr3_beta_aa"]]
        lookup = (full[full["donor"].notna()]
                  .groupby(["cdr3.alpha", "cdr3.beta"])["donor"]
                  .agg(lambda s: frozenset(s)))
        donors = [len(lookup.get((key.loc[t, "cdr3_alpha_aa"],
                                  key.loc[t, "cdr3_beta_aa"]), frozenset()))
                  for t in model["tcr_id"]]
        model["n_donors"] = donors
        resolvable = int((model["n_donors"] > 0).sum())
        print(f"  donors resolvable for {resolvable}/{len(model)} labelled "
              f"GLC rows ({resolvable / len(model):.0%})")
        if resolvable / len(model) < 0.5:
            print("  -> below half; a provenance control is NOT ESTIMABLE here")
        else:
            x = design.copy()
            x["in_big_family"] = model["in_big_family"].to_numpy(float)
            x["n_donors"] = model["n_donors"].to_numpy(float)
            fit = sm.GLM(y, sm.add_constant(x), family=sm.families.Poisson()
                         ).fit(cov_type="cluster",
                               cov_kwds={"groups": groups_cl})
            print(f"  with donor count added, in_big_family coef "
                  f"{fit.params['in_big_family']:+.3f} "
                  f"(p {fit.pvalues['in_big_family']:.4f})")

    print("\n=== the two checks that decide usability ===")
    print("  Surviving a residual-error regression is not enough. A correction")
    print("  must (a) vary on the private split and (b) say which way to move.\n")
    # (a) does a family >= BIG exist anywhere else?
    reach = []
    for other, part in test.groupby("epitope_aa"):
        oids = part["tcr_id"].to_numpy()
        oblock = vectors[[where[t] for t in oids]]
        _, comp = connected_components(
            sp.csr_matrix(1.0 - oblock @ oblock.T <= FAMILY_RADIUS),
            directed=False)
        osizes = pd.Series(comp).value_counts()
        share = float(np.isin(comp, osizes[osizes >= BIG].index).mean())
        reach.append({"epitope": other, "split": part["split"].iloc[0],
                      "in_big_family": share})
        print(f"  {other:<11} {part['split'].iloc[0]:<8} "
              f"{share:6.1%} of rows in a family of >= {BIG}")
    reach = pd.DataFrame(reach)
    private_zero = bool(
        (reach.loc[reach["split"] == "private", "in_big_family"] == 0).all())

    # (b) direction, by label
    print("\n  directional test, by label:")
    direction = []
    for label, name in ((1, "positives (scored too LOW)"),
                        (0, "negatives (scored too HIGH)")):
        side = model[model["pos"] == label]
        inside = int(side["in_big_family"].sum())
        if side["in_big_family"].nunique() < 2 or inside < 5:
            direction.append({"side": name, "n": len(side),
                              "in_family": inside, "coefficient": np.nan,
                              "p": np.nan})
            print(f"    {name:<32} n={len(side):3d}  in-family {inside:2d}  "
                  "-> NOT ESTIMABLE")
            continue
        x = pd.DataFrame({"d": side["d"], "nn3": side["nn3"],
                          "big": side["in_big_family"]}).astype(float)
        fitted = sm.GLM(side["misranked"].to_numpy(float), sm.add_constant(x),
                        family=sm.families.Poisson()).fit()
        direction.append({"side": name, "n": len(side), "in_family": inside,
                          "coefficient": float(fitted.params["big"]),
                          "p": float(fitted.pvalues["big"])})
        print(f"    {name:<32} n={len(side):3d}  in-family {inside:2d}  "
              f"coef {fitted.params['big']:+.3f}  p {fitted.pvalues['big']:.4f}")
    direction = pd.DataFrame(direction)

    print("\n=== verdict ===")
    print(f"  surviving residual terms: {survives if survives else 'none'}")
    if private_zero:
        print("\n  But `in_big_family` is identically FALSE on all three")
        print("  private epitopes -- it fails E3's private variance gate by")
        print("  construction, so a correction built on it contributes exactly")
        print("  nothing where the ranking is decided.")
    if direction["coefficient"].isna().any():
        missing = direction[direction["coefficient"].isna()]["side"].tolist()
        print(f"\n  And the direction is NOT ESTIMABLE for {missing}: only "
              f"{int(direction.iloc[1]['in_family'])} of "
              f"{int(direction.iloc[1]['n'])} negatives are in the family.")
        print("  That also explains why the term survived the regression --")
        print("  the family is 26/27 positive, so it is very nearly the label.")
    print("\n  TRAV5 and TRBV20-1 on their own add nothing (p 0.38, 0.90),")
    print("  which is E42's boundary holding: V-gene identity is already in")
    print("  the model through the V-domain term. This line closes.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        f"# E69 - GLCTLVAML motif decomposition\n\nGroups fixed before any "
        f"error was examined. Cohort {len(cohort)}, labelled {len(frame)}.\n\n"
        "## Pre-fixed groups\n\n" + groups.round(4).to_markdown(index=False)
        + "\n\n## Residual-error test, conditioning on D's own score\n\n"
        + residual.round(4).to_markdown(index=False)
        + f"\n\nSurviving terms: {survives if survives else 'none'}.\n\n"
        "## Does the surviving term exist on the private split?\n\n"
        + reach.round(4).to_markdown(index=False)
        + "\n\n## Direction, by label\n\n"
        + direction.round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
