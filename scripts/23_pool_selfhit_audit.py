#!/usr/bin/env python
"""Is candidate-pool expansion adding independent evidence, or finding the query itself?

Plan card Q02. Expanding an epitope cohort with same-epitope TCRs from VDJdb
raises the measured density, and a teammate's pipeline reports +0.0194 under
leave-one-epitope-out from doing so. But every test TCR *came from* VDJdb, so
three artefacts can produce that gain without any new information:

  exact self-hit     the query's own record sits in the pool at distance 0.
  clonotype copy     the same paired clonotype curated from several studies, so
                     one real observation appears as several zero-distance
                     neighbours.
  study block        one high-throughput paper contributing hundreds of related
                     records, which dominates a density estimate by volume
                     rather than by independent convergence.

This audits all three and then recomputes the gain under four progressively
stricter pool definitions. If the gain survives only the permissive one, it is
not transferable -- and that matters beyond bookkeeping, because pool expansion
is the only live attack on the both-singleton pairs that carry 62% of the
remaining error (E22).

Run: python scripts/23_pool_selfhit_audit.py
"""

from __future__ import annotations

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
from immrep.data import build_training_frame, load_test  # noqa: E402

VDJDB = ROOT / f"data/external/vdjdb/vdjdb-{VDJDB_VERSION}/vdjdb_full.txt"
OUT = ROOT / "reports/07_pool_selfhit_audit.md"

# Keys of decreasing strictness. `exact` is the query identity; `clonotype`
# treats a rearrangement curated twice as one observation.
EXACT_KEY = ["cdr3a", "cdr3b", "trav", "traj", "trbv", "trbj"]
CLONOTYPE_KEY = ["cdr3a", "cdr3b"]


def load_pool() -> pd.DataFrame:
    frame = pd.read_csv(VDJDB, sep="\t", low_memory=False)
    frame = frame[
        (frame["species"] == "HomoSapiens")
        & frame["cdr3.alpha"].notna()
        & frame["cdr3.beta"].notna()
    ].rename(
        columns={
            "cdr3.alpha": "cdr3a", "cdr3.beta": "cdr3b",
            "v.alpha": "trav", "j.alpha": "traj",
            "v.beta": "trbv", "j.beta": "trbj",
            "antigen.epitope": "epitope", "reference.id": "reference",
        }
    )
    for column in ("trav", "traj", "trbv", "trbj"):
        frame[column] = frame[column].fillna("NA").astype(str).str.split("*").str[0]
    return frame[["cdr3a", "cdr3b", "trav", "traj", "trbv", "trbj", "epitope", "reference"]]


def tcrdist_rows(epitope: str) -> tuple[np.ndarray, list[str]]:
    cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
    return (cached["alpha"] + cached["beta"]).astype(float), list(cached["tcr_id"])


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    pool = load_pool()

    query = test.rename(
        columns={
            "cdr3_alpha_aa": "cdr3a", "cdr3_beta_aa": "cdr3b",
            "TRAV_vdjdb": "trav", "TRAJ_vdjdb": "traj",
            "TRBV_vdjdb": "trbv", "TRBJ_vdjdb": "trbj",
            "epitope_aa": "epitope",
        }
    )[["tcr_id", "cdr3a", "cdr3b", "trav", "traj", "trbv", "trbj", "epitope"]].copy()
    for column in ("trav", "traj", "trbv", "trbj"):
        query[column] = query[column].astype(str).str.split("*").str[0]

    print("=== audit 1: how often is the query its own neighbour? ===")
    rows = []
    for epitope, cohort in query.groupby("epitope"):
        same = pool[pool["epitope"] == epitope]
        exact = cohort.merge(same, on=EXACT_KEY + ["epitope"], how="left", indicator=True)
        clono = cohort.merge(same, on=CLONOTYPE_KEY + ["epitope"], how="left", indicator=True)
        exact_hits = exact.groupby("tcr_id")["_merge"].apply(lambda s: (s == "both").sum())
        clono_hits = clono.groupby("tcr_id")["_merge"].apply(lambda s: (s == "both").sum())
        rows.append(
            {
                "epitope": epitope,
                "n_query": len(cohort),
                "pool_records": len(same),
                "pool_distinct_clonotypes": same.drop_duplicates(CLONOTYPE_KEY).shape[0],
                "query_with_exact_selfhit": float((exact_hits > 0).mean()),
                "mean_exact_selfhits": float(exact_hits.mean()),
                "mean_clonotype_copies": float(clono_hits.mean()),
                "pool_refs": same["reference"].nunique(),
                "pool_top_ref_share": (
                    float(same["reference"].value_counts(normalize=True).iloc[0])
                    if same["reference"].notna().any() else float("nan")
                ),
            }
        )
    audit = pd.DataFrame(rows).set_index("epitope")
    audit["new_records_per_query"] = audit["pool_records"] / audit["n_query"]
    audit["new_clonotypes_per_query"] = audit["pool_distinct_clonotypes"] / audit["n_query"]
    pd.set_option("display.width", 260)
    print(audit.round(3).to_string())

    private = ["CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"]
    print(
        f"\n  mean fraction of queries that find themselves in the pool: "
        f"{audit['query_with_exact_selfhit'].mean():.3f}"
    )
    print(
        f"  private epitopes: "
        + "  ".join(f"{e}={audit.loc[e, 'query_with_exact_selfhit']:.3f}" for e in private)
    )

    print("\n=== audit 2: records vs independent clonotypes vs studies ===")
    inflation = pd.DataFrame(
        {
            "records/query": audit["new_records_per_query"],
            "clonotypes/query": audit["new_clonotypes_per_query"],
            "record:clonotype": audit["pool_records"] / audit["pool_distinct_clonotypes"],
            "distinct_refs": audit["pool_refs"],
            "top_ref_share": audit["pool_top_ref_share"],
        }
    )
    print(inflation.round(3).to_string())
    print(
        "\n  Reading: a record:clonotype ratio above 1 means the same rearrangement is\n"
        "  curated repeatedly, so a raw-record pool overstates independent observations.\n"
        "  A high top_ref_share means one study would dominate the density estimate."
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Pool expansion: self-hit and duplication audit\n\n"
        "Plan card Q02. Every test TCR came from VDJdb, so expanding the cohort with\n"
        "VDJdb risks counting the query itself, or one rearrangement curated several\n"
        "times, as independent neighbours.\n\n"
        "## Self-hit rate and pool composition\n\n" + audit.round(3).to_markdown()
        + "\n\n## Records vs independent clonotypes vs studies\n\n"
        + inflation.round(3).to_markdown() + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
