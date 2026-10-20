#!/usr/bin/env python
"""Decide the final submission: density alone, or density plus structure/specificity.

The decision rule is fixed **before** looking at the numbers, because choosing
components by their leave-one-epitope-out score is itself fitting on five folds
-- the thing that already cost 0.006 in scripts/07_loeo_fusion.py.

  Include a component only if rank-averaging it with density improves the
  LOEO mean by more than 0.02, i.e. about two standard errors of the
  five-epitope estimate (SE ~0.011). Otherwise keep density alone.

Fusion is a plain rank average with no fitted weights, for the same reason.

Run: python scripts/12_final_fusion.py
"""

from __future__ import annotations

import itertools
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.metric import per_epitope_auc, score_report  # noqa: E402

MIN_GAIN_TO_INCLUDE = 0.02

DENSITY_SPACES = ("default", "small")
SIGMA = 0.1
STRUCTURE_FILE = ROOT / "data/interim/boltz_features.csv"
SPECIFICITY_FILE = ROOT / "data/interim/specificity_scores.csv"

# Structure features worth testing: the organisers' statistic, Boltz's own
# interface confidence, and the per-chain splits since the alpha chain carries
# more of the sequence signal than beta.
STRUCTURE_COLUMNS = (
    "pae_tcr_pep_min",
    "pae_a_pep_p10",
    "pae_b_pep_p10",
    "iptm_tcr_pep_max",
    "iptm_tcr_pep_mean",
    "iptm_a_pep",
    "iptm_b_pep",
    "iptm",
    "confidence_score",
    "plddt_peptide",
)


def density_scores(test: pd.DataFrame) -> pd.Series:
    """Rank-averaged kernel density over the chosen SCEPTR spaces."""
    ranks = []
    for name in DENSITY_SPACES:
        cached = np.load(ROOT / f"data/interim/sceptr/{name}.npz", allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        vectors = cached["vectors"]
        values = pd.Series(index=test["tcr_id"], dtype=float)
        for _, cohort in test.groupby("epitope_aa"):
            rows = [position[t] for t in cohort["tcr_id"]]
            block = vectors[rows]
            block = block / np.linalg.norm(block, axis=1, keepdims=True)
            distances = 1.0 - block @ block.T
            np.fill_diagonal(distances, np.inf)
            values.loc[cohort["tcr_id"]] = np.exp(-distances / SIGMA).mean(axis=1)
        ranks.append(rank_within_epitope(values, test))
    return sum(ranks) / len(ranks)


def rank_within_epitope(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    """Percentile rank inside each epitope -- the only scale the metric reads."""
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["value"] = frame["tcr_id"].map(values)
    frame["rank"] = frame.groupby("epitope_aa")["value"].rank(pct=True)
    return frame.set_index("tcr_id")["rank"]


def main() -> int:
    test = load_test()
    labelled = build_training_frame()

    components: dict[str, pd.Series] = {"density": density_scores(test)}

    if STRUCTURE_FILE.exists():
        structure = pd.read_csv(STRUCTURE_FILE).set_index("tcr_id")
        coverage = structure.index.intersection(test["tcr_id"])
        print(f"structure features: {len(coverage)}/{len(test)} complexes folded")
        if len(coverage) == len(test):
            for column in STRUCTURE_COLUMNS:
                if column in structure.columns:
                    components[f"struct:{column}"] = rank_within_epitope(
                        structure[column], test
                    )
        else:
            print("  -> incomplete, structure held out of the decision")
    else:
        print("structure features: not yet available")

    if SPECIFICITY_FILE.exists():
        spec = pd.read_csv(SPECIFICITY_FILE, index_col=0).iloc[:, 0]
        components["specificity"] = rank_within_epitope(spec, test)
        print("specificity score: loaded")

    # Standalone scores first.
    print("\n=== standalone (per-epitope AUC on the five public epitopes) ===")
    standalone = {}
    for name, series in components.items():
        aucs = per_epitope_auc(labelled, labelled["tcr_id"].map(series).to_numpy())
        standalone[name] = aucs.mean()
        print(f"  {name:<28} {aucs.mean():.4f}   " +
              "  ".join(f"{k[:5]}={v:.3f}" for k, v in aucs.items()))

    baseline = standalone["density"]
    print(f"\ndensity alone = {baseline:.4f}; a component must add "
          f">{MIN_GAIN_TO_INCLUDE} to be included")

    # Pairwise and triple rank averages against density.
    others = [n for n in components if n != "density"]
    print("\n=== rank-average with density ===")
    results = []
    for size in (1, 2):
        for combo in itertools.combinations(others, size):
            blended = sum(
                [components["density"]] + [components[c] for c in combo]
            ) / (1 + len(combo))
            aucs = per_epitope_auc(labelled, labelled["tcr_id"].map(blended).to_numpy())
            results.append(
                {"with": " + ".join(combo), "mean": aucs.mean(),
                 "gain": aucs.mean() - baseline, **aucs.round(3).to_dict()}
            )
    if results:
        table = pd.DataFrame(results).sort_values("mean", ascending=False)
        pd.set_option("display.width", 240)
        print(table.head(15).to_string(index=False))

        best = table.iloc[0]
        print(f"\nbest blend: density + {best['with']}  "
              f"mean={best['mean']:.4f}  gain={best['gain']:+.4f}")
        verdict = (
            f"INCLUDE {best['with']}" if best["gain"] > MIN_GAIN_TO_INCLUDE
            else "KEEP DENSITY ALONE -- no component clears the pre-set threshold"
        )
        print(f"decision: {verdict}")
    else:
        print("  no additional components available yet")

    print("\n=== density alone, full report ===")
    print(score_report(labelled, labelled["tcr_id"].map(components["density"]).to_numpy())
          .round(4).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
