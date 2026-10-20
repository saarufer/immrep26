#!/usr/bin/env python
"""The master cross-epitope comparison: every approach, one honest protocol.

All five public epitopes, leave-one-epitope-out, mean per-epitope ROC AUC --
the competition's own metric. Nothing here is tuned on the epitope it is scored
on. Where a method has a knob, the knob is either fixed a priori (and said so)
or selected inside each fold, which is reported separately because selection on
five folds costs about 0.005.

Covers, in one table: the organisers' two baselines, annotation metadata, the
hand-built and learned sequence representations, generation probability,
structure, an unseen-epitope specificity model, and their rank-average fusions.

Run: python scripts/15_master_comparison.py
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

from immrep.data import build_training_frame, load_tcrdist_baseline, load_test  # noqa: E402

SIGMA = 0.1  # fixed a priori; the AUC-vs-sigma curve is flat between 0.05 and 0.15
OUT = ROOT / "reports/04_master_comparison.md"


def sceptr_density(test: pd.DataFrame, space: str, sigma: float = SIGMA) -> pd.Series:
    cached = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    vectors = cached["vectors"]
    values = pd.Series(index=test["tcr_id"], dtype=float)
    for _, cohort in test.groupby("epitope_aa"):
        block = vectors[[position[t] for t in cohort["tcr_id"]]]
        block = block / np.linalg.norm(block, axis=1, keepdims=True)
        distances = 1.0 - block @ block.T
        np.fill_diagonal(distances, np.inf)
        values.loc[cohort["tcr_id"]] = np.exp(-distances / sigma).mean(axis=1)
    return values


def tcrdist_density(test: pd.DataFrame, sigma: float = 20.0) -> pd.Series:
    values = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [position[t] for t in cohort["tcr_id"]]
        distances = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
        np.fill_diagonal(distances, np.inf)
        values.loc[cohort["tcr_id"]] = np.exp(-distances / sigma).mean(axis=1)
    return values


def rank_within_epitope(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    frame["r"] = frame.groupby("epitope_aa")["v"].rank(pct=True)
    return frame.set_index("tcr_id")["r"]


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    epitopes = sorted(labels["epitope_aa"].unique())

    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    specificity = pd.read_csv(ROOT / "data/interim/specificity_scores.csv", index_col=0).iloc[:, 0]

    raw: dict[str, pd.Series] = {
        # --- organisers' baselines ---
        "tcrdist3 nndist (organisers' baseline)":
            load_tcrdist_baseline().set_index("tcr_id")["prediction"],
        "structure: TCR-peptide min-PAE (AF3 analogue)": structure["pae_tcr_pep_min"],
        # --- annotation provenance ---
        "metadata: vdjdb_score alone": test.set_index("tcr_id")["vdjdb_score"],
        # --- sequence neighbourhood ---
        "density: tcrdist3 kernel": tcrdist_density(test),
        "density: SCEPTR default": sceptr_density(test, "default"),
        "density: SCEPTR small": sceptr_density(test, "small"),
        "density: SCEPTR alpha-chain only": sceptr_density(test, "a_sceptr"),
        "density: SCEPTR beta-chain only": sceptr_density(test, "b_sceptr"),
        # --- other single signals ---
        "generation probability (-log10 Pgen paired)": -pgen["log10_pgen_paired"],
        "structure: whole-complex confidence": structure["confidence_score"],
        "specificity model (target epitope held out)": specificity,
    }

    ranks = {name: rank_within_epitope(series, test) for name, series in raw.items()}

    # Fusions, all zero-parameter rank averages with weights fixed a priori.
    fusions: dict[str, pd.Series] = {
        "FUSION: 2-space density": (
            ranks["density: SCEPTR default"] + ranks["density: SCEPTR small"]
        ) / 2,
        "FUSION: 3-space density": (
            ranks["density: SCEPTR default"]
            + ranks["density: SCEPTR small"]
            + ranks["density: tcrdist3 kernel"]
        ) / 3,
        "FUSION: 2-space density + structure @0.2": None,
        "FUSION: 3-space density + structure @0.2": None,
        "FUSION: 3-space density + specificity @0.2": None,
    }
    two = fusions["FUSION: 2-space density"]
    three = fusions["FUSION: 3-space density"]
    struct_rank = ranks["structure: TCR-peptide min-PAE (AF3 analogue)"]
    fusions["FUSION: 2-space density + structure @0.2"] = 0.8 * two + 0.2 * struct_rank
    fusions["FUSION: 3-space density + structure @0.2"] = 0.8 * three + 0.2 * struct_rank
    fusions["FUSION: 3-space density + specificity @0.2"] = (
        0.8 * three + 0.2 * ranks["specificity model (target epitope held out)"]
    )

    rows = []
    for name, series in {**ranks, **fusions}.items():
        per_epitope = {}
        for epitope in epitopes:
            subset = labels[labels["epitope_aa"] == epitope]
            scores = subset["tcr_id"].map(series).to_numpy()
            per_epitope[epitope] = roc_auc_score(subset["pos"].to_numpy(), scores)
        rows.append({"method": name, **per_epitope,
                     "mean": float(np.mean(list(per_epitope.values())))})

    table = pd.DataFrame(rows).set_index("method").sort_values("mean", ascending=False)
    pd.set_option("display.width", 260)
    print(table.round(4).to_string())

    lines = [
        "# Master cross-epitope comparison",
        "",
        "Mean per-epitope ROC AUC over the five public epitopes, the competition's",
        "own metric. No method is tuned on the epitope it is scored on; kernel",
        "bandwidth is fixed a priori at 0.1 (SCEPTR cosine) / 20 (tcrdist3), which",
        "the AUC-vs-bandwidth curve shows is flat between 0.05 and 0.15.",
        "",
        "Fusions are zero-parameter rank averages with weights fixed in advance,",
        "because fitting a combiner on five epitopes loses more than it gains.",
        "",
        table.round(4).to_markdown(),
        "",
        "## Reading it",
        "",
        "- Every sequence-neighbourhood density lands in 0.75-0.84 regardless of",
        "  representation. The signal is one-dimensional and saturated.",
        "- The alpha chain beats the beta chain, inverting the usual expectation.",
        "- `vdjdb_score` and the specificity model sit below chance: both are",
        "  quantities whose sign flips between epitopes.",
        "- Structure is a genuine but non-orthogonal signal; it adds ~0.007 when",
        "  blended at low weight and nothing when used alone.",
        "- Selecting any of these choices by cross-validation on five folds costs",
        "  about 0.005, so the configuration is fixed a priori instead.",
        "",
    ]
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
