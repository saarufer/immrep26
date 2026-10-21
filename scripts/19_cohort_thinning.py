#!/usr/bin/env python
"""Does density survive a sparser cohort? A robustness test aimed at the private split.

Density is computed *within* each epitope cohort, so its quality depends on how
well that cohort covers the relevant TCR space. The private cohorts are both
smaller and markedly sparser than the public ones: the fraction of TCRs with a
neighbour inside the clustering radius is 2.5% for CINGVCWTV and 4.4% for
LLWNGPMAV, against 4.7-55% across the public five. CINGVCWTV is therefore
sparser than anything the method has been validated on, which makes the
leave-one-epitope-out estimate an extrapolation rather than an interpolation.

Thinning the public cohorts simulates that regime directly, with real labels.

Two variants, because the naive version confounds two effects:

  thin_both       drop TCRs from the cohort and score only the survivors.
                  Mixes a genuinely worse density estimate with a noisier AUC
                  estimate computed on fewer rows.
  thin_reference  drop TCRs from the *reference set* only, but still score every
                  labelled row. Isolates the effect of reference coverage, which
                  is the thing in question.

The payoff is the second table: AUC against the close-neighbour fraction that
thinning actually produced, which can be read at each private epitope's observed
sparsity to get a calibrated expectation instead of an assumption.

Run: python scripts/19_cohort_thinning.py
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

from immrep.data import build_training_frame, load_test  # noqa: E402

SIGMA = 0.1
CLUSTER_RADIUS = 0.12
FRACTIONS = (1.0, 0.8, 0.6, 0.45, 0.3, 0.2, 0.12)
REPEATS = 40
STRUCTURE_WEIGHT = 0.2
SEED = 0
OUT = ROOT / "reports/05_cohort_thinning.md"

PRIVATE_SPARSITY = {"CINGVCWTV": 0.025, "LLWNGPMAV": 0.044, "LTDEMIAQY": 0.146}


def density_from(block: np.ndarray, reference: np.ndarray, sigma: float) -> np.ndarray:
    """Kernel density of each row of `block` against `reference`.

    Self-matches are removed by dropping the zero-distance entry per row, which
    is what the within-cohort version does by masking the diagonal.
    """
    distances = 1.0 - block @ reference.T
    # Each scored row appears in the reference, so its own zero distance is the
    # row minimum; sorting and dropping column 0 removes exactly one self-match.
    distances = np.sort(distances, axis=1)[:, 1:]
    return np.exp(-distances / sigma).mean(axis=1)


def close_neighbour_fraction(block: np.ndarray, reference: np.ndarray) -> float:
    distances = 1.0 - block @ reference.T
    nearest = np.sort(distances, axis=1)[:, 1]
    return float((nearest <= CLUSTER_RADIUS).mean())


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    vectors = cached["vectors"]
    normed = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")

    rng = np.random.default_rng(SEED)
    rows = []

    for epitope, cohort in test.groupby("epitope_aa"):
        labelled = labels[labels["epitope_aa"] == epitope]
        if labelled.empty:
            continue
        cohort_ids = cohort["tcr_id"].to_numpy()
        cohort_rows = normed[[position[t] for t in cohort_ids]]
        is_labelled = np.isin(cohort_ids, labelled["tcr_id"].to_numpy())
        label_lookup = dict(zip(labelled["tcr_id"], labelled["pos"]))
        y_all = np.array([label_lookup.get(t, -1) for t in cohort_ids])
        structure_rank = pd.Series(
            structure.loc[cohort_ids, "pae_tcr_pep_min"].to_numpy(), index=cohort_ids
        ).rank(pct=True).to_numpy()

        for fraction in FRACTIONS:
            keep_n = max(12, int(round(fraction * len(cohort_ids))))
            for repeat in range(1 if fraction == 1.0 else REPEATS):
                keep = (
                    np.ones(len(cohort_ids), dtype=bool)
                    if fraction == 1.0
                    else np.isin(np.arange(len(cohort_ids)), rng.choice(len(cohort_ids), keep_n, replace=False))
                )
                reference = cohort_rows[keep]

                for variant, scored_mask in (
                    ("thin_both", keep & is_labelled),
                    ("thin_reference", is_labelled),
                ):
                    y = y_all[scored_mask]
                    if len(np.unique(y)) < 2:
                        continue
                    block = cohort_rows[scored_mask]
                    density = density_from(block, reference, SIGMA)
                    # Rank within this cohort, then blend structure exactly as
                    # the final submission does.
                    density_rank = pd.Series(density).rank(pct=True).to_numpy()
                    blended = (
                        (1 - STRUCTURE_WEIGHT) * density_rank
                        + STRUCTURE_WEIGHT
                        * pd.Series(structure_rank[scored_mask]).rank(pct=True).to_numpy()
                    )
                    rows.append({
                        "epitope": epitope,
                        "variant": variant,
                        "fraction": fraction,
                        "repeat": repeat,
                        "n_reference": int(keep.sum()),
                        "n_scored": int(scored_mask.sum()),
                        "close_frac": close_neighbour_fraction(block, reference),
                        "auc_density": roc_auc_score(y, density),
                        "auc_blend": roc_auc_score(y, blended),
                    })

    frame = pd.DataFrame(rows)
    pd.set_option("display.width", 260)

    print("=== AUC vs thinning, averaged over repeats (density only) ===")
    pivot = (
        frame[frame.variant == "thin_reference"]
        .pivot_table(index="epitope", columns="fraction", values="auc_density", aggfunc="mean")
    )
    print(pivot.round(4).to_string())

    print("\n=== same, but thinning the scored rows too (the naive design) ===")
    pivot_both = (
        frame[frame.variant == "thin_both"]
        .pivot_table(index="epitope", columns="fraction", values="auc_density", aggfunc="mean")
    )
    print(pivot_both.round(4).to_string())

    print("\n=== mean over epitopes, both variants, density and blend ===")
    summary = frame.pivot_table(
        index="fraction", columns="variant", values=["auc_density", "auc_blend"], aggfunc="mean"
    )
    print(summary.round(4).to_string())

    print("\n=== spread across repeats (thin_reference, density) ===")
    spread = (
        frame[frame.variant == "thin_reference"]
        .pivot_table(index="fraction", values="auc_density", aggfunc=["mean", "std", "min"])
    )
    print(spread.round(4).to_string())

    # The calibration that matters: AUC as a function of achieved sparsity.
    print("\n=== calibration: AUC binned by the close-neighbour fraction thinning produced ===")
    sub = frame[frame.variant == "thin_reference"].copy()
    edges = [0, 0.02, 0.04, 0.07, 0.12, 0.20, 0.35, 1.01]
    sub["bin"] = pd.cut(sub["close_frac"], bins=edges)
    calibration = sub.groupby("bin", observed=True).agg(
        n=("auc_density", "size"),
        mean_close_frac=("close_frac", "mean"),
        auc_density=("auc_density", "mean"),
        auc_blend=("auc_blend", "mean"),
        auc_density_sd=("auc_density", "std"),
    )
    print(calibration.round(4).to_string())

    print("\n=== read-off for the private epitopes at their observed sparsity ===")
    for epitope, close in PRIVATE_SPARSITY.items():
        hit = calibration[
            [close in interval for interval in calibration.index]
        ]
        if len(hit):
            row = hit.iloc[0]
            print(
                f"  {epitope:<11} close_frac={close:.3f} -> expected density AUC "
                f"{row.auc_density:.3f} (sd {row.auc_density_sd:.3f}), blend {row.auc_blend:.3f}"
            )
        else:
            print(f"  {epitope:<11} close_frac={close:.3f} -> outside the thinned range")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Cohort thinning: does density survive a sparser repertoire?\n\n"
        "Public cohorts are randomly thinned to simulate the private split's lower\n"
        "candidate coverage, with real labels retained for scoring. `thin_reference`\n"
        "drops TCRs from the reference set only and still scores every labelled row,\n"
        "isolating coverage from the noise of scoring fewer rows.\n\n"
        "## AUC vs retained fraction (thin_reference, density only)\n\n"
        + pivot.round(4).to_markdown()
        + "\n\n## Mean over epitopes\n\n"
        + summary.round(4).to_markdown()
        + "\n\n## Calibration by achieved close-neighbour fraction\n\n"
        + calibration.round(4).to_markdown()
        + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
