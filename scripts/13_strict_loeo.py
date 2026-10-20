#!/usr/bin/env python
"""Strictly unbiased leave-one-epitope-out, including the cost of choosing a configuration.

The numbers reported earlier in this project are optimistic in a way that is
easy to miss. Kernel density is unsupervised, so no label touches the score --
but sigma, which representations to average, and which structure feature to
blend at what weight were all picked by looking at the public-epitope result.
That is selection on five epitopes, and it inflates the estimate exactly like
fitting would.

This script measures three things that are usually conflated:

  oracle        best fixed configuration, scored on all five epitopes.
                Upward-biased: this is the number you get by tuning on your
                evaluation set, and it is what "0.8415" was.

  nested LOEO   per fold, the configuration is re-selected using only the four
                training epitopes, then applied to the held-out one. Unbiased
                estimate of "tune by cross-validation, then submit".

  a priori      a configuration fixed on principled grounds before looking,
                evaluated by plain LOEO. Also unbiased, and with five folds it
                can beat nested selection, because selecting among many
                near-identical options on four epitopes adds more variance
                than it removes.

The gap between oracle and nested is the selection cost. If a priori beats
nested, the right move is to make fewer decisions, not better-cross-validated
ones.

Run: python scripts/13_strict_loeo.py
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

SPACES = ("default", "small", "large", "a_sceptr")
SIGMAS = (0.05, 0.1, 0.2)
STRUCTURE_FEATURES = ("pae_tcr_pep_min", "pae_a_mhc_min", "confidence_score", "iptm_tcr_pep_max")
STRUCTURE_WEIGHTS = (0.0, 0.1, 0.2, 0.3, 0.4)


def cohort_density(vectors: np.ndarray, sigma: float) -> np.ndarray:
    normed = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    distances = 1.0 - normed @ normed.T
    np.fill_diagonal(distances, np.inf)
    return np.exp(-distances / sigma).mean(axis=1)


def build_component_ranks(test: pd.DataFrame) -> pd.DataFrame:
    """Within-epitope percentile rank of every candidate component.

    Ranking per epitope is not a modelling choice -- the metric only reads
    within-epitope order -- so doing it up front costs nothing and makes every
    component commensurable for averaging.
    """
    frame = test[["tcr_id", "epitope_aa"]].copy()

    for space in SPACES:
        cached = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        vectors = cached["vectors"]
        for sigma in SIGMAS:
            values = pd.Series(index=test["tcr_id"], dtype=float)
            for _, cohort in test.groupby("epitope_aa"):
                rows = [position[t] for t in cohort["tcr_id"]]
                values.loc[cohort["tcr_id"]] = cohort_density(vectors[rows], sigma)
            frame[f"dens:{space}:{sigma}"] = frame["tcr_id"].map(values)

    # tcrdist3 lives on a different scale; its sigma grid is scaled to match.
    for sigma in (10, 20, 40):
        values = pd.Series(index=test["tcr_id"], dtype=float)
        for epitope, cohort in test.groupby("epitope_aa"):
            cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
            position = {t: i for i, t in enumerate(cached["tcr_id"])}
            rows = [position[t] for t in cohort["tcr_id"]]
            distances = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
            np.fill_diagonal(distances, np.inf)
            values.loc[cohort["tcr_id"]] = np.exp(-distances / sigma).mean(axis=1)
        frame[f"dens:tcrdist:{sigma}"] = frame["tcr_id"].map(values)

    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    for feature in STRUCTURE_FEATURES:
        frame[f"struct:{feature}"] = frame["tcr_id"].map(structure[feature])

    value_columns = [c for c in frame.columns if c not in ("tcr_id", "epitope_aa")]
    frame[value_columns] = frame.groupby("epitope_aa")[value_columns].rank(pct=True)
    return frame


def candidate_configs(columns: list[str]) -> list[tuple[tuple[str, ...], str, float]]:
    """(density components to average, structure feature, structure weight)."""
    density = [c for c in columns if c.startswith("dens:")]
    singles = [(c,) for c in density]
    # Pairs and triples of the best-known families, not the full power set --
    # the point is to measure selection cost, not to maximise the grid.
    pairs = [
        ("dens:default:0.1", "dens:small:0.1"),
        ("dens:default:0.1", "dens:tcrdist:20"),
        ("dens:default:0.1", "dens:a_sceptr:0.1"),
        ("dens:default:0.1", "dens:small:0.1", "dens:tcrdist:20"),
    ]
    configs = []
    for combo in singles + [p for p in pairs if all(c in columns for c in p)]:
        for feature in STRUCTURE_FEATURES:
            for weight in STRUCTURE_WEIGHTS:
                if weight == 0.0 and feature != STRUCTURE_FEATURES[0]:
                    continue  # weight 0 is the same config regardless of feature
                configs.append((combo, feature, weight))
    return configs


def score(frame: pd.DataFrame, config: tuple, rows: np.ndarray) -> np.ndarray:
    combo, feature, weight = config
    density = frame.loc[rows, list(combo)].mean(axis=1)
    if weight == 0.0:
        return density.to_numpy()
    return ((1 - weight) * density + weight * frame.loc[rows, f"struct:{feature}"]).to_numpy()


def mean_epitope_auc(frame, labels, config, epitopes: list[str]) -> float:
    aucs = []
    for epitope in epitopes:
        rows = (labels["epitope_aa"] == epitope).to_numpy()
        y = labels.loc[rows, "pos"].to_numpy()
        if len(np.unique(y)) < 2:
            continue
        aucs.append(roc_auc_score(y, score(frame, config, labels.index[rows])))
    return float(np.mean(aucs))


def main() -> int:
    test = load_test()
    ranks = build_component_ranks(test).set_index("tcr_id")

    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]].set_index("tcr_id")
    labels = labels.join(ranks.drop(columns=["epitope_aa"]))
    labels = labels.reset_index().set_index("tcr_id")

    frame = labels
    epitopes = sorted(labels["epitope_aa"].unique())
    configs = candidate_configs([c for c in ranks.columns if c != "epitope_aa"])
    print(f"{len(configs)} candidate configurations, {len(epitopes)} epitopes\n")

    def fold_auc(config, epitope) -> float:
        rows = (frame["epitope_aa"] == epitope).to_numpy()
        return roc_auc_score(
            frame.loc[rows, "pos"].to_numpy(), score(frame, config, frame.index[rows])
        )

    # --- oracle: tune on all five, report on all five (upward biased) ---
    scored = [(mean_epitope_auc(frame, frame, c, epitopes), c) for c in configs]
    oracle_value, oracle_config = max(scored, key=lambda t: t[0])

    # --- nested: choose on four, apply to the fifth ---
    nested_rows = []
    for held_out in epitopes:
        train = [e for e in epitopes if e != held_out]
        best = max(configs, key=lambda c: mean_epitope_auc(frame, frame, c, train))
        nested_rows.append(
            {
                "held_out": held_out,
                "chosen": f"{'+'.join(s.split(':')[1] for s in best[0])}"
                f" | {best[1] if best[2] else '-'}@{best[2]}",
                "auc": fold_auc(best, held_out),
            }
        )
    nested = pd.DataFrame(nested_rows)

    # --- a priori: fixed before looking, on stated grounds ---
    priors = {
        "density only (default, sigma=0.1)": (("dens:default:0.1",), STRUCTURE_FEATURES[0], 0.0),
        "two-space density average": (
            ("dens:default:0.1", "dens:small:0.1"), STRUCTURE_FEATURES[0], 0.0,
        ),
        "two-space + structure @0.2": (
            ("dens:default:0.1", "dens:small:0.1"), "pae_tcr_pep_min", 0.2,
        ),
        "three-space + structure @0.2": (
            ("dens:default:0.1", "dens:small:0.1", "dens:tcrdist:20"), "pae_tcr_pep_min", 0.2,
        ),
    }

    print("=== a priori configurations, plain LOEO (unbiased: nothing was selected) ===")
    prior_table = {}
    for name, config in priors.items():
        per_fold = {e: fold_auc(config, e) for e in epitopes}
        prior_table[name] = {**per_fold, "mean": float(np.mean(list(per_fold.values())))}
    prior_frame = pd.DataFrame(prior_table).T.sort_values("mean", ascending=False)
    pd.set_option("display.width", 240)
    print(prior_frame.round(4).to_string())

    print("\n=== nested LOEO: configuration re-chosen inside every fold ===")
    print(nested.round(4).to_string(index=False))
    print(f"  nested mean = {nested['auc'].mean():.4f}")

    print(f"\n=== the three numbers ===")
    best_prior_name = prior_frame.index[0]
    best_prior = prior_frame.loc[best_prior_name, "mean"]
    print(f"  oracle (tuned on the evaluation set)  {oracle_value:.4f}   <- biased, not achievable")
    print(f"  nested LOEO (tune by CV, then submit) {nested['auc'].mean():.4f}")
    print(f"  best a priori config, plain LOEO      {best_prior:.4f}   ({best_prior_name})")
    print(f"\n  selection cost (oracle - nested)      {oracle_value - nested['auc'].mean():+.4f}")
    print(f"  CV selection vs fixing a priori       {nested['auc'].mean() - best_prior:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
