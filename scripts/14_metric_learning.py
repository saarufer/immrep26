#!/usr/bin/env python
"""Learn the metric that density is computed in, under strict leave-one-epitope-out.

Fitting a classifier on top of density already failed (scripts/07): with five
epitopes a combiner loses more to variance than it gains. But that attempt left
the *space* fixed and learned weights over summary features. The signal in this
task is "is this TCR inside a convergent cluster", so the thing worth learning
is the geometry that decides what counts as a cluster.

Parameterisation is deliberately minimal, because four training epitopes cannot
support much:

  scale-only   one parameter, the kernel bandwidth sigma.
  diagonal     64 parameters, one per SCEPTR dimension, shrunk towards the
               identity so the learned metric degrades to plain cosine when the
               data does not support a deviation.

The objective is a direct differentiable surrogate for the competition metric:
for every positive-negative pair *within* an epitope, penalise the pair being
ordered wrongly. That is within-epitope pairwise ranking, which is what mean
per-epitope AUC measures, rather than a classification loss that would also
chase the between-epitope offsets the metric ignores.

Evaluation is strictly nested: the metric is fitted on four epitopes -- with its
own inner leave-one-epitope-out for early stopping and shrinkage strength -- and
only then applied to the fifth.

Run: python scripts/14_metric_learning.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
SPACE = "default"
MAX_STEPS = 400
SEED = 0


def cohort_tensors(test: pd.DataFrame, labels: pd.DataFrame) -> dict[str, dict]:
    """Per epitope: the labelled rows' vectors, their labels, and the full
    cohort's vectors.

    Density must be computed against the *whole* cohort, including the 252
    ambiguous rows, because at inference time on a private epitope there is no
    way to know which rows would have been ambiguous. Only the labelled rows
    contribute to the loss.
    """
    cached = np.load(ROOT / f"data/interim/sceptr/{SPACE}.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    vectors = cached["vectors"]

    label_of = dict(zip(labels["tcr_id"], labels["pos"]))
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].tolist()
        rows = [position[t] for t in ids]
        block = torch.tensor(vectors[rows], dtype=torch.float32, device=DEVICE)
        scored = np.array([t in label_of for t in ids])
        out[epitope] = {
            "cohort": block,
            "scored": torch.tensor(scored, device=DEVICE),
            "y": torch.tensor(
                [float(label_of[t]) for t in ids if t in label_of],
                dtype=torch.float32, device=DEVICE,
            ),
            "ids": ids,
        }
    return out


def density(cohort: torch.Tensor, log_scale: torch.Tensor, log_sigma: torch.Tensor) -> torch.Tensor:
    """Kernel density inside one cohort under a diagonal-scaled cosine metric."""
    scaled = cohort * torch.exp(log_scale)
    scaled = scaled / scaled.norm(dim=1, keepdim=True).clamp_min(1e-8)
    distances = 1.0 - scaled @ scaled.T
    distances = distances + torch.eye(len(cohort), device=cohort.device) * 1e9
    return torch.exp(-distances / torch.exp(log_sigma)).mean(dim=1)


def pairwise_rank_loss(scores: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """Softplus surrogate for within-epitope AUC."""
    positives, negatives = scores[y > 0.5], scores[y < 0.5]
    if len(positives) == 0 or len(negatives) == 0:
        return torch.zeros((), device=scores.device)
    margin = positives.unsqueeze(1) - negatives.unsqueeze(0)
    # Scores are O(0.1); scale the margin so softplus is in its informative range.
    return torch.nn.functional.softplus(-margin * 50.0).mean()


def fit(
    data: dict, train_epitopes: list[str], *, diagonal: bool, shrinkage: float, steps: int
) -> tuple[torch.Tensor, torch.Tensor]:
    torch.manual_seed(SEED)
    dim = data[train_epitopes[0]]["cohort"].shape[1]
    log_scale = torch.zeros(dim, device=DEVICE, requires_grad=diagonal)
    log_sigma = torch.tensor(np.log(0.1), device=DEVICE, requires_grad=True)

    params = [log_sigma] + ([log_scale] if diagonal else [])
    optimiser = torch.optim.Adam(params, lr=0.02)

    for _ in range(steps):
        optimiser.zero_grad()
        loss = torch.zeros((), device=DEVICE)
        for epitope in train_epitopes:
            entry = data[epitope]
            scores = density(entry["cohort"], log_scale, log_sigma)
            loss = loss + pairwise_rank_loss(scores[entry["scored"]], entry["y"])
        loss = loss / len(train_epitopes)
        if diagonal:
            loss = loss + shrinkage * (log_scale**2).mean()
        loss.backward()
        optimiser.step()
    return log_scale.detach(), log_sigma.detach()


def evaluate(data: dict, epitope: str, log_scale: torch.Tensor, log_sigma: torch.Tensor) -> float:
    entry = data[epitope]
    with torch.no_grad():
        scores = density(entry["cohort"], log_scale, log_sigma)[entry["scored"]]
    return roc_auc_score(entry["y"].cpu().numpy(), scores.cpu().numpy())


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    data = cohort_tensors(test, labels)
    epitopes = sorted(labels["epitope_aa"].unique())

    print(f"device={DEVICE}  space={SPACE}\n")

    baseline = {}
    for epitope in epitopes:
        zero = torch.zeros(data[epitope]["cohort"].shape[1], device=DEVICE)
        baseline[epitope] = evaluate(
            data, epitope, zero, torch.tensor(np.log(0.1), device=DEVICE)
        )
    print("=== reference: plain cosine density, sigma=0.1 (no learning) ===")
    print("  " + "  ".join(f"{e[:5]}={baseline[e]:.3f}" for e in epitopes))
    print(f"  mean = {np.mean(list(baseline.values())):.4f}\n")

    rows = []
    for variant, diagonal in (("scale-only (1 param)", False), ("diagonal (64 params)", True)):
        fold_aucs = {}
        for held_out in epitopes:
            train = [e for e in epitopes if e != held_out]
            if diagonal:
                # Inner LOEO over the four training epitopes picks the shrinkage,
                # so the held-out epitope never influences any choice.
                best, best_shrinkage = -1.0, None
                for shrinkage in (0.03, 0.1, 0.3, 1.0):
                    inner = []
                    for inner_held in train:
                        inner_train = [e for e in train if e != inner_held]
                        s, g = fit(
                            data, inner_train, diagonal=True,
                            shrinkage=shrinkage, steps=MAX_STEPS,
                        )
                        inner.append(evaluate(data, inner_held, s, g))
                    if np.mean(inner) > best:
                        best, best_shrinkage = float(np.mean(inner)), shrinkage
                scale, sigma = fit(
                    data, train, diagonal=True, shrinkage=best_shrinkage, steps=MAX_STEPS
                )
                note = f"shrink={best_shrinkage}"
            else:
                scale, sigma = fit(data, train, diagonal=False, shrinkage=0.0, steps=MAX_STEPS)
                note = f"sigma={float(torch.exp(sigma)):.3f}"
            fold_aucs[held_out] = evaluate(data, held_out, scale, sigma)
            print(f"  {variant:<22} hold out {held_out:<11} auc={fold_aucs[held_out]:.4f}  {note}")
        rows.append({"variant": variant, **fold_aucs, "mean": float(np.mean(list(fold_aucs.values())))})
        print()

    table = pd.DataFrame(rows).set_index("variant")
    table.loc["plain cosine, sigma=0.1"] = {**baseline, "mean": np.mean(list(baseline.values()))}
    pd.set_option("display.width", 240)
    print("=== strict nested LOEO ===")
    print(table.round(4).to_string())
    best = table["mean"].idxmax()
    print(f"\nbest: {best} at {table.loc[best, 'mean']:.4f}")
    print(f"learning gain over plain cosine: "
          f"{table['mean'].max() - table.loc['plain cosine, sigma=0.1', 'mean']:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
