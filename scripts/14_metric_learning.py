#!/usr/bin/env python
"""Learn the metric that density is computed in, under strict leave-one-epitope-out.

Fitting a classifier on top of density already failed (scripts/07): with five
epitopes a combiner loses more to variance than it gains. That attempt left the
*space* fixed and learned weights over summary features. Since the signal in
this task is "is this TCR inside a convergent cluster", the thing worth
learning is the geometry that decides what counts as a cluster.

Parameterisation stays minimal, because four training epitopes support little:

  scale-only   one parameter, the kernel bandwidth sigma.
  diagonal     64 parameters, one per SCEPTR dimension, shrunk towards the
               identity so the metric degrades to plain cosine when the data
               does not support a deviation.

## The surrogate has to be scale-free

A first version of this script minimised `softplus(-margin * 50)` over
positive-negative pairs and drove sigma to 0.22-0.40, where held-out AUC is
0.8167, when the true optimum is 0.07-0.10 at 0.8397. The cause is that
softplus penalises negative margins roughly linearly, so the loss is dominated
by its worst pairs. Small sigma makes density spiky -- nearest-neighbour
dominated -- which produces a few badly ordered pairs with large negative
margins, while large sigma pushes every density towards 1 and every margin
towards 0, parking the loss at log 2. Flattening the scores was cheaper than
ordering them, so the optimiser flattened them.

The fix is to standardise density *within each epitope* before the pairwise
loss. Then shrinking the spread cannot reduce the loss, and only the ordering
can. Standardising is also the right thing on its own terms: the competition
metric reads nothing but within-epitope order.

Verification criterion: with the fix, scale-only learning should recover
sigma near 0.07-0.10 and score near 0.8397. If it does not, the surrogate is
still wrong and the diagonal result cannot be interpreted.

Evaluation is strictly nested: the metric is fitted on four epitopes -- with an
inner leave-one-epitope-out for the shrinkage strength -- and only then applied
to the fifth.

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
STEPS = 500
SEED = 0
INIT_SIGMA = 0.1


def cohort_tensors(test: pd.DataFrame, labels: pd.DataFrame) -> dict[str, dict]:
    """Per epitope: full-cohort vectors, a mask of labelled rows, their labels.

    Density is computed against the *whole* cohort including the ambiguous rows,
    because on a private epitope there is no way to know which rows would have
    been ambiguous. Only labelled rows enter the loss.
    """
    cached = np.load(ROOT / f"data/interim/sceptr/{SPACE}.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    vectors = cached["vectors"]
    label_of = dict(zip(labels["tcr_id"], labels["pos"]))

    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].tolist()
        block = torch.tensor(
            vectors[[position[t] for t in ids]], dtype=torch.float32, device=DEVICE
        )
        scored = np.array([t in label_of for t in ids])
        out[epitope] = {
            "cohort": block,
            "scored": torch.tensor(scored, device=DEVICE),
            "y": torch.tensor(
                [float(label_of[t]) for t in ids if t in label_of],
                dtype=torch.float32, device=DEVICE,
            ),
        }
    return out


def density(cohort: torch.Tensor, log_scale: torch.Tensor, log_sigma: torch.Tensor) -> torch.Tensor:
    scaled = cohort * torch.exp(log_scale)
    scaled = scaled / scaled.norm(dim=1, keepdim=True).clamp_min(1e-8)
    distances = 1.0 - scaled @ scaled.T
    distances = distances + torch.eye(len(cohort), device=cohort.device) * 1e9
    return torch.exp(-distances / torch.exp(log_sigma)).mean(dim=1)


def ranking_loss(scores: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """Scale-free within-epitope pairwise ranking loss.

    Standardising `scores` first is what stops the optimiser from reducing the
    loss by flattening the density distribution instead of ordering it.
    """
    z = (scores - scores.mean()) / scores.std().clamp_min(1e-6)
    positives, negatives = z[y > 0.5], z[y < 0.5]
    if len(positives) == 0 or len(negatives) == 0:
        return torch.zeros((), device=scores.device)
    margin = positives.unsqueeze(1) - negatives.unsqueeze(0)
    return torch.nn.functional.softplus(-margin).mean()


def fit(data, train_epitopes, *, diagonal: bool, shrinkage: float):
    torch.manual_seed(SEED)
    dim = data[train_epitopes[0]]["cohort"].shape[1]
    log_scale = torch.zeros(dim, device=DEVICE, requires_grad=diagonal)
    log_sigma = torch.tensor(float(np.log(INIT_SIGMA)), device=DEVICE, requires_grad=True)

    optimiser = torch.optim.Adam([log_sigma] + ([log_scale] if diagonal else []), lr=0.02)
    for _ in range(STEPS):
        optimiser.zero_grad()
        loss = torch.zeros((), device=DEVICE)
        for epitope in train_epitopes:
            entry = data[epitope]
            scores = density(entry["cohort"], log_scale, log_sigma)
            loss = loss + ranking_loss(scores[entry["scored"]], entry["y"])
        loss = loss / len(train_epitopes)
        if diagonal:
            loss = loss + shrinkage * (log_scale**2).mean()
        loss.backward()
        optimiser.step()
    return log_scale.detach(), log_sigma.detach()


def evaluate(data, epitope, log_scale, log_sigma) -> float:
    entry = data[epitope]
    with torch.no_grad():
        scores = density(entry["cohort"], log_scale, log_sigma)[entry["scored"]]
    return roc_auc_score(entry["y"].cpu().numpy(), scores.cpu().numpy())


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    data = cohort_tensors(test, labels)
    epitopes = sorted(labels["epitope_aa"].unique())
    zeros = torch.zeros(data[epitopes[0]]["cohort"].shape[1], device=DEVICE)

    print(f"device={DEVICE}  space={SPACE}  steps={STEPS}\n")

    baseline = {
        e: evaluate(data, e, zeros, torch.tensor(float(np.log(0.1)), device=DEVICE))
        for e in epitopes
    }
    print("=== reference: plain cosine density, sigma=0.1 (nothing learned) ===")
    print("  " + "  ".join(f"{e[:5]}={baseline[e]:.3f}" for e in epitopes))
    print(f"  mean = {np.mean(list(baseline.values())):.4f}\n")

    rows = []
    for variant, diagonal in (("scale-only (1 param)", False), ("diagonal (64 params)", True)):
        fold = {}
        for held_out in epitopes:
            train = [e for e in epitopes if e != held_out]
            if diagonal:
                best, best_shrinkage = -1.0, None
                for shrinkage in (0.03, 0.1, 0.3, 1.0):
                    inner = [
                        evaluate(
                            data, inner_held,
                            *fit(data, [e for e in train if e != inner_held],
                                 diagonal=True, shrinkage=shrinkage),
                        )
                        for inner_held in train
                    ]
                    if np.mean(inner) > best:
                        best, best_shrinkage = float(np.mean(inner)), shrinkage
                scale, sigma = fit(data, train, diagonal=True, shrinkage=best_shrinkage)
                note = f"shrink={best_shrinkage}  sigma={float(torch.exp(sigma)):.3f}"
            else:
                scale, sigma = fit(data, train, diagonal=False, shrinkage=0.0)
                note = f"sigma={float(torch.exp(sigma)):.3f}"
            fold[held_out] = evaluate(data, held_out, scale, sigma)
            print(f"  {variant:<22} hold out {held_out:<11} auc={fold[held_out]:.4f}  {note}")
        rows.append({"variant": variant, **fold, "mean": float(np.mean(list(fold.values())))})
        print()

    table = pd.DataFrame(rows).set_index("variant")
    table.loc["plain cosine, sigma=0.1"] = {
        **baseline, "mean": float(np.mean(list(baseline.values())))
    }
    pd.set_option("display.width", 240)
    print("=== strict nested LOEO ===")
    print(table.round(4).to_string())

    reference = table.loc["plain cosine, sigma=0.1", "mean"]
    scale_only = table.loc["scale-only (1 param)", "mean"]
    print(f"\nsurrogate check: scale-only should land near {reference:.4f} with sigma near "
          f"0.07-0.10 if the loss is now aligned with AUC.")
    print(f"  scale-only = {scale_only:.4f}  ({scale_only - reference:+.4f} vs plain cosine)")
    best = table["mean"].idxmax()
    print(f"\nbest: {best} at {table.loc[best, 'mean']:.4f}  "
          f"(learning gain {table['mean'].max() - reference:+.4f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
