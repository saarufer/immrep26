"""The official IMMREP26 metric, reverse-engineered and confirmed.

The competition ships `tcrdist_predictions.csv`, whose public-leaderboard score
is 0.824. Scoring that file against the published TCRvdb labels under every
plausible metric gives exactly one match:

    mean per-epitope ROC AUC, ambiguous rows dropped   0.824728  -> 0.824
    pooled ROC AUC                                     0.857734
    pooled average precision                           0.860230
    mean per-epitope average precision                 0.833810
    n-weighted mean per-epitope ROC AUC                0.846627
    mean per-epitope AUC, ambiguous as negative        0.776939
    mean per-epitope AUC, ambiguous as positive        0.790829

Two consequences drive every modelling decision downstream:

1. Only the ranking *within* an epitope matters. Scores are never compared
   across epitopes, so the 28.7%-66.3% spread in validation rate between
   epitopes needs no calibration, and rank-normalising per epitope is free.
2. Every epitope carries equal weight regardless of size, so the 130-row
   LTDEMIAQY cohort counts as much as the 243-row CINGVCWTV one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# Reproduced with scripts/02_confirm_metric.py.
TCRDIST3_PUBLIC_SCORE = 0.824


def per_epitope_auc(
    labels: pd.DataFrame,
    predictions: pd.Series | np.ndarray,
    *,
    epitope_col: str = "epitope_aa",
    label_col: str = "pos",
) -> pd.Series:
    """ROC AUC within each epitope. Rows must already exclude ambiguous calls.

    Epitopes with a single class present score NaN rather than raising, which
    happens in leave-one-epitope-out folds built from tiny subsets.
    """
    frame = labels[[epitope_col, label_col]].copy()
    frame["_pred"] = np.asarray(predictions, dtype=float)

    out: dict[str, float] = {}
    for epitope, group in frame.groupby(epitope_col, sort=True):
        y = group[label_col].to_numpy()
        out[epitope] = (
            float(roc_auc_score(y, group["_pred"].to_numpy()))
            if len(np.unique(y)) == 2
            else float("nan")
        )
    return pd.Series(out, name="auc")


def official_score(
    labels: pd.DataFrame,
    predictions: pd.Series | np.ndarray,
    *,
    epitope_col: str = "epitope_aa",
    label_col: str = "pos",
) -> float:
    """The leaderboard number: unweighted mean of the per-epitope AUCs."""
    return float(
        per_epitope_auc(
            labels, predictions, epitope_col=epitope_col, label_col=label_col
        ).mean()
    )


def hanley_mcneil_se(auc: float, n_pos: int, n_neg: int) -> float:
    """Standard error of an AUC estimate (Hanley & McNeil 1982).

    Used to set expectations rather than to test anything: with roughly 213,
    201 and 114 scored rows on the three private epitopes, the standard error
    of their mean is about 0.023, so private-leaderboard gaps much below 0.05
    are not distinguishable from noise.
    """
    if n_pos < 1 or n_neg < 1:
        return float("nan")
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (
        auc * (1.0 - auc)
        + (n_pos - 1) * (q1 - auc * auc)
        + (n_neg - 1) * (q2 - auc * auc)
    ) / (n_pos * n_neg)
    return float(np.sqrt(max(var, 0.0)))


def score_report(
    labels: pd.DataFrame,
    predictions: pd.Series | np.ndarray,
    *,
    epitope_col: str = "epitope_aa",
    label_col: str = "pos",
) -> pd.DataFrame:
    """Per-epitope AUC with its standard error, plus a `mean` row.

    Always report this rather than the scalar score: a single mean hides an
    epitope that collapsed, which is the main failure mode when extrapolating
    to unseen epitopes.
    """
    aucs = per_epitope_auc(
        labels, predictions, epitope_col=epitope_col, label_col=label_col
    )
    rows = []
    for epitope, auc in aucs.items():
        group = labels[labels[epitope_col] == epitope]
        n_pos = int((group[label_col] == 1).sum())
        n_neg = int((group[label_col] == 0).sum())
        rows.append(
            {
                "epitope": epitope,
                "n": n_pos + n_neg,
                "n_pos": n_pos,
                "pos_rate": n_pos / (n_pos + n_neg) if n_pos + n_neg else float("nan"),
                "auc": auc,
                "se": hanley_mcneil_se(auc, n_pos, n_neg),
            }
        )
    report = pd.DataFrame(rows).set_index("epitope")
    report.loc["mean"] = {
        "n": report["n"].sum(),
        "n_pos": report["n_pos"].sum(),
        "pos_rate": float("nan"),
        "auc": report["auc"].mean(),
        # Independent epitopes, so the mean's SE shrinks by 1/k.
        "se": float(np.sqrt((report["se"] ** 2).sum()) / len(report)),
    }
    return report
