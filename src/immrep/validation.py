"""Leave-one-epitope-out validation.

Random splits would leak badly here: validated TCRs sit in tight clusters, so a
near-identical clonotype on the other side of a random split hands the model
the answer. Holding out whole epitopes is also the only split that matches what
the private leaderboard actually asks -- generalise to an epitope never seen in
training.

Five folds from five public epitopes is a small, noisy experiment. The standard
error of the mean per-epitope AUC is about 0.011, so this harness can separate
approaches that differ by ~0.03 and cannot be used to tune hyperparameters
finely. Per-fold numbers are always returned, never just the mean: a single
average hides an epitope that collapsed, which is the characteristic failure
of extrapolating to a new epitope.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# A model is fitted on the training epitopes and must return one score per row
# of the held-out frame. Higher means more likely to validate.
Fitter = Callable[[pd.DataFrame, np.ndarray, pd.DataFrame], np.ndarray]


def leave_one_epitope_out(
    features: pd.DataFrame,
    labels: pd.Series | np.ndarray,
    epitopes: pd.Series,
    fitter: Fitter,
) -> pd.DataFrame:
    """Fit on all epitopes but one, score the held-out one, repeat.

    Returns one row per fold with its AUC, plus a `mean` row.
    """
    labels = np.asarray(labels)
    epitopes = pd.Series(epitopes).reset_index(drop=True)
    features = features.reset_index(drop=True)

    rows = []
    for held_out in sorted(epitopes.unique()):
        test_mask = (epitopes == held_out).to_numpy()
        train_mask = ~test_mask

        predictions = fitter(
            features.loc[train_mask], labels[train_mask], features.loc[test_mask]
        )
        y = labels[test_mask]
        rows.append(
            {
                "held_out": held_out,
                "n": int(test_mask.sum()),
                "pos_rate": float(y.mean()),
                "auc": float(roc_auc_score(y, predictions))
                if len(np.unique(y)) == 2
                else float("nan"),
            }
        )

    report = pd.DataFrame(rows).set_index("held_out")
    report.loc["mean"] = {
        "n": report["n"].sum(),
        "pos_rate": float("nan"),
        "auc": report["auc"].mean(),
    }
    return report


def single_feature_fitter(column: str, *, sign: float = 1.0) -> Fitter:
    """Use one feature directly, ignoring the training labels.

    The reference any fitted model has to beat: the unsupervised density
    features already reach 0.84, and a model that cannot clear that is adding
    variance rather than information.

    Cluster-shape features are NaN for singletons, so missing values are filled
    from the training fold's median -- the same treatment the fitted models get,
    which keeps the comparison honest and avoids reading the held-out cohort.
    """

    def fit(train_x: pd.DataFrame, _train_y: np.ndarray, test_x: pd.DataFrame) -> np.ndarray:
        fallback = train_x[column].median()
        values = test_x[column].fillna(fallback).to_numpy(dtype=float)
        return sign * np.nan_to_num(values, nan=0.0)

    return fit


def sklearn_fitter(make_model: Callable[[], object], columns: list[str] | None = None) -> Fitter:
    """Fit any sklearn-style classifier, imputing NaNs from the training fold.

    Cluster-shape features are NaN for singletons by construction, so the
    imputation value must come from the training epitopes only -- taking it from
    the full frame would leak the held-out cohort's distribution.
    """

    def fit(train_x: pd.DataFrame, train_y: np.ndarray, test_x: pd.DataFrame) -> np.ndarray:
        used = columns or [c for c in train_x.columns if c not in ("tcr_id", "epitope_aa")]
        medians = train_x[used].median()
        model = make_model()
        model.fit(train_x[used].fillna(medians).to_numpy(dtype=float), train_y)
        filled = test_x[used].fillna(medians).to_numpy(dtype=float)
        if hasattr(model, "predict_proba"):
            return model.predict_proba(filled)[:, 1]
        return model.decision_function(filled)

    return fit
