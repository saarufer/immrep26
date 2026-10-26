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


# ---------------------------------------------------------------------------
# Standing gates. A feature passes both before its LOEO score is believed.


def private_variance_gate(test: pd.DataFrame, series: pd.Series) -> tuple[bool, str]:
    """E3's rule: a feature constant inside a private epitope contributes
    exactly zero where the ranking is decided, however well it scores publicly.

    Returns (passes, per-epitope level counts).
    """
    detail, ok = [], True
    for epitope, cohort in test[test["split"] == "private"].groupby("epitope_aa"):
        levels = cohort["tcr_id"].map(series).nunique()
        detail.append(f"{epitope[:5]}:{levels}")
        if levels <= 1:
            ok = False
    return ok, " ".join(detail)


CHARGE = {**{a: 1.0 for a in "KRH"}, **{a: -1.0 for a in "DE"}}
MASS = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103,
        "L": 113, "I": 113, "N": 114, "D": 115, "Q": 128, "K": 128, "E": 129,
        "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186}


def typicality(test: pd.DataFrame) -> pd.Series:
    """How ordinary a TCR's own chain features are, within its own cohort.

    This exists because it is the gravity well of this project. Four features
    that looked like different ideas all turned out to re-encode it: the
    neighbourhood density itself, Boltz interface confidence (E14), the
    whole-V-domain alignment distance (E26, E42) and alpha-beta pairing
    compatibility (E45) -- the last of which it outscored, 0.6378 to 0.6336,
    while involving no pairing at all.

    Deliberately crude: the summed negative z-score of CDR3 length, charge and
    mean residue mass on both chains, within epitope. A new feature that cannot
    beat this is not carrying new information, whatever it is called.
    """
    out = pd.Series(index=test["tcr_id"], dtype=float)
    columns = {}
    for side, cdr3 in (("a", "cdr3_alpha_aa"), ("b", "cdr3_beta_aa")):
        seqs = test[cdr3].astype(str)
        columns[f"len_{side}"] = seqs.str.len().to_numpy(dtype=float)
        columns[f"charge_{side}"] = np.array(
            [sum(CHARGE.get(c, 0.0) for c in s) for s in seqs])
        columns[f"mass_{side}"] = np.array(
            [np.mean([MASS.get(c, 110) for c in s]) for s in seqs])
    frame = pd.DataFrame(columns, index=test.index)
    for _, cohort in test.groupby("epitope_aa"):
        block = frame.loc[cohort.index]
        score = np.zeros(len(block))
        for column in block.columns:
            values = block[column].to_numpy()
            score -= np.abs((values - values.mean()) / (values.std() or 1.0))
        out.loc[cohort["tcr_id"]] = score
    return out


def typicality_gate(test: pd.DataFrame, labels: pd.DataFrame,
                    series: pd.Series) -> dict:
    """Is this feature carrying anything typicality does not already?

    Reports the feature's mean per-epitope AUC, typicality's on the same rows,
    and their within-epitope rank correlation. A feature that scores no higher
    than typicality is a fifth re-encoding of it, and the correlation says how
    directly.
    """
    from scipy.stats import spearmanr

    reference = typicality(test)

    def mean_auc(values: pd.Series) -> float:
        per = []
        for _, block in labels.groupby("epitope_aa"):
            mapped = block["tcr_id"].map(values).to_numpy()
            if len(np.unique(mapped)) > 1:
                per.append(roc_auc_score(block["pos"], mapped))
        return float(np.mean(per)) if per else float("nan")

    feature_auc, reference_auc = mean_auc(series), mean_auc(reference)
    rho = float(np.mean([
        spearmanr(cohort["tcr_id"].map(series),
                  cohort["tcr_id"].map(reference)).statistic
        for _, cohort in test.groupby("epitope_aa")]))
    beats = feature_auc > reference_auc
    return {"feature_auc": feature_auc, "typicality_auc": reference_auc,
            "rho_with_typicality": rho, "beats_typicality": beats,
            "verdict": ("carries something typicality does not"
                        if beats else
                        "does NOT beat typicality -- likely a re-encoding")}
