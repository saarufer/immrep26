"""Building and validating submission files.

The competition accepts any finite real score per `tcr_id`; only the ranking
within each epitope is read. Rows must cover all of test.csv exactly once.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .constants import N_TEST_ROWS
from .data import load_test

ROOT = Path(__file__).resolve().parents[2]
SUBMISSION_DIR = ROOT / "submissions"


class SubmissionError(ValueError):
    """Raised when a submission would be rejected or silently mis-scored."""


def validate(submission: pd.DataFrame, *, test: pd.DataFrame | None = None) -> None:
    """Fail loudly on anything that would cost a submission slot.

    Checks column names, row count, id coverage and ordering, and that every
    prediction is finite. Ordering is checked because test.csv order is the
    documented requirement, even though scoring joins on `tcr_id`.
    """
    test = load_test() if test is None else test

    if list(submission.columns)[:2] != ["tcr_id", "prediction"]:
        raise SubmissionError(
            f"columns must start with ['tcr_id', 'prediction'], got {list(submission.columns)}"
        )
    if len(submission) != N_TEST_ROWS:
        raise SubmissionError(f"expected {N_TEST_ROWS} rows, got {len(submission)}")
    if submission["tcr_id"].duplicated().any():
        dupes = submission.loc[submission["tcr_id"].duplicated(), "tcr_id"].tolist()
        raise SubmissionError(f"duplicate tcr_id: {dupes[:5]}")

    missing = set(test["tcr_id"]) - set(submission["tcr_id"])
    extra = set(submission["tcr_id"]) - set(test["tcr_id"])
    if missing or extra:
        raise SubmissionError(
            f"id mismatch: {len(missing)} missing e.g. {sorted(missing)[:3]}, "
            f"{len(extra)} unexpected e.g. {sorted(extra)[:3]}"
        )
    if not submission["tcr_id"].tolist() == test["tcr_id"].tolist():
        raise SubmissionError("row order differs from test.csv")

    predictions = submission["prediction"].to_numpy(dtype=float)
    if not np.isfinite(predictions).all():
        bad = int((~np.isfinite(predictions)).sum())
        raise SubmissionError(f"{bad} non-finite predictions (NaN/inf)")

    # A constant score inside an epitope makes that epitope's AUC 0.5, which is
    # a silent third of the private score rather than an error.
    joined = test[["tcr_id", "epitope_aa"]].merge(submission, on="tcr_id")
    flat = [
        epitope
        for epitope, group in joined.groupby("epitope_aa")
        if group["prediction"].nunique() <= 1
    ]
    if flat:
        raise SubmissionError(f"prediction is constant within epitope(s): {flat}")


def build(
    scores: pd.Series | dict[str, float],
    name: str,
    *,
    test: pd.DataFrame | None = None,
    rank_within_epitope: bool = True,
) -> Path:
    """Write a validated submission to submissions/<name>.csv.

    `scores` maps tcr_id -> score. Because the metric only reads within-epitope
    ranking, scores are converted to within-epitope ranks by default: it cannot
    change the score, and it keeps one epitope's scale from being mistaken for
    signal if the file is ever pooled or blended with another.
    """
    test = load_test() if test is None else test
    series = pd.Series(scores, dtype=float) if isinstance(scores, dict) else scores.astype(float)

    submission = test[["tcr_id", "epitope_aa"]].copy()
    submission["prediction"] = submission["tcr_id"].map(series)
    unmapped = int(submission["prediction"].isna().sum())
    if unmapped:
        raise SubmissionError(f"{unmapped} test rows have no score")

    if rank_within_epitope:
        submission["prediction"] = submission.groupby("epitope_aa")["prediction"].rank(
            method="average", pct=True
        )

    submission = submission[["tcr_id", "prediction"]]
    validate(submission, test=test)

    SUBMISSION_DIR.mkdir(parents=True, exist_ok=True)
    path = SUBMISSION_DIR / f"{name}.csv"
    submission.to_csv(path, index=False)
    return path
