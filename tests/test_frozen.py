"""Guards on the frozen predictor and the invariants everything else assumes.

The predictor was frozen on 2026-10-25. These tests exist so that a later edit
to the metric, the data loader or the submission builder cannot silently change
what was frozen without something going red.

They are deliberately cheap: no GPU, no folding, no MD. The expensive stages
have their own self-checks inside their scripts.

TCRvdb is licensed and gitignored, so every test that needs labels skips rather
than fails when the data is absent. A skipped run still checks what it can.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep import constants as K  # noqa: E402

SUBMISSION = ROOT / "submissions/density3_struct03.csv"
FROZEN_SHA256 = "523a4692df02467b3d9b4fe0e51136a4aa4df3af7ddae2262cf8ed282ac5608c"
FROZEN_PUBLIC = 0.8498

needs_data = pytest.mark.skipif(
    not (ROOT / "data/external/tcrvdb").exists(),
    reason="TCRvdb is licensed and not committed",
)


def _labels():
    from immrep.data import build_training_frame

    return build_training_frame()


# --------------------------------------------------------------------------
# The frozen artifact


def test_frozen_submission_exists():
    assert SUBMISSION.exists(), "the frozen submission file is missing"


def test_frozen_submission_is_byte_identical():
    """If this fails the frozen file was edited, whatever the score says."""
    digest = hashlib.sha256(SUBMISSION.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256, (
        f"frozen submission changed\n  expected {FROZEN_SHA256}\n  found    {digest}")


@needs_data
def test_frozen_submission_still_scores_what_was_recorded():
    from immrep.metric import official_score

    frame = pd.read_csv(SUBMISSION)
    series = frame.set_index(frame.columns[0])[frame.columns[1]]
    labels = _labels()
    score = official_score(labels, labels["tcr_id"].map(series).to_numpy())
    assert score == pytest.approx(FROZEN_PUBLIC, abs=5e-5)


@needs_data
def test_kaggle_truncates_rather_than_rounds():
    """The leaderboard showed 0.849 for a 0.8498 file. Guard the convention the
    audit depends on, since rounding would have shown 0.850."""
    assert f"{int(FROZEN_PUBLIC * 1000) / 1000:.3f}" == "0.849"


# --------------------------------------------------------------------------
# The metric


@needs_data
def test_metric_reproduces_the_organisers_baseline():
    """E1: this is how the metric was identified in the first place."""
    from immrep.data import load_tcrdist_baseline
    from immrep.metric import official_score

    baseline = load_tcrdist_baseline()
    column = "prediction" if "prediction" in baseline else baseline.columns[1]
    series = baseline.set_index("tcr_id")[column]
    labels = _labels()
    score = official_score(labels, labels["tcr_id"].map(series).to_numpy())
    assert score == pytest.approx(K.BASELINES["tcrdist3"], abs=1.5e-3)


def test_metric_is_a_mean_of_per_epitope_aucs_not_a_pooled_auc():
    """The distinction E1 turned on: equal weight per epitope regardless of size.
    Built so a pooled AUC would give a visibly different answer."""
    from immrep.metric import official_score, per_epitope_auc

    labels = pd.DataFrame({
        "epitope_aa": ["A"] * 100 + ["B"] * 4,
        "pos": [1, 0] * 50 + [1, 1, 0, 0],
    })
    # Perfect ordering on the big epitope, inverted on the tiny one.
    predictions = np.r_[np.tile([1.0, 0.0], 50), [0.0, 0.0, 1.0, 1.0]]
    per = per_epitope_auc(labels, predictions)
    assert per["A"] == pytest.approx(1.0)
    assert per["B"] == pytest.approx(0.0)
    assert official_score(labels, predictions) == pytest.approx(0.5)


def test_single_class_epitope_scores_nan_rather_than_raising():
    """LOEO folds built from tiny subsets hit this."""
    from immrep.metric import per_epitope_auc

    labels = pd.DataFrame({"epitope_aa": ["A", "A"], "pos": [1, 1]})
    assert np.isnan(per_epitope_auc(labels, np.array([0.2, 0.8]))["A"])


# --------------------------------------------------------------------------
# Data invariants the whole project rests on


@needs_data
def test_test_set_shape_and_split():
    from immrep.data import load_test

    test = load_test()
    assert len(test) == K.N_TEST_ROWS
    assert test["tcr_id"].is_unique
    sizes = test.groupby("split").size().to_dict()
    assert sizes["public"] == K.N_PUBLIC
    assert sizes["private"] == K.N_PRIVATE


@needs_data
def test_private_rows_carry_no_label():
    """If this ever fails, the private leaderboard has leaked."""
    from immrep.data import load_test

    test = load_test()
    private = set(test.loc[test["split"] == "private", "tcr_id"])
    assert not (set(_labels()["tcr_id"]) & private)


@needs_data
def test_scored_rows_are_the_public_split_minus_ambiguous():
    from immrep.data import load_labels

    ambiguous = len(load_labels(drop_ambiguous=False)) - len(
        load_labels(drop_ambiguous=True))
    assert len(_labels()) + ambiguous == K.N_PUBLIC


def test_epitope_counts_sum_to_the_test_set():
    assert K.N_PUBLIC + K.N_PRIVATE == K.N_TEST_ROWS
    assert len(K.PUBLIC_EPITOPES) == 5 and len(K.PRIVATE_EPITOPES) == 3
    assert not set(K.PUBLIC_EPITOPES) & set(K.PRIVATE_EPITOPES)


# --------------------------------------------------------------------------
# The submission guard itself


@needs_data
@pytest.mark.parametrize("break_it", ["row_count", "nan", "constant", "columns"])
def test_validate_rejects_broken_submissions(break_it):
    """submission.validate is the last line of defence before a submission;
    check it actually catches each failure it claims to."""
    from immrep.data import load_test
    from immrep.submission import validate

    test = load_test()
    frame = pd.read_csv(SUBMISSION)
    if break_it == "row_count":
        frame = frame.iloc[:-1]
    elif break_it == "nan":
        frame.iloc[0, 1] = np.nan
    elif break_it == "constant":
        frame.iloc[:, 1] = 1.0
    elif break_it == "columns":
        frame = frame.rename(columns={frame.columns[1]: "wrong_name"})

    with pytest.raises(Exception):
        validate(frame, test=test)


@needs_data
def test_validate_accepts_the_frozen_submission():
    from immrep.data import load_test
    from immrep.submission import validate

    validate(pd.read_csv(SUBMISSION), test=load_test())
