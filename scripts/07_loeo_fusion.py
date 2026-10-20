#!/usr/bin/env python
"""Does fusing density with cluster shape beat the best single density feature?

Everything is scored by leave-one-epitope-out, which is the only split that
matches what the private leaderboard asks. The bar to clear is the best
unsupervised feature at 0.8396; a fitted model that lands below that is adding
variance, not information.

Run: python scripts/07_loeo_fusion.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame  # noqa: E402
from immrep.features import build_feature_matrix, load_space  # noqa: E402
from immrep.validation import (  # noqa: E402
    leave_one_epitope_out,
    single_feature_fitter,
    sklearn_fitter,
)

SPACE_DIR = ROOT / "data/interim/sceptr"
# `default` leads, `a_sceptr` is kept because the alpha chain carries more of
# this signal than beta and is the only space that helps GLCTLVAML.
SPACE_NAMES = ("default", "a_sceptr", "small")


def main() -> int:
    frame = build_training_frame()
    spaces = [load_space(SPACE_DIR / f"{n}.npz", n) for n in SPACE_NAMES]
    features = build_feature_matrix(frame, spaces)

    if features["tcr_id"].tolist() != frame["tcr_id"].tolist():
        raise ValueError("feature matrix lost row alignment with the training frame")

    labels = frame["pos"].to_numpy()
    epitopes = frame["epitope_aa"]
    value_columns = [c for c in features.columns if c not in ("tcr_id", "epitope_aa")]

    results: dict[str, pd.DataFrame] = {}

    # Unsupervised references, no labels involved.
    for column in ("default__kern0.1", "a_sceptr__ntopk15", "clu_entropy_b"):
        results[f"[single] {column}"] = leave_one_epitope_out(
            features, labels, epitopes, single_feature_fitter(column)
        )

    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    density_only = [c for c in value_columns if "__" in c]
    cluster_only = [c for c in value_columns if c.startswith(("clu_", "in_cluster"))]

    # Heavy regularisation on purpose: five folds cannot support model selection,
    # so the fitted models are kept simple enough that they mostly cannot overfit.
    for name, columns in (
        ("density only", density_only),
        ("cluster shape only", cluster_only),
        ("density + cluster", value_columns),
    ):
        results[f"[logreg C=0.1] {name}"] = leave_one_epitope_out(
            features,
            labels,
            epitopes,
            sklearn_fitter(
                lambda: make_pipeline(
                    StandardScaler(), LogisticRegression(C=0.1, max_iter=2000)
                ),
                columns,
            ),
        )

    results["[hgb depth=2] density + cluster"] = leave_one_epitope_out(
        features,
        labels,
        epitopes,
        sklearn_fitter(
            lambda: HistGradientBoostingClassifier(
                max_depth=2, max_iter=150, learning_rate=0.05, l2_regularization=1.0
            ),
            value_columns,
        ),
    )

    table = pd.DataFrame(
        {name: report["auc"] for name, report in results.items()}
    ).T.sort_values("mean", ascending=False)
    print("bar to clear: 0.8396 (best unsupervised single feature, full-data estimate)\n")
    print(table.round(4).to_string())

    out = ROOT / "reports/03_loeo_fusion.md"
    out.write_text(
        "# LOEO: density vs cluster shape vs fusion\n\n"
        "Mean per-epitope AUC, leave-one-epitope-out. Columns are the held-out "
        "epitope.\n\n" + table.round(4).to_markdown() + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
