#!/usr/bin/env python
"""Confirm the official metric by scoring the organisers' tcrdist3 baseline
(public LB = 0.824) against the published labels under every candidate metric."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immrep.data import load_labels, load_tcrdist_baseline  # noqa: E402
from immrep.metric import TCRDIST3_PUBLIC_SCORE, official_score, score_report  # noqa: E402


def main() -> int:
    raw = load_labels(drop_ambiguous=False).merge(
        load_tcrdist_baseline()[["tcr_id", "prediction"]], on="tcr_id"
    )
    scored = raw[raw.ambiguous == 0].copy()
    scored["pos"] = scored["pos"].astype(int)

    def mean_ep(df, y):
        return (
            df.assign(_y=y)
            .groupby("epitope_aa")
            .apply(lambda g: roc_auc_score(g._y, g.prediction), include_groups=False)
            .mean()
        )

    cands = {
        "mean per-epitope ROC AUC (ambiguous dropped)": official_score(scored, scored.prediction),
        "pooled ROC AUC": roc_auc_score(scored.pos, scored.prediction),
        "pooled average precision": average_precision_score(scored.pos, scored.prediction),
        "mean per-epitope average precision": scored.groupby("epitope_aa")
        .apply(lambda g: average_precision_score(g.pos, g.prediction), include_groups=False)
        .mean(),
        "n-weighted mean per-epitope ROC AUC": (
            lambda a, n: (a * n).sum() / n.sum()
        )(
            scored.groupby("epitope_aa").apply(
                lambda g: roc_auc_score(g.pos, g.prediction), include_groups=False
            ),
            scored.groupby("epitope_aa").size(),
        ),
        "mean per-epitope AUC, ambiguous as negative": mean_ep(raw, (raw.pos == 1).astype(int)),
        "mean per-epitope AUC, ambiguous as positive": mean_ep(
            raw, ((raw.pos == 1) | (raw.ambiguous == 1)).astype(int)
        ),
    }
    print(f"target public-leaderboard score: {TCRDIST3_PUBLIC_SCORE}\n")
    for name, value in sorted(cands.items(), key=lambda kv: abs(kv[1] - TCRDIST3_PUBLIC_SCORE)):
        hit = "  <== MATCH" if abs(value - TCRDIST3_PUBLIC_SCORE) < 0.0015 else ""
        print(f"  {value:.6f}  {name}{hit}")
    print("\nper-epitope detail for the baseline:")
    print(score_report(scored, scored.prediction).round(4).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
