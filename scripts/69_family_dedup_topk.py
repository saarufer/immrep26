#!/usr/bin/env python
"""E77 -- de-duplicating the top-k neighbourhood by family.

Externally proposed: a TCR's eight nearest neighbours may represent only two or
three independent recognition modes, so capping each family's contribution
should sharpen the ranking. Uses family rather than donor, so unlike E57's
machinery it is computable on the private split.

Run: .venv/bin/python scripts/69_family_dedup_topk.py
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/44_family_dedup.md"
RADIUS, K = 0.12, 8


def main() -> int:
    test, labels = load_test(), build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    standard, deduped, effective = {}, {}, {}
    for _, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, family = connected_components(
            sp.csr_matrix(d <= RADIUS), directed=False)
        np.fill_diagonal(d, np.inf)
        for i, tcr in enumerate(ids):
            order = np.argsort(d[i])
            standard[tcr] = -d[i][order[:K]].mean()
            # One representative per family, taken nearest-first. The query's
            # own family is not excluded: keeping one member of it is the
            # point, since the proposal is to cap a family rather than drop it.
            seen, picked = set(), []
            for j in order:
                if family[j] in seen:
                    continue
                seen.add(family[j])
                picked.append(d[i, j])
                if len(picked) == K:
                    break
            deduped[tcr] = -float(np.mean(picked))
            effective[tcr] = len(set(family[order[:K]]))

    rows = []
    for epitope, cohort in test.groupby("epitope_aa"):
        values = np.array([effective[t] for t in cohort["tcr_id"]])
        rows.append({"epitope": epitope, "split": cohort["split"].iloc[0],
                     "effective_families": float(values.mean()),
                     "all_one_family": float(np.mean(values == 1))})
        print(f"  {epitope:<11} {cohort['split'].iloc[0]:<8} effective families "
              f"{values.mean():.2f}/{K}   all one family "
              f"{np.mean(values == 1):.0%}")
    shape = pd.DataFrame(rows)

    def macro(series):
        return float(np.mean([roc_auc_score(b["pos"], b["tcr_id"].map(series))
                              for _, b in labels.groupby("epitope_aa")]))
    std, ded = pd.Series(standard), pd.Series(deduped)
    eff = pd.Series(effective)
    rho = float(np.mean([spearmanr(b["tcr_id"].map(std), b["tcr_id"].map(ded)
                                   ).statistic
                         for _, b in test.groupby("epitope_aa")]))
    print(f"\n  standard {macro(std):.4f}   de-duplicated {macro(ded):.4f}   "
          f"delta {macro(ded) - macro(std):+.4f}")
    print(f"  effective family count scored alone: {macro(eff):.4f}")
    print(f"  within-epitope Spearman between the two: {rho:+.4f}")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E77 - family-deduplicated top-k\n\n" + shape.round(4).to_markdown(index=False)
        + f"\n\n- standard **{macro(std):.4f}**, de-duplicated **{macro(ded):.4f}**\n"
        f"- effective family count alone: **{macro(eff):.4f}**\n"
        f"- within-epitope Spearman: {rho:+.4f}\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
