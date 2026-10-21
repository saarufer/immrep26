#!/usr/bin/env python
"""Six remaining sequence-side cards, all sharing the density machinery.

D04 chemical substitution cost     -- is edit distance the wrong metric for CDR3?
D05 multi-scale neighbourhood      -- is one bandwidth throwing away structure?
D10 applicability map              -- where does the model have no evidence?
S01 residual learnable information -- can 1,782 labels add anything at all?
S02 singleton rescue model         -- a model fitted only where the error is
S07 abstention layer               -- say which predictions are unsupported

Each is judged the same way: strict nested leave-one-epitope-out against the
submitted blend, with the pre-registered +0.02 inclusion threshold. S07 is the
exception -- it is not scored by AUC at all, since the competition metric has no
way to reward declining to answer.

Run: python scripts/27_remaining_sequence_cards.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import parasail
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.metric import per_epitope_auc  # noqa: E402

SIGMA = 0.1
STRUCTURE_WEIGHT = 0.3
CLUSTER_RADIUS = 0.12


def rank_within(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    return frame.assign(r=frame.groupby("epitope_aa")["v"].rank(pct=True)).set_index("tcr_id")["r"]


def sceptr_distances(space: str, test: pd.DataFrame) -> dict[str, tuple[np.ndarray, list[str]]]:
    cached = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].tolist()
        block = vectors[[position[t] for t in ids]]
        out[epitope] = (1.0 - block @ block.T, ids)
    return out


def density_from(matrices: dict, test: pd.DataFrame, sigma: float) -> pd.Series:
    out = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, (distances, ids) in matrices.items():
        masked = distances.copy()
        np.fill_diagonal(masked, np.inf)
        out.loc[ids] = np.exp(-masked / sigma).mean(axis=1)
    return out


def cdr3_distance(cohort: pd.DataFrame, matrix, gap_open: int, gap_extend: int) -> np.ndarray:
    """Alignment distance on the two CDR3 loops only, under a given cost matrix."""
    n = len(cohort)
    total = np.zeros((n, n), dtype=np.float32)
    for column in ("cdr3_alpha_aa", "cdr3_beta_aa"):
        seqs = cohort[column].tolist()
        self_score = np.array(
            [parasail.nw_scan_16(s, s, gap_open, gap_extend, matrix).score for s in seqs],
            dtype=float,
        )
        block = np.zeros((n, n), dtype=np.float32)
        for i in range(n):
            for j in range(i + 1, n):
                score = parasail.nw_scan_16(seqs[i], seqs[j], gap_open, gap_extend, matrix).score
                denominator = np.sqrt(self_score[i] * self_score[j])
                value = 1.0 - score / denominator if denominator > 0 else 1.0
                block[i, j] = block[j, i] = value
        total += block
    return total


def main() -> int:
    test = load_test()
    labels = build_training_frame()
    epitopes = sorted(labels["epitope_aa"].unique())
    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    contacts = pd.read_csv(ROOT / "data/interim/structure_contacts.csv").set_index("tcr_id")

    default = sceptr_distances("default", test)
    small = sceptr_distances("small", test)
    tcrdist = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz", allow_pickle=True)
        pos = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [pos[t] for t in cohort["tcr_id"]]
        tcrdist[epitope] = (
            (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)],
            cohort["tcr_id"].tolist(),
        )

    blend = (
        rank_within(density_from(default, test, SIGMA), test)
        + rank_within(density_from(small, test, SIGMA), test)
        + rank_within(density_from(tcrdist, test, 20.0), test)
    ) / 3
    submitted = (1 - STRUCTURE_WEIGHT) * blend + STRUCTURE_WEIGHT * rank_within(
        structure["pae_tcr_pep_min"], test
    )
    base_per = {e: roc_auc_score(
        labels[labels.epitope_aa == e].pos,
        labels[labels.epitope_aa == e].tcr_id.map(submitted)) for e in epitopes}
    base = float(np.mean(list(base_per.values())))
    pd.set_option("display.width", 250)
    print(f"submitted blend = {base:.4f}\n")

    # ---------------- D04 ----------------
    print("=== D04: chemical substitution cost on the CDR3 loops ===")
    print("   same alignment, three pre-declared cost matrices, nothing fitted\n")
    schemes = {
        "identity (edit-like)": (parasail.matrix_create("ACDEFGHIKLMNPQRSTVWY", 1, 0), 2, 1),
        "BLOSUM62": (parasail.blosum62, 11, 1),
        "BLOSUM45 (permissive)": (parasail.blosum45, 11, 1),
    }
    rows = []
    for name, (matrix, go, ge) in schemes.items():
        values = pd.Series(index=test["tcr_id"], dtype=float)
        for epitope, cohort in test.groupby("epitope_aa"):
            cohort = cohort.reset_index(drop=True)
            distances = cdr3_distance(cohort, matrix, go, ge)
            np.fill_diagonal(distances, np.inf)
            values.loc[cohort["tcr_id"]] = np.exp(-distances / 0.12).mean(axis=1)
        aucs = per_epitope_auc(labels, labels["tcr_id"].map(values).to_numpy())
        rows.append({"scheme": name, "mean": aucs.mean(), **aucs.round(3).to_dict()})
    print(pd.DataFrame(rows).round(4).to_string(index=False))

    # ---------------- D05 ----------------
    print("\n=== D05: multi-scale neighbourhood ===")
    scales = (0.03, 0.1, 0.3)
    parts = [rank_within(density_from(default, test, s), test) for s in scales]
    rows = [{"variant": f"single sigma={s}",
             "mean": per_epitope_auc(labels, labels["tcr_id"].map(p).to_numpy()).mean()}
            for s, p in zip(scales, parts)]
    multi = sum(parts) / len(parts)
    rows.append({"variant": "rank-average of all three",
                 "mean": per_epitope_auc(labels, labels["tcr_id"].map(multi).to_numpy()).mean()})
    print(pd.DataFrame(rows).round(4).to_string(index=False))
    for weight in (0.2, 0.3):
        combined = (1 - weight) * submitted + weight * multi
        score = per_epitope_auc(labels, labels["tcr_id"].map(combined).to_numpy()).mean()
        print(f"  blended into the submission at w={weight}: {score:.4f}  gain {score - base:+.4f}")

    # ---------------- D10 ----------------
    print("\n=== D10: applicability map -- where is there no evidence? ===")
    strata = []
    for epitope, (distances, ids) in default.items():
        masked = distances.copy()
        np.fill_diagonal(masked, np.inf)
        nearest = masked.min(axis=1)
        strata.append(pd.DataFrame({"tcr_id": ids, "epitope_aa": epitope, "nn1": nearest}))
    strata = pd.concat(strata, ignore_index=True)
    strata["band"] = pd.cut(strata["nn1"], [0, 0.06, 0.12, 0.25, 0.5, 2.0],
                            labels=["<0.06", "0.06-0.12", "0.12-0.25", "0.25-0.5", ">0.5"])
    joined = labels.merge(strata[["tcr_id", "nn1", "band"]], on="tcr_id")
    joined["score"] = joined["tcr_id"].map(submitted)
    table = joined.groupby("band", observed=True).apply(
        lambda g: pd.Series({
            "n": len(g), "pos_rate": g.pos.mean(),
            "auc": roc_auc_score(g.pos, g.score) if g.pos.nunique() == 2 else np.nan,
        }), include_groups=False)
    print(table.round(3).to_string())
    private = strata[strata["epitope_aa"].isin(("CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"))]
    print("\n  private-side occupancy of those bands:")
    print(pd.crosstab(private["epitope_aa"], private["band"], normalize="index").round(3).to_string())

    # ---------------- S01 / S02 ----------------
    print("\n=== S01: can the 1,782 labels add anything on top of the frozen score? ===")
    print("   one extra feature, L2-heavy logistic regression, nested LOEO\n")
    extras = {
        "contacts_a_pep": rank_within(contacts["contacts_a_pep"], test),
        "iptm": rank_within(structure["iptm"], test),
        "multi-scale density": multi,
    }
    for name, extra in extras.items():
        folds = []
        for held_out in epitopes:
            train = labels[labels.epitope_aa != held_out]
            test_fold = labels[labels.epitope_aa == held_out]
            X_train = np.column_stack([train.tcr_id.map(submitted), train.tcr_id.map(extra)])
            X_test = np.column_stack([test_fold.tcr_id.map(submitted), test_fold.tcr_id.map(extra)])
            model = LogisticRegression(C=0.05, max_iter=2000).fit(X_train, train.pos)
            folds.append(roc_auc_score(test_fold.pos, model.predict_proba(X_test)[:, 1]))
        print(f"  + {name:<22} nested {np.mean(folds):.4f}   gain {np.mean(folds) - base:+.4f}")

    print("\n=== S02: a model fitted only on singletons, spliced into the global order ===")
    folds, folds_singleton = [], []
    for held_out in epitopes:
        pieces = []
        for epitope in epitopes:
            distances, ids = default[epitope]
            _, family = connected_components(
                sp.csr_matrix(distances <= CLUSTER_RADIUS), directed=False)
            sizes = pd.Series(family).value_counts()
            pieces.append(pd.DataFrame({
                "tcr_id": ids,
                "singleton": ~np.isin(family, sizes[sizes >= 2].index)}))
        flags = pd.concat(pieces, ignore_index=True).set_index("tcr_id")["singleton"]
        frame = labels.assign(
            singleton=labels.tcr_id.map(flags),
            base=labels.tcr_id.map(submitted),
            nn1=labels.tcr_id.map(strata.set_index("tcr_id")["nn1"]),
            contacts=labels.tcr_id.map(contacts["contacts_a_pep"]))
        train = frame[(frame.epitope_aa != held_out) & frame.singleton]
        target = frame[frame.epitope_aa == held_out]
        model = LogisticRegression(C=0.1, max_iter=2000).fit(
            train[["nn1", "contacts"]], train.pos)
        adjusted = target.base.to_numpy().copy()
        mask = target.singleton.to_numpy()
        if mask.any():
            # Re-rank only the singletons, inside the score range they occupy.
            fitted = model.predict_proba(target.loc[mask, ["nn1", "contacts"]])[:, 1]
            lo, hi = adjusted[mask].min(), adjusted[mask].max()
            adjusted[mask] = lo + (hi - lo) * pd.Series(fitted).rank(pct=True).to_numpy()
        folds.append(roc_auc_score(target.pos, adjusted))
        if mask.any() and target.pos[mask].nunique() == 2:
            folds_singleton.append(roc_auc_score(target.pos[mask], adjusted[mask]))
    print(f"  global macro after singleton re-ranking: {np.mean(folds):.4f}   "
          f"gain {np.mean(folds) - base:+.4f}")
    print(f"  within-singleton AUC after re-ranking:   {np.mean(folds_singleton):.4f}")

    # ---------------- S07 ----------------
    print("\n=== S07: abstention layer (not scored by the competition metric) ===")
    joined["band"] = joined["band"].astype(str)
    joined["confident"] = joined["nn1"] <= 0.12
    for name, mask in (("has a close neighbour", joined.confident),
                       ("isolated", ~joined.confident)):
        subset = joined[mask]
        auc = roc_auc_score(subset.pos, subset.score) if subset.pos.nunique() == 2 else np.nan
        print(f"  {name:<22} n={len(subset):5d}  pos_rate={subset.pos.mean():.3f}  AUC={auc:.3f}")
    print("\n  The metric cannot reward abstention -- every row is scored -- so this is")
    print("  a reporting layer for downstream use, not a submission component.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
