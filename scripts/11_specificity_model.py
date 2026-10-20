#!/usr/bin/env python
"""Epitope-conditioned specificity: does this TCR look like a binder of *this* peptide?

Every signal tried so far counts neighbours inside the cohort and ignores the
peptide entirely, and all of them saturate near 0.84
(reports/02_signal_landscape.md). This model reads the peptide, so it is
neighbour-free and can in principle rank the 34-83% of each cohort that sits in
no cluster.

The design constraint that makes or breaks it: **the target epitope is removed
from training entirely.** The test TCRs are themselves VDJdb entries annotated
to that epitope, so a model that has seen them would score every one of them
positive and the feature would be constant. Holding the epitope out turns this
into genuine unseen-epitope specificity prediction, which the literature puts
at AUC 0.55-0.65 -- so the honest expectation here is low.

Negatives are the same TCRs paired with *wrong* epitopes. That is what forces
the model to learn compatibility: TCR identity alone carries no label
information when each TCR appears as both a positive and a negative.

Run: python scripts/11_specificity_model.py
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

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.metric import per_epitope_auc  # noqa: E402

VDJDB = ROOT / "data/external/vdjdb/vdjdb-2026-06-03/vdjdb_full.txt"
OUT = ROOT / "data/interim/specificity_scores.csv"

# The two largest epitopes hold 44k of 76k paired TCRs; without a cap the model
# would mostly learn those two peptides.
MAX_PER_EPITOPE = 1500
MIN_PER_EPITOPE = 15
NEGATIVES_PER_POSITIVE = 4
MAX_PEPTIDE_LEN = 12
SEED = 0

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def blosum62() -> dict[str, np.ndarray]:
    """Per-residue BLOSUM62 rows, so an unseen peptide still lands somewhere
    sensible -- one-hot positions cannot generalise across substitutions."""
    from Bio.Align import substitution_matrices

    matrix = substitution_matrices.load("BLOSUM62")
    return {
        aa: np.array([float(matrix[aa][other]) for other in AMINO_ACIDS], dtype=np.float32)
        for aa in AMINO_ACIDS
    }


def encode_peptides(peptides: pd.Series, table: dict[str, np.ndarray]) -> np.ndarray:
    """Left-aligned BLOSUM encoding padded to MAX_PEPTIDE_LEN, plus length.

    Anchor residues sit at both ends of an MHC-I epitope, so the last three
    positions are also encoded right-aligned to keep the C-terminal anchor at a
    fixed offset regardless of peptide length.
    """
    width = MAX_PEPTIDE_LEN * len(AMINO_ACIDS)
    out = np.zeros((len(peptides), width + 3 * len(AMINO_ACIDS) + 1), dtype=np.float32)
    for row, peptide in enumerate(peptides):
        for i, aa in enumerate(peptide[:MAX_PEPTIDE_LEN]):
            if aa in table:
                out[row, i * len(AMINO_ACIDS) : (i + 1) * len(AMINO_ACIDS)] = table[aa]
        for k, aa in enumerate(peptide[-3:]):
            if aa in table:
                start = width + k * len(AMINO_ACIDS)
                out[row, start : start + len(AMINO_ACIDS)] = table[aa]
        out[row, -1] = len(peptide)
    return out


def load_vdjdb_positives() -> pd.DataFrame:
    frame = pd.read_csv(VDJDB, sep="\t", low_memory=False)
    frame = frame[
        (frame["species"] == "HomoSapiens")
        & frame["cdr3.alpha"].notna()
        & frame["cdr3.beta"].notna()
        & (frame["mhc.class"] == "MHCI")
    ].copy()
    frame = frame.rename(
        columns={
            "cdr3.alpha": "cdr3a",
            "cdr3.beta": "cdr3b",
            "v.alpha": "trav",
            "j.alpha": "traj",
            "v.beta": "trbv",
            "j.beta": "trbj",
            "antigen.epitope": "epitope",
            "mhc.a": "hla",
        }
    )
    frame = frame[
        frame["epitope"].str.len().between(8, MAX_PEPTIDE_LEN)
        & frame["epitope"].str.fullmatch(f"[{AMINO_ACIDS}]+")
    ]
    frame = frame.drop_duplicates(["cdr3a", "cdr3b", "epitope"])

    counts = frame["epitope"].value_counts()
    keep = counts[counts >= MIN_PER_EPITOPE].index
    frame = frame[frame["epitope"].isin(keep)]
    # Cap, deterministically.
    frame = (
        frame.groupby("epitope", group_keys=False)
        .apply(lambda g: g.sample(min(len(g), MAX_PER_EPITOPE), random_state=SEED))
        .reset_index(drop=True)
    )
    return frame[["cdr3a", "cdr3b", "trav", "traj", "trbv", "trbj", "epitope", "hla"]]


def sceptr_vectors(frame: pd.DataFrame) -> np.ndarray:
    """SCEPTR embeddings for arbitrary (cdr3, V, J) rows."""
    import sceptr
    from libtcrlm.schema.tcr import TravGene, TrbvGene

    model_input = pd.DataFrame(
        {
            "TRAV": frame["trav"],
            "CDR3A": frame["cdr3a"],
            "TRAJ": frame["traj"],
            "TRBV": frame["trbv"],
            "CDR3B": frame["cdr3b"],
            "TRBJ": frame["trbj"],
        }
    )
    for column, vocab in (("TRAV", {g.name for g in TravGene}), ("TRBV", {g.name for g in TrbvGene})):
        base = model_input[column].str.split("*").str[0]
        model_input.loc[~base.isin(vocab), column] = None
    return sceptr.calc_vector_representations(model_input)


def build_design(
    positives: pd.DataFrame,
    tcr_vectors: np.ndarray,
    peptide_table: dict[str, np.ndarray],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Positives plus shuffled-epitope negatives, as one design matrix."""
    epitope_pool = positives["epitope"].unique()
    # sqrt-frequency weighting: uniform sampling would draw almost only the
    # long tail, frequency sampling almost only the two giants.
    weights = np.sqrt(positives["epitope"].value_counts().reindex(epitope_pool).to_numpy())
    weights = weights / weights.sum()

    tcr_rows = np.arange(len(positives))
    peptides = [positives["epitope"].to_numpy()]
    rows = [tcr_rows]
    labels = [np.ones(len(positives), dtype=np.int8)]

    for _ in range(NEGATIVES_PER_POSITIVE):
        drawn = rng.choice(epitope_pool, size=len(positives), p=weights)
        # Resample the rare self-collisions rather than mislabel them.
        collision = drawn == positives["epitope"].to_numpy()
        while collision.any():
            drawn[collision] = rng.choice(epitope_pool, size=int(collision.sum()), p=weights)
            collision = drawn == positives["epitope"].to_numpy()
        peptides.append(drawn)
        rows.append(tcr_rows)
        labels.append(np.zeros(len(positives), dtype=np.int8))

    peptide_encoding = encode_peptides(pd.Series(np.concatenate(peptides)), peptide_table)
    design = np.hstack([tcr_vectors[np.concatenate(rows)], peptide_encoding])
    return design, np.concatenate(labels)


def main() -> int:
    rng = np.random.default_rng(SEED)
    peptide_table = blosum62()

    positives = load_vdjdb_positives()
    print(f"VDJdb training pool: {len(positives)} TCR-epitope pairs, "
          f"{positives['epitope'].nunique()} epitopes")

    print("embedding training TCRs ...")
    train_vectors = sceptr_vectors(positives)

    test = load_test()
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    test_vectors = cached["vectors"][[position[t] for t in test["tcr_id"]]]
    test_peptides = encode_peptides(test["epitope_aa"], peptide_table)

    from lightgbm import LGBMClassifier

    scores = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope in sorted(test["epitope_aa"].unique()):
        mask = (positives["epitope"] != epitope).to_numpy()
        design, labels = build_design(
            positives.loc[mask].reset_index(drop=True), train_vectors[mask], peptide_table, rng
        )
        model = LGBMClassifier(
            n_estimators=400,
            learning_rate=0.05,
            num_leaves=63,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            verbose=-1,
            random_state=SEED,
        )
        model.fit(design, labels)

        rows = (test["epitope_aa"] == epitope).to_numpy()
        held_out = np.hstack([test_vectors[rows], test_peptides[rows]])
        scores.loc[test.loc[rows, "tcr_id"]] = model.predict_proba(held_out)[:, 1]
        print(f"  {epitope:<11} trained on {int(mask.sum())} positives, scored {int(rows.sum())}")

    scores.rename("specificity").to_csv(OUT)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

    labelled = build_training_frame()
    aucs = per_epitope_auc(labelled, labelled["tcr_id"].map(scores).to_numpy())
    print("\nper-epitope AUC of the held-out specificity score:")
    print(aucs.round(4).to_string())
    print(f"  mean = {aucs.mean():.4f}   (density reference 0.8396)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
