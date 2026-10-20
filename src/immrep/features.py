"""Feature construction for validation-reliability ranking.

Every feature is computed **within an epitope cohort and without labels**, so
the same code runs unchanged on the three private epitopes. Two families:

*Density* -- how crowded is this TCR's neighbourhood. Saturates around 0.84
mean per-epitope AUC on its own, in any representation (see
reports/02_signal_landscape.md).

*Cluster shape* -- what the neighbourhood looks like once you are in one. These
exist because density conflates two things it should not: a convergent,
antigen-driven motif family and a tight knot of near-duplicate sequences that
entered VDJdb together. On the public epitopes, tighter clusters validate
*less* (AUC 0.275 for internal tightness) while clusters with more
per-position CDR3-beta diversity validate *more* (0.773), which is the
convergent-recombination signature separating the two.

Features are rank-normalised within epitope before fusion. The metric only
reads within-epitope ranking, and raw scales differ enough between cohorts that
an unnormalised feature would make a fitted model cohort-specific.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

# Cosine radius defining a cluster in SCEPTR space. 0.12 puts roughly a quarter
# of each cohort in a cluster of two or more.
CLUSTER_RADIUS = 0.12
KERNEL_SIGMAS = (0.05, 0.1, 0.2)
TOPK = (1, 3, 8, 15)


@dataclass(frozen=True)
class Space:
    """A named representation plus its per-TCR vectors, indexed by tcr_id."""

    name: str
    vectors: np.ndarray
    index: dict[str, int]

    def cohort(self, tcr_ids: list[str]) -> np.ndarray:
        rows = [self.index[t] for t in tcr_ids]
        block = self.vectors[rows]
        return block / np.linalg.norm(block, axis=1, keepdims=True)


def load_space(path, name: str) -> Space:
    cached = np.load(path, allow_pickle=True)
    return Space(
        name=name,
        vectors=cached["vectors"],
        index={t: i for i, t in enumerate(cached["tcr_id"])},
    )


def _mean_position_entropy(sequences: pd.Series) -> float:
    """Mean per-position amino-acid entropy over a set of CDR3s.

    High entropy means many different sequences share this neighbourhood, i.e.
    independent recombination events converging on a motif. Low entropy means
    near-duplicates.
    """
    if len(sequences) < 2:
        return np.nan
    length = int(np.median(sequences.str.len()))
    entropies = []
    for i in range(length):
        column = [s[i] for s in sequences if len(s) > i]
        if not column:
            continue
        p = pd.Series(column).value_counts(normalize=True)
        entropies.append(float(-(p * np.log(p)).sum()))
    return float(np.mean(entropies)) if entropies else np.nan


def density_features(distances: np.ndarray, space_name: str) -> dict[str, np.ndarray]:
    """Kernel and k-nearest-neighbour reductions of one cohort's distances."""
    off_diagonal = distances.copy()
    np.fill_diagonal(off_diagonal, np.inf)
    sorted_distances = np.sort(off_diagonal, axis=1)

    out: dict[str, np.ndarray] = {}
    for sigma in KERNEL_SIGMAS:
        out[f"{space_name}__kern{sigma}"] = np.exp(-off_diagonal / sigma).mean(axis=1)
    for k in TOPK:
        # Negated so that larger always means "denser neighbourhood".
        out[f"{space_name}__ntopk{k}"] = -sorted_distances[:, :k].mean(axis=1)
    return out


def cluster_features(
    distances: np.ndarray, cohort: pd.DataFrame, normed: np.ndarray
) -> dict[str, np.ndarray]:
    """Shape of the single-linkage cluster each TCR belongs to.

    Singletons get NaN for every shape feature, which is honest: a cluster of
    one has no shape. Downstream models receive an explicit `in_cluster` flag
    and the NaNs are imputed per fold.
    """
    adjacency = sp.csr_matrix(distances <= CLUSTER_RADIUS)
    _, labels = connected_components(adjacency, directed=False)
    sizes = pd.Series(labels).value_counts()
    biggest = np.where(labels == sizes.index[0])[0]

    n = len(cohort)
    names = [
        "clu_size",
        "clu_tightness",
        "clu_vb_purity",
        "clu_va_purity",
        "clu_entropy_b",
        "clu_entropy_a",
        "clu_centrality",
        "clu_near_biggest",
        "in_cluster",
    ]
    out = {name: np.full(n, np.nan) for name in names}

    cohort_centroid = normed.mean(axis=0)
    cohort_centroid /= np.linalg.norm(cohort_centroid)

    for label in np.unique(labels):
        members = np.where(labels == label)[0]
        out["clu_size"][members] = len(members)
        out["in_cluster"][members] = float(len(members) >= 2)
        if len(members) < 2:
            continue

        block = distances[np.ix_(members, members)]
        upper = np.triu_indices(len(members), 1)
        # Negated: larger means looser, which is the direction that validates.
        out["clu_tightness"][members] = -float(block[upper].mean())

        out["clu_vb_purity"][members] = (
            cohort["TRBV_IMGT"].iloc[members].value_counts(normalize=True).iloc[0]
        )
        out["clu_va_purity"][members] = (
            cohort["TRAV_IMGT"].iloc[members].value_counts(normalize=True).iloc[0]
        )
        out["clu_entropy_b"][members] = _mean_position_entropy(
            cohort["cdr3_beta_aa"].iloc[members]
        )
        out["clu_entropy_a"][members] = _mean_position_entropy(
            cohort["cdr3_alpha_aa"].iloc[members]
        )
        # How central this cluster sits in the cohort. The failure mode it
        # targets: a tight knot of co-annotated artifacts forms off to one side
        # of the real response, dense but peripheral.
        cluster_centroid = normed[members].mean(axis=0)
        out["clu_centrality"][members] = float(
            cluster_centroid @ cohort_centroid / np.linalg.norm(cluster_centroid)
        )
        out["clu_near_biggest"][members] = (
            0.0 if label == sizes.index[0] else -float(distances[np.ix_(members, biggest)].min())
        )
    return out


def build_feature_matrix(
    frame: pd.DataFrame, spaces: list[Space], *, rank_normalise: bool = True
) -> pd.DataFrame:
    """Per-TCR features for every row of `frame`, grouped by epitope.

    `frame` needs `tcr_id`, `epitope_aa`, the two CDR3 columns and the V-gene
    columns. Cluster shape is derived from the first space in `spaces`, which
    should be the best-performing representation.
    """
    blocks = []
    for epitope, cohort in frame.groupby("epitope_aa", sort=True):
        cohort = cohort.reset_index(drop=True)
        ids = cohort["tcr_id"].tolist()
        columns: dict[str, np.ndarray] = {}

        primary_distances = None
        primary_normed = None
        for space in spaces:
            normed = space.cohort(ids)
            distances = 1.0 - normed @ normed.T
            np.fill_diagonal(distances, 0.0)
            if primary_distances is None:
                primary_distances, primary_normed = distances, normed
            columns.update(density_features(distances, space.name))

        columns.update(cluster_features(primary_distances, cohort, primary_normed))

        block = pd.DataFrame(columns)
        block.insert(0, "tcr_id", ids)
        block.insert(1, "epitope_aa", epitope)
        blocks.append(block)

    features = pd.concat(blocks, ignore_index=True)

    if rank_normalise:
        value_columns = [c for c in features.columns if c not in ("tcr_id", "epitope_aa")]
        features[value_columns] = features.groupby("epitope_aa")[value_columns].rank(pct=True)

    # Preserve the caller's row order so downstream joins stay aligned.
    order = {t: i for i, t in enumerate(frame["tcr_id"])}
    return features.sort_values("tcr_id", key=lambda s: s.map(order)).reset_index(drop=True)
