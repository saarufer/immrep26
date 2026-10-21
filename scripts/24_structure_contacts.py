#!/usr/bin/env python
"""Interface contacts, docking geometry and physical QC from the folded complexes.

Covers plan cards D07, C01, C02 and C06 in one pass over the 2,636 Boltz
structures, because they all need the same atom coordinates.

**D07 -- contact fingerprints.** Every confidence read tried so far (E14) is a
scalar summarising a whole interface, and all of them turned out to track how
canonical the TCR sequence is rather than what it touches. A fingerprint indexed
by *peptide position* is different: it says which part of the antigen each TCR
engages. The peptide is identical for every complex within an epitope, so
peptide position is a common coordinate and no IMGT numbering is needed --
which also avoids depending on a numbering tool for the one thing that has to be
exactly right.

Per complex, for each peptide position: heavy-atom contact counts and minimum
distances from TCR alpha and TCR beta, and how much of that contact comes from
CDR3 rather than the germline loops. CDR3 is located by exact substring match in
the assembled domain, which succeeds for all 2,636 rows.

**C01 -- geometry QC.** Identity was checked in scripts/21; this checks physics:
inter-chain heavy-atom clashes, whether the peptide actually sits in the MHC
groove, and whether the TCR actually touches the pMHC. A confident-looking
structure that fails these should not contribute a contact feature.

**C06 -- docking geometry.** The alpha-to-beta centroid vector against the
peptide N-to-C vector approximates the crossing angle, plus how far the TCR
centroid sits above the groove.

Run: python scripts/24_structure_contacts.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

PREDICTIONS = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
OUT_FEATURES = ROOT / "data/interim/structure_contacts.csv"
OUT_FINGERPRINTS = ROOT / "data/interim/contact_fingerprints.npz"

CONTACT_CUTOFF = 4.5   # heavy-atom contact, standard for interface definition
CLASH_CUTOFF = 2.2     # closer than any real non-bonded heavy-atom pair
MAX_PEPTIDE_LEN = 10


def read_atoms(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Heavy-atom coordinates and residue indices per chain."""
    coords: dict[str, list] = {}
    resid: dict[str, list] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip() or line[12:16].strip()[:1]
        if element == "H":
            continue
        chain = line[21]
        coords.setdefault(chain, []).append(
            (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        )
        resid.setdefault(chain, []).append(int(line[22:26]))
    return {
        c: (np.asarray(coords[c], dtype=np.float32), np.asarray(resid[c], dtype=np.int32))
        for c in coords
    }


def contact_profile(
    source: tuple[np.ndarray, np.ndarray],
    peptide: tuple[np.ndarray, np.ndarray],
    cdr3_residues: set[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-peptide-residue contact count, minimum distance, and CDR3 share."""
    src_xyz, src_res = source
    pep_xyz, pep_res = peptide
    positions = np.unique(pep_res)
    counts = np.zeros(MAX_PEPTIDE_LEN, dtype=np.float32)
    minima = np.full(MAX_PEPTIDE_LEN, 99.0, dtype=np.float32)
    cdr3_share = np.zeros(MAX_PEPTIDE_LEN, dtype=np.float32)

    tree = cKDTree(src_xyz)
    for slot, position in enumerate(positions[:MAX_PEPTIDE_LEN]):
        atoms = pep_xyz[pep_res == position]
        if not len(atoms):
            continue
        distances, indices = tree.query(atoms, k=min(64, len(src_xyz)))
        distances = np.atleast_2d(distances)
        indices = np.atleast_2d(indices)
        close = distances <= CONTACT_CUTOFF
        counts[slot] = float(close.sum())
        minima[slot] = float(distances.min())
        if close.any():
            touching = src_res[indices[close]]
            cdr3_share[slot] = float(np.isin(touching, list(cdr3_residues)).mean())
    return counts, minima, cdr3_share


def main() -> int:
    test = load_test().set_index("tcr_id")
    rows = []
    fingerprints: dict[str, np.ndarray] = {}

    ids = sorted(p.name for p in PREDICTIONS.iterdir() if p.is_dir())
    for index, tcr_id in enumerate(ids, 1):
        pdb = PREDICTIONS / tcr_id / f"{tcr_id}_model_0.pdb"
        if not pdb.exists():
            continue
        atoms = read_atoms(pdb)
        if not {"A", "B", "M", "P"} <= set(atoms):
            continue
        row = test.loc[tcr_id]

        # CDR3 residue numbers: located by substring position in the domain.
        cdr3 = {}
        for chain, full, loop in (
            ("A", row.full_seq_alpha_aa, row.cdr3_alpha_aa),
            ("B", row.full_seq_beta_aa, row.cdr3_beta_aa),
        ):
            start = full.index(loop)
            first = int(atoms[chain][1].min())
            cdr3[chain] = set(range(first + start, first + start + len(loop)))

        record: dict[str, object] = {"tcr_id": tcr_id}
        profile = []
        for chain in ("A", "B"):
            counts, minima, share = contact_profile(atoms[chain], atoms["P"], cdr3[chain])
            tag = "a" if chain == "A" else "b"
            record[f"contacts_{tag}_pep"] = float(counts.sum())
            record[f"mindist_{tag}_pep"] = float(minima.min())
            record[f"pep_positions_touched_{tag}"] = float((counts > 0).sum())
            record[f"cdr3_share_{tag}"] = float(
                np.average(share, weights=counts) if counts.sum() else 0.0
            )
            profile.append(np.concatenate([counts, np.minimum(minima, 20.0), share]))
        fingerprints[tcr_id] = np.concatenate(profile).astype(np.float32)

        # C01 geometry QC.
        clashes = 0
        chain_list = ["A", "B", "M", "P"]
        for i, first in enumerate(chain_list):
            tree = cKDTree(atoms[first][0])
            for second in chain_list[i + 1:]:
                clashes += int(sum(len(x) for x in tree.query_ball_point(atoms[second][0], CLASH_CUTOFF)))
        record["interchain_clashes"] = float(clashes)

        mhc_tree = cKDTree(atoms["M"][0])
        record["peptide_mhc_contacts"] = float(
            sum(len(x) for x in mhc_tree.query_ball_point(atoms["P"][0], CONTACT_CUTOFF))
        )
        record["tcr_touches_pmhc"] = float(
            record["contacts_a_pep"] + record["contacts_b_pep"] > 0
        )

        # C06 docking geometry.
        alpha_centre = atoms["A"][0].mean(axis=0)
        beta_centre = atoms["B"][0].mean(axis=0)
        peptide_xyz, peptide_res = atoms["P"]
        order = np.argsort(peptide_res)
        n_term = peptide_xyz[order][: max(1, len(order) // 5)].mean(axis=0)
        c_term = peptide_xyz[order][-max(1, len(order) // 5):].mean(axis=0)
        ab_vector = beta_centre - alpha_centre
        pep_vector = c_term - n_term
        cosine = float(
            ab_vector @ pep_vector
            / max(np.linalg.norm(ab_vector) * np.linalg.norm(pep_vector), 1e-6)
        )
        record["crossing_angle_deg"] = float(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
        tcr_centre = np.vstack([atoms["A"][0], atoms["B"][0]]).mean(axis=0)
        record["tcr_height_over_peptide"] = float(
            np.linalg.norm(tcr_centre - peptide_xyz.mean(axis=0))
        )
        rows.append(record)

        if index % 400 == 0:
            print(f"  {index}/{len(ids)} ...")

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_FEATURES, index=False)
    np.savez_compressed(
        OUT_FINGERPRINTS,
        tcr_id=np.array(list(fingerprints), dtype=object),
        fingerprint=np.vstack([fingerprints[t] for t in fingerprints]),
    )
    print(f"\nwrote {OUT_FEATURES.relative_to(ROOT)}  {len(frame)} complexes")
    print(f"wrote {OUT_FINGERPRINTS.relative_to(ROOT)}  dim={len(next(iter(fingerprints.values())))}")
    print("\n=== C01 geometry QC ===")
    print(f"  complexes with any inter-chain clash < {CLASH_CUTOFF} A: "
          f"{int((frame.interchain_clashes > 0).sum())}/{len(frame)}")
    print(f"  median clashes: {frame.interchain_clashes.median():.0f}")
    print(f"  peptide-MHC contacts: median {frame.peptide_mhc_contacts.median():.0f}, "
          f"min {frame.peptide_mhc_contacts.min():.0f}")
    print(f"  complexes where the TCR touches the peptide at all: "
          f"{int(frame.tcr_touches_pmhc.sum())}/{len(frame)}")
    print("\n=== C06 docking geometry ===")
    print(frame[["crossing_angle_deg", "tcr_height_over_peptide",
                 "contacts_a_pep", "contacts_b_pep", "cdr3_share_a", "cdr3_share_b"]]
          .describe().loc[["min", "mean", "max", "std"]].round(2).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
