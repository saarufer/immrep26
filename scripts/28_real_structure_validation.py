#!/usr/bin/env python
"""Do the predicted contacts match crystallographic ones? (plan card D09)

E28 built contact fingerprints from Boltz predictions and E31 used them to argue
that the discordant negatives differ at peptide-contacting positions. Both rest
on predicted interfaces, and E31 said so as its main limitation: a predicted
contact is not a measured one. This checks it against experiment.

The PDB holds solved TCR-pMHC complexes for four of the competition's epitopes.
Taking the TCR, peptide and MHC sequences straight out of those crystals,
folding them through the *same* Boltz pipeline used for the competition, and
comparing the resulting fingerprint to the one measured in the crystal gives a
direct read on whether the fingerprints mean anything. 15 complexes is a small
sample, so the question asked is deliberately coarse: do predicted and measured
fingerprints agree better than chance, and is the agreement good enough that
differences *between* TCRs are interpretable.

Stage 1 writes Boltz inputs; run scripts/29_fold_real.sh; stage 2 compares.

Run: python scripts/28_real_structure_validation.py [--compare]
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PDB_DIR = ROOT / "data/external/pdb/cif"
INPUT_DIR = ROOT / "data/interim/realval/inputs"
PRED_DIR = ROOT / "data/interim/realval/out/boltz_results_inputs/predictions"
TRUTH = ROOT / "data/interim/realval/crystal_fingerprints.json"
OUT = ROOT / "reports/09_real_structure_validation.md"

CONTACT_CUTOFF = 4.5
MAX_PEPTIDE_LEN = 10
EPITOPES = {
    "GILGFVFTL", "YLQPRTFLL", "TTDPSFLGRY", "NLVPMVATV",
    "GLCTLVAML", "CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY",
    # 9x7r, released 2025-10, the only TCR-pMHC complex found that postdates the
    # structure model. Its peptide is not a competition epitope; it is included
    # purely as a single out-of-distribution anchor, since every other complex
    # here predates Boltz-2 and so tests memorisation rather than prediction.
    "AMDLGIHKV",
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def parse_cif_atoms(path: Path) -> dict[str, dict]:
    """Heavy atoms per author chain, with per-residue sequence and coordinates."""
    text = gzip.open(path, "rt", errors="ignore").read()
    # Locate the atom_site loop and its column order.
    match = re.search(r"loop_\s*((?:_atom_site\.\S+\s*)+)", text)
    if not match:
        return {}
    columns = re.findall(r"_atom_site\.(\S+)", match.group(1))
    index = {name: i for i, name in enumerate(columns)}
    needed = ("group_PDB", "type_symbol", "label_comp_id", "auth_asym_id",
              "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z")
    if any(k not in index for k in needed):
        return {}

    start = match.end()
    chains: dict[str, dict] = {}
    for line in text[start:].splitlines():
        if line.startswith(("#", "loop_", "_")):
            if line.startswith("#") and chains:
                break
            continue
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < len(columns) or fields[index["group_PDB"]] != "ATOM":
            continue
        if fields[index["type_symbol"]] == "H":
            continue
        chain = fields[index["auth_asym_id"]]
        residue = fields[index["label_comp_id"]]
        try:
            number = int(fields[index["auth_seq_id"]])
            xyz = (float(fields[index["Cartn_x"]]), float(fields[index["Cartn_y"]]),
                   float(fields[index["Cartn_z"]]))
        except ValueError:
            continue
        entry = chains.setdefault(chain, {"xyz": [], "resnum": [], "residues": {}})
        entry["xyz"].append(xyz)
        entry["resnum"].append(number)
        entry["residues"].setdefault(number, residue)

    for entry in chains.values():
        entry["xyz"] = np.asarray(entry["xyz"], dtype=np.float32)
        entry["resnum"] = np.asarray(entry["resnum"], dtype=np.int32)
        entry["seq"] = "".join(
            THREE_TO_ONE.get(entry["residues"][n], "X") for n in sorted(entry["residues"])
        )
    return chains


def classify(seq: str) -> str:
    if seq in EPITOPES:
        return "peptide"
    if "GSHSMRY" in seq[:30] and len(seq) >= 250:
        return "mhc"
    if 100 <= len(seq) <= 260 and re.search(r"Y[FLIVM]C", seq) and re.search(r"FG.G[TS]", seq):
        return "tcr"
    return "other"


def crystal_fingerprint(chains: dict, tcr_ids: list[str], peptide_id: str) -> np.ndarray:
    """The same 60-dim layout as scripts/24: per-peptide-position contact count
    and minimum distance, for each of two TCR chains, plus a CDR3-share slot
    left at zero because crystal CDR3 bounds are not resolved here."""
    peptide = chains[peptide_id]
    positions = sorted(set(peptide["resnum"].tolist()))
    blocks = []
    for chain_id in tcr_ids[:2]:
        tree = cKDTree(chains[chain_id]["xyz"])
        counts = np.zeros(MAX_PEPTIDE_LEN, dtype=np.float32)
        minima = np.full(MAX_PEPTIDE_LEN, 20.0, dtype=np.float32)
        for slot, number in enumerate(positions[:MAX_PEPTIDE_LEN]):
            atoms = peptide["xyz"][peptide["resnum"] == number]
            if not len(atoms):
                continue
            near = tree.query_ball_point(atoms, CONTACT_CUTOFF)
            counts[slot] = float(sum(len(x) for x in near))
            distances, _ = tree.query(atoms, k=1)
            minima[slot] = float(np.min(distances))
        blocks.append(np.concatenate([counts, minima, np.zeros(MAX_PEPTIDE_LEN, np.float32)]))
    while len(blocks) < 2:
        blocks.append(np.concatenate([np.zeros(MAX_PEPTIDE_LEN), np.full(MAX_PEPTIDE_LEN, 20.0),
                                      np.zeros(MAX_PEPTIDE_LEN)]).astype(np.float32))
    return np.concatenate(blocks)


def stage_inputs() -> int:
    classes = pd.read_csv(ROOT / "data/external/pdb/entry_classes.csv")
    complexes = classes[(classes["n_TCR"] >= 2) & classes["has_mhc"] & classes["has_pep"]]
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    truth, written = {}, 0

    for row in complexes.itertuples(index=False):
        chains = parse_cif_atoms(PDB_DIR / f"{row.pdb}.cif.gz")
        if not chains:
            print(f"  {row.pdb}: could not parse atoms")
            continue
        roles = {c: classify(d["seq"]) for c, d in chains.items()}
        peptides = [c for c, r in roles.items() if r == "peptide"]
        tcrs = [c for c, r in roles.items() if r == "tcr"]
        mhcs = [c for c, r in roles.items() if r == "mhc"]
        if not peptides or len(tcrs) < 2 or not mhcs:
            print(f"  {row.pdb}: roles incomplete {sorted(set(roles.values()))}")
            continue

        # Pick the copy whose TCR chains actually touch the chosen peptide, since
        # several entries contain multiple copies of the complex.
        peptide_id = None
        best = -1
        for candidate in peptides:
            tree = cKDTree(chains[candidate]["xyz"])
            touch = sum(
                sum(len(x) for x in tree.query_ball_point(chains[t]["xyz"], CONTACT_CUTOFF))
                for t in tcrs
            )
            if touch > best:
                best, peptide_id = touch, candidate
        if best <= 0:
            print(f"  {row.pdb}: no TCR-peptide contact in any copy")
            continue
        peptide_tree = cKDTree(chains[peptide_id]["xyz"])
        contacting = sorted(
            tcrs,
            key=lambda t: -sum(
                len(x) for x in peptide_tree.query_ball_point(chains[t]["xyz"], CONTACT_CUTOFF)
            ),
        )[:2]

        name = f"{row.pdb}_{row.ep}"
        truth[name] = {
            "fingerprint": crystal_fingerprint(chains, contacting, peptide_id).tolist(),
            "epitope": row.ep,
            "tcr_chains": contacting,
            "peptide_chain": peptide_id,
            "peptide_seq": chains[peptide_id]["seq"],
        }
        # Boltz input in the project's fixed chain order: TCR a, TCR b, MHC, peptide.
        mhc_seq = chains[max(mhcs, key=lambda m: len(chains[m]["seq"]))]["seq"][:182]
        yaml = "version: 1\nsequences:\n"
        for chain_id, label in zip(contacting, ("A", "B")):
            yaml += (f"  - protein:\n      id: {label}\n"
                     f"      sequence: {chains[chain_id]['seq']}\n      msa: empty\n")
        yaml += (f"  - protein:\n      id: M\n      sequence: {mhc_seq}\n      msa: empty\n"
                 f"  - protein:\n      id: P\n"
                 f"      sequence: {chains[peptide_id]['seq']}\n      msa: empty\n")
        (INPUT_DIR / f"{name}.yaml").write_text(yaml)
        written += 1
        print(f"  {name}: TCR {contacting} peptide {peptide_id} "
              f"({len(chains[peptide_id]['seq'])} aa), {best} crystal contacts")

    TRUTH.parent.mkdir(parents=True, exist_ok=True)
    TRUTH.write_text(json.dumps(truth, indent=1))
    print(f"\nwrote {written} Boltz inputs and {len(truth)} crystal fingerprints")
    print("next: ./scripts/29_fold_real.sh   then rerun with --compare")
    return 0


def stage_compare() -> int:
    truth = json.loads(TRUTH.read_text())
    rows = []
    for name, record in truth.items():
        pdb = PRED_DIR / name / f"{name}_model_0.pdb"
        if not pdb.exists():
            continue
        coords: dict[str, list] = {}
        resnum: dict[str, list] = {}
        for line in pdb.read_text().splitlines():
            if not line.startswith("ATOM"):
                continue
            element = line[76:78].strip() or line[12:16].strip()[:1]
            if element == "H":
                continue
            chain = line[21]
            coords.setdefault(chain, []).append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54])))
            resnum.setdefault(chain, []).append(int(line[22:26]))
        if not {"A", "B", "P"} <= set(coords):
            continue
        chains = {
            c: {"xyz": np.asarray(coords[c], np.float32),
                "resnum": np.asarray(resnum[c], np.int32),
                "residues": {}}
            for c in coords
        }
        positions = sorted(set(chains["P"]["resnum"].tolist()))
        blocks = []
        for chain_id in ("A", "B"):
            tree = cKDTree(chains[chain_id]["xyz"])
            counts = np.zeros(MAX_PEPTIDE_LEN, np.float32)
            minima = np.full(MAX_PEPTIDE_LEN, 20.0, np.float32)
            for slot, number in enumerate(positions[:MAX_PEPTIDE_LEN]):
                atoms = chains["P"]["xyz"][chains["P"]["resnum"] == number]
                counts[slot] = float(sum(len(x) for x in tree.query_ball_point(atoms, CONTACT_CUTOFF)))
                distances, _ = tree.query(atoms, k=1)
                minima[slot] = float(np.min(distances))
            blocks.append(np.concatenate([counts, minima, np.zeros(MAX_PEPTIDE_LEN, np.float32)]))
        predicted = np.concatenate(blocks)
        measured = np.asarray(record["fingerprint"], dtype=float)

        # Correlate each block separately. Pooling contact counts with minimum
        # distances into one Spearman would inflate it: the two blocks are
        # inversely related and occupy different value ranges, so they separate
        # in the rank order of *both* vectors regardless of whether the
        # per-position pattern agrees. Only the within-block numbers are valid.
        counts_slots = np.r_[0:MAX_PEPTIDE_LEN, 30:30 + MAX_PEPTIDE_LEN]
        dist_slots = np.r_[MAX_PEPTIDE_LEN:2 * MAX_PEPTIDE_LEN,
                           30 + MAX_PEPTIDE_LEN:30 + 2 * MAX_PEPTIDE_LEN]
        rho_counts = spearmanr(predicted[counts_slots], measured[counts_slots]).correlation
        rho_dist = spearmanr(predicted[dist_slots], measured[dist_slots]).correlation
        touched_pred = set(np.where(predicted[counts_slots] > 0)[0])
        touched_real = set(np.where(measured[counts_slots] > 0)[0])
        jaccard = (len(touched_pred & touched_real) / len(touched_pred | touched_real)
                   if (touched_pred | touched_real) else np.nan)
        rows.append({
            "complex": name, "epitope": record["epitope"],
            "rho_counts": rho_counts, "rho_mindist": rho_dist,
            "jaccard_touched_positions": jaccard,
            "pred_contacts": float(predicted[counts_slots].sum()),
            "real_contacts": float(measured[counts_slots].sum()),
        })

    frame = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print("=== predicted vs crystallographic contact fingerprints ===")
    print(frame.round(3).to_string(index=False))
    print(f"\n  median Spearman (contact counts per peptide position): {frame.rho_counts.median():.3f}")
    print(f"  median Spearman (minimum distance per position):      {frame.rho_mindist.median():.3f}")
    print(f"  median Jaccard of touched positions: {frame.jaccard_touched_positions.median():.3f}")
    print(f"  predicted vs real total contacts: "
          f"{frame.pred_contacts.mean():.0f} vs {frame.real_contacts.mean():.0f}")
    print("\n  cross-complex check: does the prediction order complexes the same way?")
    if len(frame) >= 4:
        rho = spearmanr(frame.pred_contacts, frame.real_contacts).correlation
        print(f"    Spearman(total predicted contacts, total real contacts) = {rho:+.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Predicted against crystallographic contacts\n\n"
        "TCR, peptide and MHC sequences taken from solved complexes and folded\n"
        "through the same Boltz pipeline used for the competition, then compared to\n"
        "the interface measured in the crystal.\n\n" + frame.round(3).to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(stage_compare() if "--compare" in sys.argv else stage_inputs())
