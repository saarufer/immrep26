#!/usr/bin/env python
"""Docking consistency across diffusion samples, from coordinates.

Confidence turned out not to be the orthogonal signal it looked like: every
Boltz confidence read, including `iptm_mhc_pep` at an interface containing no
TCR at all, scores ~0.65-0.76, which says confidence mostly reports how
canonical the TCR sequence is (reports/00_experiment_log.md, E14).

Consistency is a different quantity. A TCR-pMHC pair with one real binding mode
should dock to the same pose every time the diffusion head is sampled; a pair
with no real mode should land differently each time. Boltz runs its trunk once
and samples only the diffusion head, so five samples cost 1.36x one sample.

Measurement: superpose each pair of samples on the MHC and peptide C-alphas --
the invariant part of the complex -- then take the RMSD of the TCR C-alphas.
That asks "does the TCR sit in the same place relative to the peptide groove",
which is the geometric question, rather than "is the whole complex
reproducible", which the conserved folds would answer trivially yes.

Run: python scripts/17_boltz_consistency.py
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PREDICTIONS = ROOT / "data/interim/boltz/out5/boltz_results_inputs/predictions"
OUT = ROOT / "data/interim/boltz_consistency.csv"

TCR_CHAINS = ("A", "B")
ANCHOR_CHAINS = ("M", "P")


def read_ca(path: Path) -> dict[str, np.ndarray]:
    """C-alpha coordinates per chain, in file order."""
    coords: dict[str, list[list[float]]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        coords.setdefault(line[21], []).append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    return {chain: np.asarray(xyz, dtype=np.float64) for chain, xyz in coords.items()}


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation and translation that best superpose `mobile` onto `target`."""
    mobile_centre, target_centre = mobile.mean(0), target.mean(0)
    correlation = (mobile - mobile_centre).T @ (target - target_centre)
    u, _, vt = np.linalg.svd(correlation)
    sign = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, sign]) @ u.T
    return rotation, target_centre - rotation @ mobile_centre


def pair_rmsd(first: dict, second: dict) -> tuple[float, float]:
    """TCR and peptide RMSD between two samples, after anchoring on MHC+peptide."""
    anchor_a = np.vstack([first[c] for c in ANCHOR_CHAINS])
    anchor_b = np.vstack([second[c] for c in ANCHOR_CHAINS])
    rotation, translation = kabsch(anchor_a, anchor_b)

    tcr_a = np.vstack([first[c] for c in TCR_CHAINS]) @ rotation.T + translation
    tcr_b = np.vstack([second[c] for c in TCR_CHAINS])
    tcr_rmsd = float(np.sqrt(((tcr_a - tcr_b) ** 2).sum(axis=1).mean()))

    # Peptide alone, anchored on the MHC only: does the peptide itself shift in
    # the groove between samples.
    rotation, translation = kabsch(first["M"], second["M"])
    peptide_a = first["P"] @ rotation.T + translation
    peptide_rmsd = float(np.sqrt(((peptide_a - second["P"]) ** 2).sum(axis=1).mean()))
    return tcr_rmsd, peptide_rmsd


def summarise(directory: Path, tcr_id: str) -> dict | None:
    models = sorted(directory.glob(f"{tcr_id}_model_*.pdb"))
    if len(models) < 2:
        return None
    structures = [read_ca(p) for p in models]
    needed = set(TCR_CHAINS) | set(ANCHOR_CHAINS)
    if any(not needed <= set(s) for s in structures):
        return None

    tcr, peptide = [], []
    for a, b in itertools.combinations(structures, 2):
        t, p = pair_rmsd(a, b)
        tcr.append(t)
        peptide.append(p)

    confidences = []
    for path in sorted(directory.glob(f"confidence_{tcr_id}_model_*.json")):
        confidences.append(json.loads(path.read_text()))

    row = {
        "tcr_id": tcr_id,
        "n_samples": len(structures),
        # Negated so that larger always means "more consistent", matching the
        # direction of every other feature in the project.
        "consist_tcr_rmsd_mean": -float(np.mean(tcr)),
        "consist_tcr_rmsd_max": -float(np.max(tcr)),
        "consist_pep_rmsd_mean": -float(np.mean(peptide)),
        "consist_pep_rmsd_max": -float(np.max(peptide)),
    }
    if confidences:
        for key in ("iptm", "confidence_score", "complex_plddt"):
            values = [c[key] for c in confidences]
            row[f"ms_{key}_mean"] = float(np.mean(values))
            row[f"ms_{key}_sd"] = -float(np.std(values))
        a_pep = [c["pair_chains_iptm"]["0"]["3"] for c in confidences]
        b_pep = [c["pair_chains_iptm"]["1"]["3"] for c in confidences]
        row["ms_tcr_pep_iptm_mean"] = float(np.mean(a_pep + b_pep))
        row["ms_tcr_pep_iptm_sd"] = -float(np.std(a_pep + b_pep))
    return row


def main() -> int:
    if not PREDICTIONS.exists():
        sys.exit(f"missing {PREDICTIONS}; run scripts/16_boltz_multisample.sh first")

    rows = []
    directories = sorted(p for p in PREDICTIONS.iterdir() if p.is_dir())
    for index, directory in enumerate(directories, 1):
        row = summarise(directory, directory.name)
        if row is not None:
            rows.append(row)
        if index % 250 == 0:
            print(f"  {index}/{len(directories)} ...")

    if not rows:
        sys.exit("no complexes with >=2 samples yet")

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}  {len(frame)} complexes")
    print(
        frame[["consist_tcr_rmsd_mean", "consist_pep_rmsd_mean", "ms_tcr_pep_iptm_sd"]]
        .describe()
        .loc[["min", "mean", "max", "std"]]
        .round(3)
        .to_string()
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
