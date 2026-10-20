#!/usr/bin/env python
"""Extract interface-confidence features from the Boltz-2 predictions.

The organisers' structural baseline is AlphaFold3 min-PAE at the TCR-peptide
interface, which reaches 0.787 alone. The same quantity is recomputed here from
Boltz's PAE matrix, alongside Boltz's own chain-pair ipTM, which is a
confidence read the AF3 baseline did not use.

Chain order is fixed by scripts/08_boltz_inputs.py -- A=TCR alpha, B=TCR beta,
M=MHC platform, P=peptide -- and the PAE matrix is laid out in that order, so
boundaries come from the per-complex chain lengths in the input YAML. Alpha and
beta lengths vary between TCRs, so they cannot be hard-coded.

Works on a partially finished run, which is how it gets validated before the
full 10-hour fold completes.

Run: python scripts/10_boltz_features.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

INPUTS = ROOT / "data/interim/boltz/inputs"
PREDICTIONS = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
OUT = ROOT / "data/interim/boltz_features.csv"

# Boltz indexes chains by position, matching the YAML order.
ALPHA, BETA, MHC, PEPTIDE = 0, 1, 2, 3


def chain_bounds(yaml_path: Path) -> list[tuple[int, int]]:
    """Half-open [start, stop) residue ranges per chain, in YAML order."""
    lengths = [len(s) for s in re.findall(r"sequence: (\S+)", yaml_path.read_text())]
    bounds, offset = [], 0
    for length in lengths:
        bounds.append((offset, offset + length))
        offset += length
    return bounds


def interface_stats(pae: np.ndarray, first: tuple[int, int], second: tuple[int, int]) -> dict:
    """PAE summaries over one inter-chain block, symmetrised.

    min-PAE is the organisers' statistic: the single most confidently placed
    residue pair across the interface. The mean and the 10th percentile are
    added because a min over a few thousand entries is a noisy order statistic.
    """
    block = np.concatenate(
        [
            pae[first[0] : first[1], second[0] : second[1]].ravel(),
            pae[second[0] : second[1], first[0] : first[1]].ravel(),
        ]
    )
    return {
        "min": float(block.min()),
        "p10": float(np.percentile(block, 10)),
        "mean": float(block.mean()),
    }


def extract(tcr_id: str) -> dict | None:
    directory = PREDICTIONS / tcr_id
    pae_file = directory / f"pae_{tcr_id}_model_0.npz"
    confidence_file = directory / f"confidence_{tcr_id}_model_0.json"
    plddt_file = directory / f"plddt_{tcr_id}_model_0.npz"
    if not (pae_file.exists() and confidence_file.exists()):
        return None

    pae = np.load(pae_file)["pae"]
    confidence = json.loads(confidence_file.read_text())
    bounds = chain_bounds(INPUTS / f"{tcr_id}.yaml")
    if pae.shape[0] != bounds[-1][1]:
        raise ValueError(f"{tcr_id}: PAE is {pae.shape[0]} but chains sum to {bounds[-1][1]}")

    row: dict[str, object] = {"tcr_id": tcr_id}

    for label, (i, j) in (
        ("a_pep", (ALPHA, PEPTIDE)),
        ("b_pep", (BETA, PEPTIDE)),
        ("a_mhc", (ALPHA, MHC)),
        ("b_mhc", (BETA, MHC)),
        ("ab", (ALPHA, BETA)),
        ("mhc_pep", (MHC, PEPTIDE)),
    ):
        stats = interface_stats(pae, bounds[i], bounds[j])
        for key, value in stats.items():
            # Negated so that larger always means "more confident interface",
            # matching the direction of every other feature in the project.
            row[f"pae_{label}_{key}"] = -value
        row[f"iptm_{label}"] = confidence["pair_chains_iptm"][str(i)][str(j)]

    # The organisers' baseline is the better of the two TCR-peptide interfaces.
    row["pae_tcr_pep_min"] = max(row["pae_a_pep_min"], row["pae_b_pep_min"])
    row["iptm_tcr_pep_max"] = max(row["iptm_a_pep"], row["iptm_b_pep"])
    row["iptm_tcr_pep_mean"] = 0.5 * (row["iptm_a_pep"] + row["iptm_b_pep"])

    for key in ("confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt"):
        row[key] = confidence[key]
    row["complex_ipde"] = -confidence["complex_ipde"]

    if plddt_file.exists():
        plddt = np.load(plddt_file)["plddt"]
        row["plddt_peptide"] = float(plddt[bounds[PEPTIDE][0] : bounds[PEPTIDE][1]].mean())
    return row


def main() -> int:
    ids = sorted(p.stem for p in INPUTS.glob("*.yaml"))
    rows = [r for r in (extract(i) for i in ids) if r is not None]
    if not rows:
        sys.exit(f"no finished predictions under {PREDICTIONS}")

    features = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}  {len(features)}/{len(ids)} complexes")
    print(f"columns: {[c for c in features.columns if c != 'tcr_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
