#!/usr/bin/env python
"""Write one Boltz-2 YAML per test TCR: a four-chain TCR-peptide-MHC complex.

Why structure at all: every sequence-neighbourhood signal saturates near 0.84
and every supervised layer on five epitopes loses to variance
(reports/02_signal_landscape.md, reports/03_loeo_fusion.md). Structure is the
one remaining source that is neighbour-free -- it reads this TCR against this
peptide directly -- so it can in principle rank the 34-83% of each cohort that
sits in no cluster. The organisers' own AlphaFold3 min-PAE baseline reaches
0.787 alone, below density, which is why the interesting question is not
whether structure wins but whether it is *complementary*.

Chain layout, fixed across all complexes so downstream indices are stable:
  A = TCR alpha variable domain   (`full_seq_alpha_aa`, ~112 aa)
  B = TCR beta variable domain    (`full_seq_beta_aa`, ~114 aa)
  M = HLA alpha1-alpha2 platform  (allele-specific, 182 aa)
  P = peptide                     (9-10 aa)

Roughly 416 residues total. Beta-2 microglobulin and the MHC alpha3 domain are
left out: they sit away from the TCR-peptide interface that carries the signal,
and dropping them cuts the complex by a third.

Predictions run single-sequence (`msa: empty`). The TCR and MHC folds are
heavily represented in training data -- the pilot reached complex_pLDDT 0.96
without an MSA -- and it keeps the full set to about 10 GPU-hours instead of
adding an MSA-server round trip per complex.

Run: python scripts/08_boltz_inputs.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

HLA_FASTA = ROOT / "data/external/hla/A_prot.fasta"
OUTDIR = ROOT / "data/interim/boltz/inputs"

# IMGT/HLA protein sequences carry a 24-residue signal peptide; the
# alpha1-alpha2 platform is the next 182 residues.
SIGNAL_PEPTIDE_LEN = 24
PLATFORM_LEN = 182

# The competition uses two HLA alleles. Four-field IMGT names pin the exact
# protein sequence rather than trusting a two-field prefix match.
ALLELE_RECORDS = {"A*02:01": "A*02:01:01:01", "A*01:01": "A*01:01:01:01"}


def load_hla_platforms(path: Path) -> dict[str, str]:
    """Allele -> alpha1-alpha2 sequence, from the IMGT/HLA protein FASTA."""
    sequences: dict[str, str] = {}
    name, buffer = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(buffer)
            name, buffer = line.split()[1], []
        else:
            buffer.append(line.strip())
    if name:
        sequences[name] = "".join(buffer)

    platforms = {}
    for allele, record in ALLELE_RECORDS.items():
        if record not in sequences:
            raise KeyError(f"{record} not found in {path}")
        full = sequences[record]
        platforms[allele] = full[SIGNAL_PEPTIDE_LEN : SIGNAL_PEPTIDE_LEN + PLATFORM_LEN]
    return platforms


def yaml_for(row: pd.Series, platform: str) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        f"  - protein:\n      id: A\n      sequence: {row.full_seq_alpha_aa}\n      msa: empty\n"
        f"  - protein:\n      id: B\n      sequence: {row.full_seq_beta_aa}\n      msa: empty\n"
        f"  - protein:\n      id: M\n      sequence: {platform}\n      msa: empty\n"
        f"  - protein:\n      id: P\n      sequence: {row.epitope_aa}\n      msa: empty\n"
    )


def main() -> int:
    test = load_test()
    platforms = load_hla_platforms(HLA_FASTA)

    unknown = set(test["hla_long"]) - set(platforms)
    if unknown:
        raise ValueError(f"no HLA platform sequence for {sorted(unknown)}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    written = 0
    for row in test.itertuples(index=False):
        # tcr_id is filesystem-safe (alphanumerics and underscores) and becomes
        # the prediction's name, which is how outputs are joined back.
        target = OUTDIR / f"{row.tcr_id}.yaml"
        if target.exists():
            continue
        target.write_text(yaml_for(row, platforms[row.hla_long]))
        written += 1

    total = len(list(OUTDIR.glob("*.yaml")))
    print(f"wrote {written} new YAML(s); {total} total in {OUTDIR.relative_to(ROOT)}")
    if total != len(test):
        raise ValueError(f"expected {len(test)} inputs, found {total}")
    print("\nnext: scripts/09_boltz_run.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
