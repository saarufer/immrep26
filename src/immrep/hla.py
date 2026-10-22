"""Allele-specific HLA class I alpha1-alpha2 platform sequences.

Shared by the structure pipelines so every complex uses the same groove for a
given allele. IMGT/HLA protein records carry a 24-residue signal peptide; the
platform is the next 182 residues.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FASTA = ROOT / "data/external/hla/A_prot.fasta"
SIGNAL_PEPTIDE_LEN = 24
PLATFORM_LEN = 182
ALLELE_RECORDS = {"A*02:01": "A*02:01:01:01", "A*01:01": "A*01:01:01:01"}


def load_hla_platforms(path: Path | None = None) -> dict[str, str]:
    path = path or FASTA
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
    out = {}
    for allele, record in ALLELE_RECORDS.items():
        if record not in sequences:
            raise KeyError(f"{record} not found in {path}")
        full = sequences[record]
        out[allele] = full[SIGNAL_PEPTIDE_LEN : SIGNAL_PEPTIDE_LEN + PLATFORM_LEN]
    return out
