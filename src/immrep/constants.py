"""IMMREP26 competition constants.

Counts are the per-epitope TCR totals stated in the competition data description.
Public epitopes have functional-validation labels published via TCRvdb
(2026-10-01 release); the private three are held back for the final ranking.
"""

from __future__ import annotations

PUBLIC_EPITOPES: dict[str, dict[str, object]] = {
    "GILGFVFTL": {"hla": "A*02:01", "n": 670, "gene": "IAV M1"},
    "YLQPRTFLL": {"hla": "A*02:01", "n": 419, "gene": "SARS-CoV-2 S"},
    "TTDPSFLGRY": {"hla": "A*01:01", "n": 387, "gene": "SARS-CoV-2 ORF1ab"},
    "NLVPMVATV": {"hla": "A*02:01", "n": 366, "gene": "CMV pp65"},
    "GLCTLVAML": {"hla": "A*02:01", "n": 192, "gene": "EBV BMLF1"},
}

PRIVATE_EPITOPES: dict[str, dict[str, object]] = {
    "CINGVCWTV": {"hla": "A*02:01", "n": 243, "gene": "HCV NS3"},
    "LLWNGPMAV": {"hla": "A*02:01", "n": 229, "gene": "YFV NS4B"},
    "LTDEMIAQY": {"hla": "A*01:01", "n": 130, "gene": "SARS-CoV-2 ORF1ab"},
}

ALL_EPITOPES = {**PUBLIC_EPITOPES, **PRIVATE_EPITOPES}

N_TEST_ROWS = 2636
N_PUBLIC = sum(int(v["n"]) for v in PUBLIC_EPITOPES.values())   # 2034
N_PRIVATE = sum(int(v["n"]) for v in PRIVATE_EPITOPES.values()) # 602

# Positive (validated) rates reported in the bioRxiv preprint for the two
# epitopes released in the 2025-05-01 TCRvdb version. The ~29-point spread is
# why pooled-across-epitope calibration matters if the official metric turns
# out to be a single global AUC rather than a per-epitope mean.
PREPRINT_VALIDATION_RATES = {
    "YLQPRTFLL": {"screened": 421, "validated": 171},  # 40.6%
    "GLCTLVAML": {"screened": 193, "validated": 134},  # 69.4%
}

# Organiser baselines on the public leaderboard (metric not yet confirmed).
BASELINES = {"tcrdist3": 0.824, "alphafold3_min_pae": 0.787}

VDJDB_VERSION = "2026-06-03"
TCRVDB_VERSION = "2026-10-01"
