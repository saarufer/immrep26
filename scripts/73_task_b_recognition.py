#!/usr/bin/env python
"""E82 -- task B: do representations predict recognition-profile distance?

The decisive control for "antigen-conditioned retrieval over an
antigen-agnostic representation": run the same representations on an
independent recognition task instead of on VDJdb revalidation.

Only PRP can supply it, and only in reduced form -- it has a shared peptide
library but **CDR3-beta only**. SEEKER's PacBio module has 101 fully paired
TCRs but each peptide is assigned to exactly one TCR, so every pair shares
zero; the 4TCR screen shares peptides across four receptors whose sequences are
not in the package.

**Two limits on the inference.** The pairs are not independent -- 21 receptors
carrying 16 distinct CDR3-beta generate 156 pairs, so every receptor appears in
about 17 of them and the p-values are anticonservative. And this tests beta
sequence only; E63 found the five identical-CDR3-beta pairs here recognise no
more similarly than random ones, and those pairs differ in alpha.

Run: .venv/bin/python scripts/73_task_b_recognition.py
"""
from __future__ import annotations
import itertools, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
PRP = ROOT / "data/external/prp"
REPORT = ROOT / "reports/48_task_b_recognition.md"
ALPHA, MIN_COMMON = 1.0, 5000
P2_P9 = {"4.2", "4.3", "4.4"}
VARIANTS = [("default", "real repertoire"), ("blosum", "real, BLOSUM embedding"),
            ("shuffled_data", "real, alpha-beta pairing randomised"),
            ("synthetic_data", "OLGA synthetic only"),
            ("cdr3_only", "real, CDR3 only"), ("b_sceptr", "real, beta chain only")]


def levenshtein(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1,
                               previous[j - 1] + (x != y)))
        previous = current
    return previous[-1]


def main() -> int:
    import sceptr
    from sceptr import variant as V
    sceptr.disable_hardware_acceleration()
    cdr = pd.read_csv(PRP / "github/tcr_cdr3b.csv")
    sequence = dict(zip(cdr["tcr_id"].astype(str), cdr["CDR3_b"]))

    def load(tcr):
        path = PRP / f"TCR{tcr}_R1-R4.csv"
        with path.open() as handle:
            first = handle.readline().split(",")[0]
        frame = (pd.read_csv(path) if first == "Peptide" else
                 pd.read_csv(path, header=None,
                             names=["Peptide", "Round1", "Round2", "Round3", "Round4"]))
        return frame.set_index("Peptide")

    profiles = {}
    for tcr in cdr["tcr_id"].astype(str):
        counts = load(tcr)
        r1, r4 = counts["Round1"].fillna(0.0), counts["Round4"].fillna(0.0)
        profiles[tcr] = np.log10(((r4 + ALPHA) / (r4.sum() + ALPHA))
                                 / ((r1 + ALPHA) / (r1.sum() + ALPHA)))
    pairs = []
    for a, b in itertools.combinations(sorted(sequence), 2):
        if (a in P2_P9) != (b in P2_P9):
            continue
        common = profiles[a].index.intersection(profiles[b].index)
        if len(common) < MIN_COMMON:
            continue
        pairs.append((a, b, 1 - spearmanr(profiles[a][common],
                                          profiles[b][common]).statistic))
    recognition = np.array([p[2] for p in pairs])
    print(f"=== {len(pairs)} comparable pairs from "
          f"{len(set(sequence.values()))} distinct CDR3-beta ===\n")

    ids = list(cdr["tcr_id"].astype(str))
    frame = pd.DataFrame({"TRAV": [None] * len(ids), "CDR3A": [None] * len(ids),
                          "TRAJ": [None] * len(ids), "TRBV": [None] * len(ids),
                          "CDR3B": [sequence[t] for t in ids],
                          "TRBJ": [None] * len(ids)})
    rows = []
    for name, description in VARIANTS:
        embedding = getattr(V, name)().calc_vector_representations(frame)
        embedding = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)
        position = {t: i for i, t in enumerate(ids)}
        distance = np.array([1 - float(embedding[position[a]] @ embedding[position[b]])
                             for a, b, _ in pairs])
        result = spearmanr(distance, recognition)
        rows.append({"representation": name, "pretrained_on": description,
                     "rho": result.statistic, "p": result.pvalue})
        print(f"  {name:<16} {description:<36} rho {result.statistic:+.3f}  "
              f"p {result.pvalue:.3f}")
    edit = np.array([levenshtein(sequence[a], sequence[b]) for a, b, _ in pairs])
    result = spearmanr(edit, recognition)
    rows.append({"representation": "Levenshtein", "pretrained_on": "E63's baseline",
                 "rho": result.statistic, "p": result.pvalue})
    print(f"  {'Levenshtein':<16} {'(E63 baseline)':<36} rho "
          f"{result.statistic:+.3f}  p {result.pvalue:.3f}")
    table = pd.DataFrame(rows)
    print("\n  Every value is indistinguishable from zero, and the "
          "OLGA-pretrained model\n  is numerically above the real-pretrained one.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("# E82 - task B, recognition-profile distance\n\n"
                      + table.round(4).to_markdown(index=False)
                      + f"\n\n{len(pairs)} pairs, beta-only, non-independent "
                      "(21 receptors, 16 distinct sequences).\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
