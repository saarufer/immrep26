#!/usr/bin/env python
"""E63 -- is recognition distance the same thing as sequence distance?

PRP's Dryad package carries peptide counts and nothing else; the TCR sequences
live in the companion GitHub repo as `data/tcr_cdr3b.csv`, 487 bytes. Reading
it gives the finding this experiment is built on:

    CASSVAYYSTDTQYF  ->  135.7, 19.2C1, 28.1
    CASSVATYSTDTQYF  ->  19.2C2, 4.3
    CASSVGLYSTDTQYF  ->  19.2C5, 4.1
    CASSVGTYSTDTQYF  ->  19.2C4, 8.4

**Twenty-one screened receptors carry sixteen distinct CDR3β.** Five comparable
pairs have a CDR3β Hamming distance of exactly zero while being different
receptors -- the engineered 19.2 variants were built to carry another patient
TCR's β. If their recognition profiles diverge, that is the strongest possible
form of "sequence close, recognition far", because the sequence distance is not
merely small, it is zero.

**The limitation that bounds every claim here: only CDR3β is available.** No
alpha chain, no V/J. So "sequence distance zero" means "zero on the axis that
can be measured, alpha unknown", and a divergence could be an alpha effect
rather than a failure of sequence distance in general. Stated once here and
carried into the conclusion.

Two further constraints:

  * the two anchor libraries are different panels (E60: 272,888 and 151,914
    peptides sharing 256), so TCRs are compared only within their own;
  * profile distance has a noise floor. E62 was nearly misread because a
    cross-receptor correlation was compared against nothing; here the
    within-receptor split-half reliability is computed **first** and every
    distance is read against it.

Enrichment is round 4 against round 1, not against the naive library, because
E60 measured that only 8.9% of selected peptides appear in their own input
library.

Run: .venv/bin/python scripts/61_recognition_vs_sequence.py
"""

from __future__ import annotations

import itertools
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
PRP = ROOT / "data/external/prp"
REPORT = ROOT / "reports/36_recognition_vs_sequence.md"
P2_P9 = {"4.2", "4.3", "4.4"}
MIN_COMMON = 5_000
ALPHA = 1.0
SEED = 0


def load_counts(tcr: str) -> pd.DataFrame:
    path = PRP / f"TCR{tcr}_R1-R4.csv"
    with path.open() as handle:
        first = handle.readline().split(",")[0]
    if first == "Peptide":
        return pd.read_csv(path).set_index("Peptide")
    return pd.read_csv(path, header=None,
                       names=["Peptide", "Round1", "Round2", "Round3",
                              "Round4"]).set_index("Peptide")


def profile(counts: pd.DataFrame) -> pd.Series:
    r1 = counts["Round1"].fillna(0.0)
    r4 = counts["Round4"].fillna(0.0)
    return np.log10(((r4 + ALPHA) / (r4.sum() + ALPHA))
                    / ((r1 + ALPHA) / (r1.sum() + ALPHA)))


def hamming(a: str, b: str) -> float:
    if len(a) != len(b):
        return float("nan")
    return float(sum(x != y for x, y in zip(a, b)))


def levenshtein(a: str, b: str) -> int:
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1,
                               previous[j - 1] + (x != y)))
        previous = current
    return previous[-1]


def main() -> int:
    rng = np.random.default_rng(SEED)
    cdr = pd.read_csv(PRP / "github/tcr_cdr3b.csv")
    sequence = dict(zip(cdr["tcr_id"].astype(str), cdr["CDR3_b"]))
    counts = {t: load_counts(t) for t in sequence}
    profiles = {t: profile(c) for t, c in counts.items()}

    print("=== noise floor first: does a receptor agree with itself? ===")
    print("  Split a receptor's round-4 reads in two by a binomial draw and")
    print("  rebuild the profile from each half. Any profile distance between")
    print("  receptors has to be read against this.\n")
    ceiling = {}
    for tcr, frame in counts.items():
        r1 = frame["Round1"].fillna(0.0)
        r4 = frame["Round4"].fillna(0.0)
        keep = (r4 > 0) & (r1 > 0)
        draws = []
        for _ in range(5):
            left = rng.binomial(r4[keep].astype(int), 0.5)
            right = r4[keep].astype(int).to_numpy() - left
            base = (r1[keep] + ALPHA) / (r1[keep].sum() + ALPHA)
            a = np.log10(((left + ALPHA) / (left.sum() + ALPHA)) / base)
            b = np.log10(((right + ALPHA) / (right.sum() + ALPHA)) / base)
            draws.append(spearmanr(a, b).statistic)
        ceiling[tcr] = float(np.mean(draws))
    mean_ceiling = float(np.mean(list(ceiling.values())))
    print(f"  split-half rho per receptor: {min(ceiling.values()):.3f} to "
          f"{max(ceiling.values()):.3f}, mean **{mean_ceiling:.3f}**")

    print("\n=== every comparable pair ===")
    rows = []
    for a, b in itertools.combinations(sorted(sequence), 2):
        if (a in P2_P9) != (b in P2_P9):
            continue
        common = profiles[a].index.intersection(profiles[b].index)
        if len(common) < MIN_COMMON:
            continue
        rho = spearmanr(profiles[a][common], profiles[b][common]).statistic
        rows.append({
            "a": a, "b": b, "library": "P2_P9" if a in P2_P9 else "P2_P8",
            "common_peptides": len(common),
            "cdr3b_hamming": hamming(sequence[a], sequence[b]),
            "cdr3b_levenshtein": levenshtein(sequence[a], sequence[b]),
            "profile_rho": float(rho),
            "profile_distance": 1.0 - float(rho),
            "ceiling": (ceiling[a] + ceiling[b]) / 2,
        })
    table = pd.DataFrame(rows)
    table["rho_over_ceiling"] = table["profile_rho"] / table["ceiling"]
    print(f"  {len(table)} pairs with >= {MIN_COMMON:,} common peptides")

    zero = table[table["cdr3b_levenshtein"] == 0]
    rest = table[table["cdr3b_levenshtein"] > 0]
    print("\n=== the five identical-CDR3b pairs ===")
    print(zero[["a", "b", "common_peptides", "profile_rho", "ceiling",
                "rho_over_ceiling"]].round(3).to_string(index=False))
    print(f"\n  identical CDR3b: median profile rho "
          f"{zero['profile_rho'].median():.3f}")
    print(f"  everything else: median profile rho "
          f"{rest['profile_rho'].median():.3f} "
          f"(range {rest['profile_rho'].min():.3f} to "
          f"{rest['profile_rho'].max():.3f})")

    print("\n=== does sequence distance explain recognition distance at all? ===")
    for column in ("cdr3b_levenshtein", "cdr3b_hamming"):
        usable = table.dropna(subset=[column])
        rho = spearmanr(usable[column], usable["profile_distance"]).statistic
        print(f"  {column:<20} vs profile distance: rho {rho:+.3f} "
              f"(n={len(usable)})")

    print("\n=== counterexamples in both directions ===")
    close_far = table[(table["cdr3b_levenshtein"] <= 1)].nsmallest(
        3, "profile_rho")
    far_close = table[(table["cdr3b_levenshtein"] >= 4)].nlargest(
        3, "profile_rho")
    print("  sequence close, recognition far:")
    print(close_far[["a", "b", "cdr3b_levenshtein", "profile_rho",
                     "rho_over_ceiling"]].round(3).to_string(index=False))
    print("  sequence far, recognition close:")
    print(far_close[["a", "b", "cdr3b_levenshtein", "profile_rho",
                     "rho_over_ceiling"]].round(3).to_string(index=False))

    print("\n=== reading it ===")
    gap = zero["profile_rho"].median() - rest["profile_rho"].median()
    print(f"  Identical CDR3b buys {gap:+.3f} of profile correlation over a")
    print(f"  random comparable pair, against a reliability ceiling of "
          f"{mean_ceiling:.3f}.")
    if zero["profile_rho"].median() < 0.5 * mean_ceiling:
        print("  Receptors with the SAME CDR3beta recognise differently: the")
        print("  sequence axis has zero resolution here while the recognition")
        print("  axis separates them. Alpha is unobserved, so this bounds how")
        print("  much CDR3beta alone can carry, rather than proving sequence")
        print("  distance fails in general.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E63 - recognition distance against sequence distance (PRP)\n\n"
        f"21 screened receptors carry **16 distinct CDR3β**. Only CDR3β is "
        f"available -- no alpha, no V/J -- so every claim is bounded to that "
        f"axis. Within-receptor split-half reliability (the ceiling): "
        f"**{mean_ceiling:.3f}**.\n\n## The five identical-CDR3β pairs\n\n"
        + zero.round(4).to_markdown(index=False)
        + "\n\n## All comparable pairs\n\n"
        + table.round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
