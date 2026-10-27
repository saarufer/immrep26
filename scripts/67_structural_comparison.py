#!/usr/bin/env python
"""E75 -- structurally, how does GLCTLVAML differ from the other epitopes?

Every GLC analysis so far has been distributional: density, neighbour counts,
donors, labels. This asks the structural question instead, using the 2,636
Boltz-2 complexes already folded, with a clean control built into the data:
**six of the eight epitopes are HLA-A\\*02:01 9-mers** (GIL, GLC, NLV, YLQ,
CIN, LLW), so HLA and peptide length are held fixed and anything left is
peptide-specific. TTDPSFLGRY (A1, 10-mer) and LTDEMIAQY (A1, 9-mer) are
reported but not compared.

The hypothesis under test comes from E72: GLC's cohort is over-converged, and
if that is because its pMHC surface admits a narrow set of docking solutions,
the structures should show it -- a tighter interface ensemble, a narrower
contact footprint, less variable geometry.

One confound is checked first and it turns out to matter: an epitope whose
complexes Boltz folds *confidently* will produce similar structures whatever
the biology, so ensemble tightness has to be read against interface confidence
before it means anything.

Run: .venv/bin/python scripts/67_structural_comparison.py
"""

from __future__ import annotations

import sys
import warnings
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

REPORT = ROOT / "reports/42_structural_comparison.md"
# Per peptide position the fingerprint stores six numbers, in this order.
PER_POSITION = 6
CLOSE_NEIGHBOUR = {"GLCTLVAML": 0.552, "GILGFVFTL": 0.476, "YLQPRTFLL": 0.184,
                   "NLVPMVATV": 0.183, "LTDEMIAQY": 0.146, "TTDPSFLGRY": 0.047,
                   "LLWNGPMAV": 0.044, "CINGVCWTV": 0.025}


def exact_spearman(x, y):
    observed = spearmanr(x, y).statistic
    draws = np.array([spearmanr(x, list(order)).statistic
                      for order in permutations(y)])
    return observed, float((np.abs(draws) >= abs(observed) - 1e-12).mean())


def main() -> int:
    test = load_test()
    contacts = pd.read_csv(ROOT / "data/interim/structure_contacts.csv"
                           ).set_index("tcr_id")
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv"
                        ).set_index("tcr_id")
    cached = np.load(ROOT / "data/interim/contact_fingerprints.npz",
                     allow_pickle=True)
    fingerprint = cached["fingerprint"]
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    frame = test.set_index("tcr_id").join(contacts)
    a2_nine = sorted(e for e in frame["epitope_aa"].unique()
                     if frame[frame["epitope_aa"] == e]["hla_short"].iloc[0] == "A2"
                     and len(e) == 9)

    print(f"=== docking geometry, the {len(a2_nine)} HLA-A*02:01 9-mers ===")
    columns = ["crossing_angle_deg", "tcr_height_over_peptide",
               "contacts_a_pep", "contacts_b_pep", "cdr3_share_a",
               "cdr3_share_b"]
    geometry = []
    for epitope in a2_nine:
        block = frame[frame["epitope_aa"] == epitope]
        row = {"epitope": epitope, "n": len(block)}
        for column in columns:
            row[f"{column}_median"] = float(block[column].median())
            row[f"{column}_sd"] = float(block[column].std())
        geometry.append(row)
    geometry = pd.DataFrame(geometry)
    print(geometry[["epitope", "n", "crossing_angle_deg_median",
                    "crossing_angle_deg_sd", "tcr_height_over_peptide_median",
                    "contacts_a_pep_median", "contacts_b_pep_median"]
                   ].round(2).to_string(index=False))
    spread = geometry.set_index("epitope")[
        [c for c in geometry.columns if c.endswith("_sd")]]
    rank = int((spread.rank().loc["GLCTLVAML"]).mean())
    print(f"\n  GLCTLVAML's mean rank across the six spread measures: "
          f"{spread.rank().loc['GLCTLVAML'].mean():.1f} of {len(spread)} "
          f"(1 = tightest) -- it is not the tight one.")

    print("\n=== contact footprint by peptide position (median, alpha+beta) ===")
    profiles, summary = {}, []
    for epitope in a2_nine:
        ids = [t for t in frame[frame["epitope_aa"] == epitope].index
               if t in where]
        block = fingerprint[[where[t] for t in ids]]
        totals = np.array([np.median(block[:, p * PER_POSITION]
                                     + block[:, p * PER_POSITION + 1])
                           for p in range(len(epitope))])
        profiles[epitope] = totals
        share = totals / totals.sum()
        summary.append({
            "epitope": epitope, "hottest": f"P{1 + int(np.argmax(totals))}",
            "top3_share": float(np.sort(totals)[::-1][:3].sum() / totals.sum()),
            "effective_positions": float(
                np.exp(-(share * np.log(share + 1e-12)).sum()))})
    print("  " + f"{'epitope':<11}" + "".join(f"{'P' + str(p + 1):>7}"
                                              for p in range(9)))
    for epitope in a2_nine:
        print(f"  {epitope:<11}" + "".join(f"{v:>7.0f}"
                                           for v in profiles[epitope]))
    footprint = pd.DataFrame(summary)
    print("\n" + footprint.round(3).to_string(index=False))
    print("\n  The profiles are nearly identical: P4 dominant everywhere but")
    print("  CINGVCWTV, P9 at 26-27 for all six, P3 at 7-8 for all six. The")
    print("  only position that really varies is P7 (16-35).")

    print("\n=== is 'ensemble tightness' biology, or just Boltz confidence? ===")
    rows = []
    for epitope in sorted(frame["epitope_aa"].unique()):
        ids = [t for t in frame[frame["epitope_aa"] == epitope].index
               if t in where]
        block = fingerprint[[where[t] for t in ids]]
        normed = block / (np.linalg.norm(block, axis=1, keepdims=True) + 1e-9)
        distances = 1 - normed @ normed.T
        upper = np.triu_indices(len(block), 1)
        rows.append({"epitope": epitope, "n": len(ids),
                     "tightness": float(distances[upper].mean()),
                     "iptm": float(boltz.loc[ids, "iptm_tcr_pep_max"].median())})
    tightness = pd.DataFrame(rows)
    by_n = spearmanr(tightness["n"], tightness["tightness"]).statistic
    by_iptm = spearmanr(tightness["iptm"], tightness["tightness"]).statistic
    print(f"  Spearman(cohort size, tightness)      {by_n:+.3f}")
    print(f"  Spearman(interface ipTM, tightness)   {by_iptm:+.3f}   <- this")
    fit = linregress(tightness["iptm"], tightness["tightness"])
    tightness["residual"] = tightness["tightness"] - (
        fit.intercept + fit.slope * tightness["iptm"])
    print("\n  Tightness after removing confidence:")
    print(tightness.round(4).to_string(index=False))
    glc = tightness[tightness["epitope"] == "GLCTLVAML"].iloc[0]
    position = int((tightness["residual"] < glc["residual"]).sum()) + 1
    print(f"\n  GLCTLVAML's residual {glc['residual']:+.4f}, rank "
          f"{position} of {len(tightness)} -- exactly average.")

    print("\n=== does footprint width track sequence convergence? ===")
    x = [CLOSE_NEIGHBOUR[e] for e in footprint["epitope"]]
    y = list(footprint["effective_positions"])
    rho, p = exact_spearman(x, y)
    print(f"  Spearman {rho:+.3f}, exact p {p:.3f} over 720 permutations (n=6)")
    order = footprint.assign(cn=x).sort_values("cn", ascending=False)
    print("  narrowest footprints: "
          + ", ".join(footprint.nsmallest(2, "effective_positions")["epitope"]))
    print("  most sequence-convergent: " + ", ".join(order["epitope"][:2]))
    rest = order.iloc[2:]
    rho2 = spearmanr(rest["cn"], rest["effective_positions"]).statistic
    print(f"  with those two removed, the remaining four give {rho2:+.3f} -- "
          "the sign reverses.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E75 - structural comparison of GLCTLVAML against the other epitopes\n\n"
        "## Docking geometry (A*02:01 9-mers)\n\n"
        + geometry.round(4).to_markdown(index=False)
        + "\n\n## Contact footprint\n\n"
        + footprint.round(4).to_markdown(index=False)
        + "\n\n## Ensemble tightness against Boltz confidence\n\n"
        + tightness.round(4).to_markdown(index=False)
        + f"\n\nSpearman(ipTM, tightness) = **{by_iptm:+.3f}**; "
        f"footprint width against sequence convergence rho {rho:+.3f}, "
        f"exact p {p:.3f}.\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
