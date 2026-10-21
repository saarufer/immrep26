#!/usr/bin/env python
"""Are near-positive discordant negatives fine specificity, or assay dropout?

Plan cards S03, C09, C03, C04 on one pair set. E20 found 34 negatives that are
single-residue variants of validated positives, and E27 showed their damage
sits in the cross-stratum comparisons rather than inside their own cluster. What
nothing has settled is *why* they differ, and the two candidates imply opposite
things: a mutation that abolishes a peptide contact is real fine specificity and
no sequence method can be blamed for missing it, whereas a mutation at a
solvent-facing position points at the assay or the annotation.

The falsifiable test needs a null. A discordant pair's mutation is compared
against **matched same-label pairs** -- positive/positive neighbours at a
comparable sequence distance within the same epitope -- so "mutations tend to be
peptide-facing" cannot be mistaken for a property of discordance. Without that
control the analysis would only describe where mutations generally fall.

Measured per mutated position: whether it contacts the peptide in the predicted
structure (C09), the change in peptide-contact count between the two members
(C03), and the change in charge, polarity and volume (C04).

Run: python scripts/26_discordant_pair_anatomy.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact, mannwhitneyu

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

PREDICTIONS = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
CONTACT_CUTOFF = 4.5
OUT = ROOT / "reports/08_discordant_pair_anatomy.md"

CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.5}
HYDROPHOBIC = set("AVLIMFWCY")
# Zamyatnin residue volumes, cubic angstroms.
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}


def differences(a: str, b: str) -> list[int]:
    """Positions where two equal-length strings differ."""
    return [i for i, (x, y) in enumerate(zip(a, b))] if len(a) != len(b) else [
        i for i, (x, y) in enumerate(zip(a, b)) if x != y
    ]


def find_pairs(frame: pd.DataFrame, want_labels: tuple[int, int], max_diff: int = 2):
    """Pairs within an epitope sharing one CDR3 exactly, differing by <=max_diff
    on the other, with the requested label combination."""
    pairs = []
    for epitope, cohort in frame.groupby("epitope_aa"):
        rows = cohort.to_dict("records")
        for i, first in enumerate(rows):
            for second in rows[i + 1:]:
                labels = (first["pos"], second["pos"])
                if sorted(labels) != sorted(want_labels):
                    continue
                for shared, varied in (("alpha", "beta"), ("beta", "alpha")):
                    same = first[f"cdr3_{shared}_aa"] == second[f"cdr3_{shared}_aa"]
                    x, y = first[f"cdr3_{varied}_aa"], second[f"cdr3_{varied}_aa"]
                    if not same or len(x) != len(y):
                        continue
                    sites = differences(x, y)
                    if 1 <= len(sites) <= max_diff:
                        # Orient so `a` is the positive member when labels differ.
                        a, b = (first, second) if first["pos"] >= second["pos"] else (second, first)
                        pairs.append(
                            {
                                "epitope": epitope, "varied_chain": varied,
                                "a_id": a["tcr_id"], "b_id": b["tcr_id"],
                                "a_pos": a["pos"], "b_pos": b["pos"],
                                "sites": differences(
                                    a[f"cdr3_{varied}_aa"], b[f"cdr3_{varied}_aa"]
                                ),
                                "a_seq": a[f"cdr3_{varied}_aa"],
                                "b_seq": b[f"cdr3_{varied}_aa"],
                            }
                        )
                        break
    return pd.DataFrame(pairs)


def peptide_contact_profile(tcr_id: str, chain: str, test: pd.DataFrame):
    """Per-residue peptide-contact counts for one chain, keyed by offset into CDR3."""
    pdb = PREDICTIONS / tcr_id / f"{tcr_id}_model_0.pdb"
    if not pdb.exists():
        return None
    coords: dict[str, list] = {}
    resid: dict[str, list] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip() or line[12:16].strip()[:1]
        if element == "H":
            continue
        c = line[21]
        coords.setdefault(c, []).append(
            (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        )
        resid.setdefault(c, []).append(int(line[22:26]))
    key = "A" if chain == "alpha" else "B"
    if key not in coords or "P" not in coords:
        return None
    chain_xyz = np.asarray(coords[key]); chain_res = np.asarray(resid[key])
    pep_tree = cKDTree(np.asarray(coords["P"]))

    row = test.loc[tcr_id]
    full = row.full_seq_alpha_aa if chain == "alpha" else row.full_seq_beta_aa
    loop = row.cdr3_alpha_aa if chain == "alpha" else row.cdr3_beta_aa
    start = full.index(loop) + int(chain_res.min())

    counts = {}
    for offset in range(len(loop)):
        atoms = chain_xyz[chain_res == start + offset]
        counts[offset] = (
            int(sum(len(x) for x in pep_tree.query_ball_point(atoms, CONTACT_CUTOFF)))
            if len(atoms) else 0
        )
    return counts


def annotate(pairs: pd.DataFrame, test: pd.DataFrame) -> pd.DataFrame:
    records = []
    for pair in pairs.to_dict("records"):
        a_profile = peptide_contact_profile(pair["a_id"], pair["varied_chain"], test)
        b_profile = peptide_contact_profile(pair["b_id"], pair["varied_chain"], test)
        if a_profile is None or b_profile is None:
            continue
        for site in pair["sites"]:
            x, y = pair["a_seq"][site], pair["b_seq"][site]
            records.append(
                {
                    "epitope": pair["epitope"],
                    "discordant": pair["a_pos"] != pair["b_pos"],
                    "chain": pair["varied_chain"],
                    "a_res": x, "b_res": y,
                    # C09: does the mutated position touch the peptide at all
                    "peptide_facing": bool(a_profile.get(site, 0) + b_profile.get(site, 0) > 0),
                    "contacts_at_site": a_profile.get(site, 0) + b_profile.get(site, 0),
                    # C03: change in peptide contacts at the mutated position
                    "delta_contacts": abs(a_profile.get(site, 0) - b_profile.get(site, 0)),
                    # C04: local chemical change
                    "delta_charge": abs(CHARGE.get(x, 0.0) - CHARGE.get(y, 0.0)),
                    "polarity_swap": (x in HYDROPHOBIC) != (y in HYDROPHOBIC),
                    "delta_volume": abs(VOLUME.get(x, 140.0) - VOLUME.get(y, 140.0)),
                    "involves_pro_gly": bool({x, y} & {"P", "G"}),
                }
            )
    return pd.DataFrame(records)


def main() -> int:
    test = load_test().set_index("tcr_id")
    labels = build_training_frame()

    discordant = find_pairs(labels, (1, 0))
    control = find_pairs(labels, (1, 1))
    print(f"discordant positive/negative pairs: {len(discordant)}")
    print(f"matched positive/positive control pairs: {len(control)}\n")
    if discordant.empty:
        sys.exit("no discordant pairs found")

    annotated = pd.concat(
        [annotate(discordant, test), annotate(control, test)], ignore_index=True
    )
    pd.set_option("display.width", 250)

    print("=== mutated sites, discordant pairs vs same-label controls ===")
    summary = annotated.groupby("discordant").agg(
        sites=("peptide_facing", "size"),
        peptide_facing=("peptide_facing", "mean"),
        mean_contacts_at_site=("contacts_at_site", "mean"),
        mean_delta_contacts=("delta_contacts", "mean"),
        mean_delta_charge=("delta_charge", "mean"),
        polarity_swap=("polarity_swap", "mean"),
        mean_delta_volume=("delta_volume", "mean"),
        pro_gly=("involves_pro_gly", "mean"),
    )
    print(summary.round(3).to_string())

    print("\n=== the test that matters: is the discordant mutation more often peptide-facing? ===")
    table = pd.crosstab(annotated["discordant"], annotated["peptide_facing"])
    print(table.to_string())
    if table.shape == (2, 2):
        odds, p = fisher_exact(table.to_numpy())
        print(f"  Fisher exact: odds ratio {odds:.3f}, p = {p:.3f}")
    for column in ("contacts_at_site", "delta_contacts", "delta_charge", "delta_volume"):
        groups = [g[column].to_numpy() for _, g in annotated.groupby("discordant")]
        if min(len(g) for g in groups) >= 3:
            _, p = mannwhitneyu(*groups)
            print(f"  {column:<20} discordant {groups[1].mean():.2f} vs control "
                  f"{groups[0].mean():.2f}   Mann-Whitney p = {p:.3f}")

    print("\n=== per-epitope counts (small-sample warning) ===")
    print(annotated.groupby(["epitope", "discordant"]).size().unstack(fill_value=0).to_string())

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Anatomy of near-positive discordant clonotypes\n\n"
        "Single-residue CDR3 variants with opposite functional labels, compared\n"
        "against same-label (positive/positive) neighbour pairs at comparable\n"
        "distance so that general properties of CDR3 mutation are not mistaken for\n"
        "properties of discordance.\n\n"
        f"- discordant pairs: {len(discordant)}\n- control pairs: {len(control)}\n\n"
        "## Mutated-site properties\n\n" + summary.round(3).to_markdown()
        + "\n\n## Peptide-facing contingency\n\n" + table.to_markdown() + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
