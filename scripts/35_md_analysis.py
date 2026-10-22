#!/usr/bin/env python
"""Analysis for plan card C11, and a descriptor audit of E37 on its own structures.

Two parts, because the second only becomes available once the trajectories are
written and the first is worth having regardless.

`--static` re-measures E31's quantity on the relaxed coordinates E37 produced,
using two descriptors rather than one: the 4.5 A contact count E37 reported, and
the minimum heavy-atom distance from the mutated site to the peptide. Both read
the same physical quantity. The count is a thresholded sum and so is dominated
by whichever sites happen to be in contact at all; the distance is continuous
and always defined. If they disagree, the count's significance in E37 rests on a
few large values, which is a caveat E37 does not currently carry.

`--md` compares the same sites across the trajectories and applies the test the
plan demands: the discordant-minus-control difference has to exceed the
replica-to-replica spread of the same quantity, or there is no result. It also
checks whether the interface drifts systematically, which would mean the
implicit-solvent setup is losing the pose rather than sampling it, and whether
the effect is carried by a single epitope.

Run: python scripts/35_md_analysis.py --static
     python scripts/35_md_analysis.py --md
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

RELAXED = ROOT / "data/interim/relaxed"
PAIRS = ROOT / "data/interim/md_pairs.json"
CONTACTS = ROOT / "data/interim/md_residue_contacts.csv"
TRACE = ROOT / "data/interim/md_interface_trace.csv"
REPORT = ROOT / "reports/12_md_interface_stability.md"

CONTACT_CUTOFF = 4.5


def heavy_atoms(path: Path):
    """(chain, residue number, xyz) for heavy atoms, from the relaxed file."""
    chains, residues, coords = [], [], []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if (line[76:78].strip() or line[12:16].strip()[:1]) == "H":
            continue
        chains.append(line[21])
        residues.append(int(line[22:26]))
        coords.append([float(line[30 + 8 * i:38 + 8 * i]) for i in range(3)])
    return np.array(chains), np.array(residues), np.asarray(coords)


def peptide_mask(chains, residues):
    """The peptide is the 8-11 residue chain; OpenMM renamed it on write."""
    for chain in np.unique(chains):
        if 8 <= len(np.unique(residues[chains == chain])) <= 11:
            return chains == chain
    return None


def load_pairs() -> pd.DataFrame:
    pairs = json.loads(PAIRS.read_text())
    rows = []
    for group, entries in pairs.items():
        for entry in entries:
            for role in ("a", "b"):
                sites = entry[f"{role}_site_residues"]
                for key, numbers in sites.items():
                    for number in numbers:
                        rows.append({
                            "group": group, "pair": entry["pair"],
                            "epitope": entry["epitope"], "role": role,
                            "tcr_id": entry[f"{role}_id"],
                            "pos": entry[f"{role}_pos"],
                            "chain": key, "resid": number,
                        })
    return pd.DataFrame(rows)


def compare(frame: pd.DataFrame, column: str, label: str) -> dict:
    disc = frame.loc[frame["group"] == "discordant", column].dropna()
    ctl = frame.loc[frame["group"] == "control", column].dropna()
    stat = mannwhitneyu(disc, ctl, alternative="two-sided")
    out = {"descriptor": label, "discordant": disc.mean(), "control": ctl.mean(),
           "n_disc": len(disc), "n_ctl": len(ctl), "p": stat.pvalue}
    print(f"  {label:<34} disc {disc.mean():7.3f}  ctl {ctl.mean():7.3f}  "
          f"p={stat.pvalue:.3f}  (n {len(disc)}/{len(ctl)})")
    return out


def static() -> pd.DataFrame:
    sites = load_pairs()
    print(f"{sites.pair.nunique()} pairs, {len(sites)} mutated-site measurements "
          f"on relaxed coordinates\n")

    records = []
    for tcr_id, group in sites.groupby("tcr_id"):
        path = RELAXED / f"{tcr_id}.pdb"
        if not path.exists():
            continue
        chains, residues, coords = heavy_atoms(path)
        peptide = peptide_mask(chains, residues)
        if peptide is None:
            continue
        tree = cKDTree(coords[peptide])
        for row in group.to_dict("records"):
            mask = (chains == row["chain"]) & (residues == row["resid"])
            if not mask.any():
                continue
            hits = tree.query_ball_point(coords[mask], CONTACT_CUTOFF)
            records.append(row | {
                "contacts": int(sum(len(h) for h in hits)),
                "min_dist": float(tree.query(coords[mask])[0].min()),
            })
    frame = pd.DataFrame(records)

    print("=== E37's quantity, two descriptors, same structures ===")
    rows = [compare(frame, "contacts", "peptide contacts at site (count)"),
            compare(frame, "min_dist", "site-to-peptide distance (A)")]

    # Not an outlier artefact: Mann-Whitney is already a rank test, so it is
    # invariant to any monotone transform of the descriptor. A first version of
    # this audit "checked" the outliers by rank-transforming and re-running,
    # which returns the identical p-value by construction and verifies nothing.
    # The descriptors disagree for a different and more informative reason.
    zero = (frame["contacts"] == 0).mean()
    print(f"\n  Both tests above are rank tests, so neither is outlier-driven."
          f"\n  They disagree because the count is {zero:.0%} ties at zero:")
    print("   contact-count quantiles  "
          + "  ".join(f"q{int(q * 100)}={frame['contacts'].quantile(q):.0f}"
                      for q in (0.5, 0.75, 0.9, 1.0)))
    print("  With three quarters of the sites tied, the count test is in effect a\n"
          "  comparison of how often a site touches the peptide at all, while the\n"
          "  distance uses every site. So the two answer different questions.")

    from scipy.stats import fisher_exact

    facing = frame.assign(faces=frame["contacts"] > 0)
    table = pd.crosstab(facing["group"], facing["faces"])
    odds, p_fisher = fisher_exact(table.reindex(
        index=["discordant", "control"], columns=[False, True]).to_numpy())
    disc_rate = facing.loc[facing["group"] == "discordant", "faces"].mean()
    ctl_rate = facing.loc[facing["group"] == "control", "faces"].mean()
    print(f"\n  site touches the peptide at all   disc {disc_rate:.3f}  "
          f"ctl {ctl_rate:.3f}  p={p_fisher:.3f}  (Fisher, odds {odds:.2f})")
    rows.append({"descriptor": "site touches peptide at all (binary)",
                 "discordant": disc_rate, "control": ctl_rate,
                 "n_disc": int(table.loc["discordant"].sum()),
                 "n_ctl": int(table.loc["control"].sum()), "p": p_fisher})

    print("\n  Among sites that do touch, is the discordant one closer?")
    rows.append(compare(frame[frame["contacts"] > 0], "min_dist",
                        "distance | site touches (A)"))
    rows.append(compare(frame[frame["contacts"] > 0], "contacts",
                        "contacts | site touches (count)"))
    print("\n  If the binary rate is significant and neither conditional test is,\n"
          "  the effect is whether the mutated site faces the peptide, not how\n"
          "  tightly it sits there -- and E37's contact-count and peptide-facing\n"
          "  rows are then one finding reported twice, not two.")
    return frame, pd.DataFrame(rows)


def md() -> None:
    if not CONTACTS.exists():
        raise SystemExit(f"{CONTACTS.relative_to(ROOT)} not written yet")
    sites = load_pairs()
    contacts = pd.read_csv(CONTACTS)
    trajectories = contacts.groupby(["tcr_id", "replica"]).ngroups
    print(f"{contacts.tcr_id.nunique()} complexes, {trajectories} trajectories, "
          f"{contacts.n_frames.iloc[0]} frames each\n")

    at_site = sites.merge(contacts, on=["tcr_id", "chain", "resid"], how="inner")
    if at_site.empty:
        raise SystemExit("no mutated site matched a simulated residue")

    # Per complex-site: mean over replicas, and the spread between them.
    keys = ["group", "pair", "epitope", "role", "pos", "tcr_id", "chain", "resid"]
    per_site = at_site.groupby(keys, as_index=False).agg(
        count_mean=("count_mean", "mean"), count_replica_sd=("count_mean", "std"),
        dist_mean=("dist_mean", "mean"), dist_replica_sd=("dist_mean", "std"),
        frac=("frac_in_contact", "mean"), frac_replica_sd=("frac_in_contact", "std"),
        replicas=("replica", "nunique"))

    print("=== discordant versus control at the mutated site, over MD ===")
    rows = [compare(per_site, "count_mean", "contacts, time-averaged"),
            compare(per_site, "dist_mean", "distance, time-averaged (A)"),
            compare(per_site, "frac", "fraction of frames in contact")]

    print("\n=== the test the plan requires: effect against replica spread ===")
    guard = []
    for column, spread, label in (("count_mean", "count_replica_sd", "contacts"),
                                  ("dist_mean", "dist_replica_sd", "distance"),
                                  ("frac", "frac_replica_sd", "contact fraction")):
        disc = per_site.loc[per_site["group"] == "discordant", column].mean()
        ctl = per_site.loc[per_site["group"] == "control", column].mean()
        effect = abs(disc - ctl)
        noise = per_site[spread].median()
        verdict = "exceeds spread" if effect > noise else "WITHIN SPREAD -- no result"
        guard.append({"quantity": label, "effect": effect,
                      "median_replica_sd": noise, "verdict": verdict})
        print(f"  {label:<18} |effect| {effect:7.3f}   median replica sd {noise:7.3f}"
              f"   {verdict}")

    print("\n=== is the pose held, or lost? interface contact drift ===")
    if TRACE.exists():
        trace = pd.read_csv(TRACE)
        slopes = []
        for (tcr_id, replica), group in trace.groupby(["tcr_id", "replica"]):
            if len(group) < 10:
                continue
            slope = np.polyfit(group["ps"], group["interface_contacts"], 1)[0]
            first = group.nsmallest(max(len(group) // 10, 1), "frame")
            last = group.nlargest(max(len(group) // 10, 1), "frame")
            slopes.append({"tcr_id": tcr_id, "replica": replica,
                           "slope_per_ns": slope * 1000,
                           "start": first["interface_contacts"].mean(),
                           "end": last["interface_contacts"].mean()})
        drift = pd.DataFrame(slopes)
        change = (drift["end"] - drift["start"]) / drift["start"]
        print(f"  median interface change over the trajectory {change.median():+.1%}")
        print(f"  trajectories losing more than 25% of contacts: "
              f"{(change < -0.25).sum()}/{len(drift)}")
        print(f"  median slope {drift['slope_per_ns'].median():+.1f} contacts/ns")

    print("\n=== is one epitope carrying it? ===")
    by_epitope = per_site.pivot_table(
        index="epitope", columns="group",
        values=["count_mean", "dist_mean"], aggfunc="mean")
    print(by_epitope.round(3).to_string())
    counts = per_site.groupby(["epitope", "group"]).size().unstack(fill_value=0)
    print("\n  sites per epitope:")
    print(counts.to_string())

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# C11: is the interface difference dynamically stable?\n\n"
        f"{contacts.tcr_id.nunique()} complexes, {trajectories} trajectories.\n\n"
        "## Group comparison at the mutated site\n\n"
        + pd.DataFrame(rows).round(4).to_markdown(index=False)
        + "\n\n## Effect against replica-to-replica spread\n\n"
        + pd.DataFrame(guard).round(4).to_markdown(index=False)
        + "\n\n## By epitope\n\n" + by_epitope.round(3).to_markdown() + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--static", action="store_true")
    parser.add_argument("--md", action="store_true")
    args = parser.parse_args()
    if args.static:
        static()
    if args.md:
        md()
    if not (args.static or args.md):
        parser.error("choose --static or --md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
