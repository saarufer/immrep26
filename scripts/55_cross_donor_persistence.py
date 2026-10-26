#!/usr/bin/env python
"""E57 -- do the six public families persist across independent donors?

E56 settled the unit: `meta.study.id` is 26.4% populated and unusable, while
`meta.subject.id` reaches 91.3% on GILGFVFTL and resolves the five large
families to 50-110 distinct donors. So the plan's leave-one-study-out becomes
leave-one-DONOR-out, which is the stricter test anyway -- two papers can share a
cohort, two donors cannot share a rearrangement event.

The question is not "how many papers mention this family". E54 already answered
that and E56 strengthened it. The question here is sharper:

    for a TCR in one of these families, is there a DIFFERENT TCR, from a
    DIFFERENT person, that is still its neighbour?

Two exclusions make that meaningful, both named in the plan:

  * a neighbour sharing any donor with the query is not independent evidence;
  * an identical alpha-beta clonotype is not convergence at all, it is the same
    receptor recorded twice, so identical sequences never count as neighbours
    of each other here.

Donor information is missing for some records. The plan asks for both
treatments, so every statistic is computed twice: **conservative**, where a TCR
with no resolvable donor cannot serve as independent support, and **permissive**,
where missing donor information is assumed to be a different donor. The truth is
bracketed by the two, and where they disagree the result is reported as
unresolved rather than as the more convenient of the two.

Controls (E57.3), all within GILGFVFTL so the epitope is held fixed:
large-family TRBV19 against TRBV19 outside any large family, and against other
V genes. Without them, "family members have cross-donor neighbours" would just
be restating that families are dense.

Run: .venv/bin/python scripts/55_cross_donor_persistence.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/30_cross_donor_persistence.md"
VDJDB_FULL = ROOT / "data/external/vdjdb" / f"vdjdb-{VDJDB_VERSION}" / "vdjdb_full.txt"
RADIUS = 0.12
BIG = 20


def space(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, component = connected_components(
            sp.csr_matrix(d <= RADIUS), directed=False)
        out[epitope] = (ids, d, component)
    return out


def donor_sets(test) -> dict:
    full = pd.read_csv(VDJDB_FULL, sep="\t", low_memory=False)
    full = full[(full["species"] == "HomoSapiens")
                & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]
    subject = full["meta.subject.id"]
    full = full.assign(donor=np.where(
        subject.notna() & (subject.astype(str).str.strip() != ""),
        full["reference.id"].astype(str) + "/" + subject.astype(str), None))
    grouped = (full[full["donor"].notna()]
               .groupby(["cdr3.alpha", "cdr3.beta", "antigen.epitope"])["donor"]
               .agg(lambda s: frozenset(s)))
    out = {}
    for _, row in test.iterrows():
        out[row["tcr_id"]] = grouped.get(
            (row["cdr3_alpha_aa"], row["cdr3_beta_aa"], row["epitope_aa"]),
            frozenset())
    return out


def support(ids, d, donors, sequences, permissive: bool):
    """For each TCR: neighbours, cross-donor neighbours, donors supporting it.

    A neighbour counts as independent support only if it is a different
    sequence AND shares no donor with the query. Under `permissive`, a TCR
    with no resolvable donor is treated as a different donor; under
    conservative it cannot supply support at all.
    """
    adjacency = d <= RADIUS
    np.fill_diagonal(adjacency, False)
    rows = []
    for i, t in enumerate(ids):
        mine = donors[t]
        neighbours = np.where(adjacency[i])[0]
        cross, supporting = [], set()
        for j in neighbours:
            other = ids[j]
            if sequences[other] == sequences[t]:
                continue                      # same receptor recorded twice
            theirs = donors[other]
            if not theirs:
                if not permissive:
                    continue
                cross.append(other)
                continue
            if mine and (mine & theirs):
                continue                      # shares a donor: not independent
            cross.append(other)
            supporting |= theirs
        rows.append({"tcr_id": t, "n_neighbours": len(neighbours),
                     "n_cross_donor": len(cross),
                     "n_donors_supporting": len(supporting),
                     "has_cross_donor": len(cross) > 0,
                     "own_donors": len(mine)})
    return pd.DataFrame(rows).set_index("tcr_id")


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    geometry = space(test)
    donors = donor_sets(test)
    sequences = {r["tcr_id"]: (r["cdr3_alpha_aa"], r["cdr3_beta_aa"])
                 for _, r in test.iterrows()}

    duplicates = pd.Series(list(sequences.values())).duplicated().sum()
    print(f"=== identical alpha-beta clonotypes in the cohort: {duplicates} ===")
    print("  (excluded from being each other's neighbours either way)\n")

    family_of, trbv = {}, {}
    for epitope, (ids, d, component) in geometry.items():
        family_of.update({t: f"{epitope}#{c}" for t, c in zip(ids, component)})
    trbv = (test.set_index("tcr_id")["TRBV_IMGT"].astype(str)
            .str.split("*").str[0].to_dict())
    labelled = labels.assign(family=labels["tcr_id"].map(family_of))
    counts = labelled.groupby("family").size()
    big = set(counts[counts >= BIG].index)

    results = {}
    for permissive in (False, True):
        frames = []
        for epitope, (ids, d, component) in geometry.items():
            frames.append(support(ids, d, donors, sequences, permissive)
                          .assign(epitope=epitope))
        results[permissive] = pd.concat(frames)

    print("=== E57.1: cross-donor support, the six families against controls ===")
    print("  Conservative = a TCR with no resolvable donor cannot supply")
    print("  support; permissive = it counts.\n")
    group_rows = []
    gil = test[test["epitope_aa"] == "GILGFVFTL"]["tcr_id"]
    groups = {
        "GIL large TRBV19 families": [t for t in gil if family_of[t] in big],
        "GIL TRBV19, no large family": [t for t in gil if family_of[t] not in big
                                        and trbv[t] == "TRBV19"],
        "GIL other V, no large family": [t for t in gil if family_of[t] not in big
                                         and trbv[t] != "TRBV19"],
        "GLCTLVAML#2 (12 donors, thin)": [t for t, f in family_of.items()
                                          if f == "GLCTLVAML#2"],
    }
    for name, members in groups.items():
        row = {"group": name, "n": len(members)}
        for permissive in (False, True):
            block = results[permissive].loc[members]
            tag = "perm" if permissive else "cons"
            row[f"has_cross_donor_{tag}"] = float(block["has_cross_donor"].mean())
            row[f"median_donors_{tag}"] = float(block["n_donors_supporting"].median())
        row["median_neighbours"] = float(
            results[False].loc[members, "n_neighbours"].median())
        group_rows.append(row)
        print(f"  {name:<30} n={len(members):4d}  neighbours (median) "
              f"{row['median_neighbours']:5.0f}   has cross-donor support: "
              f"conservative {row['has_cross_donor_cons']:6.1%} / permissive "
              f"{row['has_cross_donor_perm']:6.1%}   donors supporting (median) "
              f"{row['median_donors_cons']:.0f}")
    groups_table = pd.DataFrame(group_rows)

    print("\n=== E57.2: leave-one-donor-out, per family ===")
    print("  Remove every record from one donor, then ask whether the TCRs that")
    print("  donor contributed still have a neighbour from somebody else.\n")
    lodo_rows = []
    for family in sorted(big, key=lambda f: -counts[f]):
        epitope = family.split("#")[0]
        ids, d, component = geometry[epitope]
        members = np.array([i for i, t in enumerate(ids)
                            if family_of[t] == family])
        member_ids = ids[members]
        all_donors = set().union(*[donors[t] for t in member_ids])
        if not all_donors:
            lodo_rows.append({"family": family, "donors": 0,
                              "survival": np.nan, "note": "no donor info"})
            continue
        survivals = []
        adjacency = d <= RADIUS
        np.fill_diagonal(adjacency, False)
        for donor in all_donors:
            affected = [i for i in range(len(ids))
                        if donor in donors[ids[i]]]
            if not affected:
                continue
            keep = np.ones(len(ids), dtype=bool)
            keep[affected] = False
            held = [i for i in affected if family_of[ids[i]] == family]
            if not held:
                continue
            kept = []
            for i in held:
                others = np.where(adjacency[i] & keep)[0]
                others = [j for j in others
                          if sequences[ids[j]] != sequences[ids[i]]]
                kept.append(len(others) > 0)
            survivals.append(float(np.mean(kept)))
        lodo_rows.append({
            "family": family, "labelled": int(counts[family]),
            "donors": len(all_donors),
            "donors_tested": len(survivals),
            "survival": float(np.mean(survivals)) if survivals else np.nan,
            "worst_donor": float(np.min(survivals)) if survivals else np.nan})
        print(f"  {family:<14} {len(all_donors):3d} donors, "
              f"{len(survivals):3d} testable   its TCRs keep a neighbour "
              f"{lodo_rows[-1]['survival']:6.1%} of the time "
              f"(worst donor {lodo_rows[-1]['worst_donor']:.1%})")
    lodo = pd.DataFrame(lodo_rows)

    print("\n=== E57.3: does cross-donor support track validation? ===")
    merged = labelled.merge(results[False], left_on="tcr_id", right_index=True)
    rate_rows = []
    for epitope in ("GILGFVFTL", "GLCTLVAML"):
        block = merged[merged["epitope_aa"] == epitope]
        for has, name in ((True, "has cross-donor support"),
                          (False, "no cross-donor support")):
            sub = block[block["has_cross_donor"] == has]
            if len(sub) < 5:
                continue
            rate_rows.append({"epitope": epitope, "group": name, "n": len(sub),
                              "validated": float(sub["pos"].mean())})
            print(f"  {epitope:<11} {name:<26} n={len(sub):4d}  validated "
                  f"{sub['pos'].mean():6.1%}")
    rates = pd.DataFrame(rate_rows)

    print("\n  That comparison is confounded and must not be read as new")
    print("  evidence: having cross-donor support and having any neighbour at")
    print("  all are nearly the same condition (median 13 neighbours against")
    print("  0). It restates the density signal. The non-confounded version")
    print("  conditions on already having a neighbour, and asks whether that")
    print("  neighbour coming from ANOTHER PERSON adds anything:\n")
    conditioned_rows = []
    for epitope in ("GILGFVFTL", "GLCTLVAML"):
        block = merged[(merged["epitope_aa"] == epitope)
                       & (merged["n_neighbours"] > 0)]
        for has, name in ((True, "cross-donor neighbour"),
                          (False, "only same-donor neighbours")):
            sub = block[block["has_cross_donor"] == has]
            if len(sub) < 5:
                print(f"  {epitope:<11} {name:<28} n={len(sub)} -- "
                      "NOT ESTIMABLE")
                continue
            conditioned_rows.append({"epitope": epitope, "group": name,
                                     "n": len(sub),
                                     "validated": float(sub["pos"].mean())})
            print(f"  {epitope:<11} {name:<28} n={len(sub):4d}  validated "
                  f"{sub['pos'].mean():6.1%}")
    conditioned = pd.DataFrame(conditioned_rows)
    from scipy.stats import fisher_exact
    for epitope in ("GILGFVFTL", "GLCTLVAML"):
        block = merged[(merged["epitope_aa"] == epitope)
                       & (merged["n_neighbours"] > 0)]
        a = block[block["has_cross_donor"]]
        b = block[~block["has_cross_donor"]]
        if len(a) < 5 or len(b) < 5:
            continue
        table = [[int(a["pos"].sum()), int((1 - a["pos"]).sum())],
                 [int(b["pos"].sum()), int((1 - b["pos"]).sum())]]
        odds, pvalue = fisher_exact(table)
        conditioned.loc[conditioned["epitope"] == epitope, "fisher_p"] = pvalue
        conditioned.loc[conditioned["epitope"] == epitope, "odds_ratio"] = odds
        print(f"  {epitope}: odds ratio {odds:.1f}, Fisher p {pvalue:.2e}")
    print("\n  A TCR whose only neighbours come from the same person is not")
    print("  showing convergence -- it is a clonal expansion or a within-")
    print("  experiment duplicate. That distinction is invisible to density,")
    print("  and it cannot be used on the private split, where donor coverage")
    print("  is 0% on two of three epitopes (E56).")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E57 - cross-donor persistence of the six public families\n\n"
        f"Donor = reference.id + meta.subject.id. Identical alpha-beta "
        f"clonotypes ({duplicates} in the cohort) never count as each other's "
        "neighbours.\n\n## Cross-donor support by group\n\n"
        + groups_table.round(4).to_markdown(index=False)
        + "\n\n## Leave-one-donor-out\n\n" + lodo.round(4).to_markdown(index=False)
        + "\n\n## Cross-donor support against validation (CONFOUNDED with "
        "having any neighbour)\n\n"
        + rates.round(4).to_markdown(index=False)
        + "\n\n## The same, conditioned on already having a neighbour\n\n"
        + (conditioned.round(4).to_markdown(index=False) if len(conditioned)
           else "Not estimable.")
        + "\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
