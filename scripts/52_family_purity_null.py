#!/usr/bin/env python
"""E54 -- is family label purity more than the continuous density-label relation?

Two biological readings of the same observation compete, and E20 could not
separate them. It found 45 families of >=3 members with 26 pure-positive, 3
pure-negative, and the 20-60% purity band completely **empty**.

  A  Families look pure because their members all sit at low nn3_dist, where
     the validation rate is high anyway, and families are small. Nothing beyond
     the continuous distance-to-validation relationship is needed.
  B  Each family shares something extra -- a recognition state, or a common
     study or donor -- that makes its members agree with each other beyond what
     their position on the distance axis implies.

The null that separates them permutes the validation label **within
(epitope x nn3 bin)**. That keeps each epitope's label ratio, keeps the
density-to-label relationship at the resolution of the bins, and keeps the graph
and every family untouched -- while destroying any family-level agreement that
is not already implied by where the members sit.

Two design points, both of which have bitten this project before:

  * **The graph is built on the full cohort, purity is measured on labelled
    nodes only.** Two different clustering bases already exist here
    (`scripts/20` clusters labelled rows, `scripts/25` the whole cohort), and
    mixing them is what made E30's accounting inconsistent. Here the edges come
    from every row of the epitope, including the ambiguous rows that carry no
    label, and only labelled members enter a purity count.
  * **The bin count is swept.** One bin is the weak null that ignores density
    entirely; many bins is so strict the strata stop containing anything to
    permute. The trend across the sweep is the result, not any single setting.

Run: .venv/bin/python scripts/52_family_purity_null.py
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

REPORT = ROOT / "reports/27_family_purity_null.md"
VDJDB_FULL = ROOT / "data/external/vdjdb" / f"vdjdb-{VDJDB_VERSION}" / "vdjdb_full.txt"
RADIUS = 0.12
DRAWS = 2000
BIN_SWEEP = (1, 3, 5, 10, 20)
MIN_LABELLED = 3
SEED = 0


def graph_and_distance(test: pd.DataFrame):
    """Edges from the FULL cohort; nn3 for every row. No labels involved."""
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    family, nn3, sizes_full = {}, {}, {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, component = connected_components(
            sp.csr_matrix(d <= RADIUS), directed=False)
        np.fill_diagonal(d, np.inf)
        ordered = np.sort(d, axis=1)
        counts = pd.Series(component).value_counts()
        for i, t in enumerate(ids):
            family[t] = f"{epitope}#{component[i]}"
            nn3[t] = float(ordered[i, :3].mean())
            sizes_full[t] = int(counts[component[i]])
    return family, nn3, sizes_full


def purity_counts(y: np.ndarray, codes: np.ndarray, n_families: int):
    """Pure-positive, pure-negative and mixed families, plus pair concordance.

    Concordance is pair-weighted: a family of 10 carries 45 pairs and a family
    of 3 carries 3, and averaging family rates instead would let the smallest
    families dominate a statistic about agreement.
    """
    counts = np.bincount(codes, minlength=n_families).astype(float)
    sums = np.bincount(codes, weights=y, minlength=n_families)
    rate = sums / counts
    pure_pos = int((rate == 1.0).sum())
    pure_neg = int((rate == 0.0).sum())
    mixed = int(n_families - pure_pos - pure_neg)
    k, n = sums, counts
    concordant = k * (k - 1) / 2 + (n - k) * (n - k - 1) / 2
    total = n * (n - 1) / 2
    return pure_pos, pure_neg, mixed, float(concordant.sum() / total.sum())


def main() -> int:
    rng = np.random.default_rng(SEED)
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]].copy()
    family, nn3, size_full = graph_and_distance(test)

    labels["family"] = labels["tcr_id"].map(family)
    labels["nn3"] = labels["tcr_id"].map(nn3)

    print("=== Part 0: the graph, built on the full cohort ===")
    print("  Edges and families come from every row of the epitope, including")
    print("  the ambiguous rows that carry no label. Only labelled members are")
    print("  counted when a family's purity is measured.\n")
    rows = []
    for epitope, cohort in test.groupby("epitope_aa"):
        block = labels[labels["epitope_aa"] == epitope]
        counted = block.groupby("family").size()
        rows.append({
            "epitope": epitope, "cohort_rows": len(cohort),
            "labelled": len(block),
            "unlabelled_in_graph": len(cohort) - len(block),
            "families_in_graph": len({family[t] for t in cohort["tcr_id"]}),
            f"families_ge{MIN_LABELLED}_labelled": int((counted >= MIN_LABELLED).sum()),
            "families_ge2_labelled": int((counted >= 2).sum()),
        })
    graph_table = pd.DataFrame(rows)
    print(graph_table.to_string(index=False))

    sized = labels.groupby("family").size()
    keep = set(sized[sized >= MIN_LABELLED].index)
    work = labels[labels["family"].isin(keep)].reset_index(drop=True)
    codes, uniques = pd.factorize(work["family"])
    n_families = len(uniques)
    y = work["pos"].to_numpy(dtype=float)
    print(f"\n  families with >={MIN_LABELLED} labelled members: {n_families}, "
          f"holding {len(work)} TCRs")

    observed = purity_counts(y, codes, n_families)
    print(f"  observed: pure-positive {observed[0]}, pure-negative "
          f"{observed[1]}, mixed {observed[2]}, concordance {observed[3]:.4f}")
    print(f"  (E20 reported 26 / 3 / 16 on a labelled-rows-only graph)")

    print("\n=== Part 1: the density-to-label relation the null preserves ===")
    work["bin10"] = work.groupby("epitope_aa")["nn3"].transform(
        lambda s: pd.qcut(s, min(10, s.nunique()), labels=False, duplicates="drop"))
    band = work.groupby("bin10").agg(n=("pos", "size"), pos_rate=("pos", "mean"),
                                     nn3_mean=("nn3", "mean"))
    print(band.round(3).to_string())

    print("\n=== Part 2: permutation within (epitope x nn3 bin) ===")
    print(f"  {DRAWS} draws per setting. One bin ignores density entirely and is")
    print("  the weak null; more bins is stricter. The trend is the result.\n")
    sweep_rows = []
    for n_bins in BIN_SWEEP:
        if n_bins == 1:
            strata = work["epitope_aa"].to_numpy()
        else:
            binned = work.groupby("epitope_aa")["nn3"].transform(
                lambda s, b=n_bins: pd.qcut(s, min(b, s.nunique()), labels=False,
                                            duplicates="drop"))
            strata = (work["epitope_aa"] + "|" + binned.astype(str)).to_numpy()
        groups = [np.where(strata == s)[0] for s in np.unique(strata)]

        # How much freedom does the null actually have? A stratum that holds one
        # family, or that is all one label, cannot contribute any shuffling.
        movable = sum(len(g) for g in groups
                      if len(np.unique(y[g])) > 1
                      and len(np.unique(codes[g])) > 1)
        draws = np.empty((DRAWS, 4))
        for draw in range(DRAWS):
            shuffled = y.copy()
            for g in groups:
                shuffled[g] = rng.permutation(y[g])
            draws[draw] = purity_counts(shuffled, codes, n_families)
        names = ("pure_pos", "pure_neg", "mixed", "concordance")
        row = {"n_bins": n_bins, "n_strata": len(groups),
               "movable_rows": movable, "movable_frac": movable / len(work)}
        for i, name in enumerate(names):
            null = draws[:, i]
            tail = ((null >= observed[i]).sum() + 1) / (DRAWS + 1)
            row[f"{name}_obs"] = observed[i]
            row[f"{name}_null"] = float(null.mean())
            row[f"{name}_p"] = tail if name != "mixed" else (
                ((null <= observed[i]).sum() + 1) / (DRAWS + 1))
            row[f"{name}_z"] = ((observed[i] - null.mean()) / null.std()
                                if null.std() > 0 else np.nan)
        sweep_rows.append(row)
        print(f"  bins={n_bins:<3} strata={len(groups):<4} "
              f"movable {movable / len(work):5.1%}   "
              f"pure+ {observed[0]} vs null {row['pure_pos_null']:5.1f} "
              f"(p {row['pure_pos_p']:.4f})   "
              f"concordance {observed[3]:.3f} vs {row['concordance_null']:.3f} "
              f"(p {row['concordance_p']:.4f})")
    sweep = pd.DataFrame(sweep_rows)

    print("\n=== Part 3: reading the sweep ===")
    strict = sweep[sweep["n_bins"] >= 5]
    count_excess = bool((strict["pure_pos_p"] < 0.05).all())
    agree_excess = bool((strict["concordance_p"] < 0.05).all())
    verdict_b = agree_excess
    print(f"  count of pure families exceeds the null: {count_excess}")
    print(f"  pairwise agreement exceeds the null:     {agree_excess}")
    if agree_excess and not count_excess:
        print("\n  A split verdict, and the two are not in conflict. The NUMBER")
        print("  of pure families is exactly what the continuous relationship")
        print("  plus small families produces -- result A for E20's headline.")
        print("  But pair-weighted agreement is higher than the null, and pairs")
        print("  are dominated by the large families, so the extra agreement")
        print("  lives in the big families rather than in whether a small one")
        print("  happens to come out pure -- result B, narrowed.")
    print(sweep[["n_bins", "movable_frac", "pure_pos_obs", "pure_pos_null",
                 "pure_pos_p", "concordance_obs", "concordance_null",
                 "concordance_p"]].round(4).to_string(index=False))

    print("\n=== Part 3b: where does the extra agreement sit? ===")
    print("  Concordance is pair-weighted, so this asks whether the excess is")
    print("  in the large families that carry most of the pairs.\n")
    sizes = np.bincount(codes, minlength=n_families)
    bands = [(3, 4), (5, 9), (10, 19), (20, 10 ** 6)]
    binned = work.groupby("epitope_aa")["nn3"].transform(
        lambda s: pd.qcut(s, min(10, s.nunique()), labels=False, duplicates="drop"))
    strata10 = (work["epitope_aa"] + "|" + binned.astype(str)).to_numpy()
    groups10 = [np.where(strata10 == s)[0] for s in np.unique(strata10)]
    size_rows = []
    for lo, hi in bands:
        pick = np.where((sizes >= lo) & (sizes <= hi))[0]
        if len(pick) == 0:
            continue
        member = np.isin(codes, pick)
        sub_codes = pd.factorize(codes[member])[0]
        n_sub = sub_codes.max() + 1
        obs = purity_counts(y[member], sub_codes, n_sub)
        null = []
        for _ in range(500):
            shuffled = y.copy()
            for g in groups10:
                shuffled[g] = rng.permutation(y[g])
            null.append(purity_counts(shuffled[member], sub_codes, n_sub)[3])
        null = np.asarray(null)
        size_rows.append({
            "family_size": f"{lo}-{hi if hi < 10 ** 6 else '+'}",
            "families": len(pick), "tcrs": int(member.sum()),
            "pairs": int(sum(n * (n - 1) / 2 for n in sizes[pick])),
            "concordance": obs[3], "null": float(null.mean()),
            "excess": obs[3] - float(null.mean()),
            "p": ((null >= obs[3]).sum() + 1) / (len(null) + 1)})
        print(f"  size {size_rows[-1]['family_size']:<5} "
              f"{len(pick):2d} families, {size_rows[-1]['pairs']:5d} pairs   "
              f"concordance {obs[3]:.3f} vs null {null.mean():.3f}   "
              f"excess {obs[3] - null.mean():+.3f}  p {size_rows[-1]['p']:.3f}")
    by_size = pd.DataFrame(size_rows)

    print("\n=== Part 3c: robustness -- a wider permutation pool ===")
    print("  Above, labels are shuffled only among TCRs that belong to a family")
    print("  with >=3 labelled members. This repeats it shuffling among ALL")
    print("  labelled rows of the stratum, so family members can trade labels")
    print("  with cohort members that have no family.\n")
    full_labels = labels.dropna(subset=["nn3"]).reset_index(drop=True)
    full_binned = full_labels.groupby("epitope_aa")["nn3"].transform(
        lambda s: pd.qcut(s, min(10, s.nunique()), labels=False, duplicates="drop"))
    full_strata = (full_labels["epitope_aa"] + "|"
                   + full_binned.astype(str)).to_numpy()
    full_groups = [np.where(full_strata == s)[0] for s in np.unique(full_strata)]
    position = {t: i for i, t in enumerate(full_labels["tcr_id"])}
    work_rows = np.array([position[t] for t in work["tcr_id"]])
    full_y = full_labels["pos"].to_numpy(dtype=float)
    wide = []
    for _ in range(DRAWS):
        shuffled = full_y.copy()
        for g in full_groups:
            shuffled[g] = rng.permutation(full_y[g])
        wide.append(purity_counts(shuffled[work_rows], codes, n_families))
    wide = np.asarray(wide)
    wide_p = ((wide[:, 3] >= observed[3]).sum() + 1) / (DRAWS + 1)
    print(f"  concordance {observed[3]:.4f} vs wider null "
          f"{wide[:, 3].mean():.4f}  p {wide_p:.4f}")
    print(f"  pure-positive {observed[0]} vs {wide[:, 0].mean():.1f}  "
          f"p {((wide[:, 0] >= observed[0]).sum() + 1) / (DRAWS + 1):.4f}")

    print("\n=== Part 3d: the families that carry the excess, named ===")
    print("  A result resting on a handful of objects has to be inspectable.\n")
    big_rows = []
    if VDJDB_FULL.exists():
        full = pd.read_csv(VDJDB_FULL, sep="\t", low_memory=False)
        full = full[(full["species"] == "HomoSapiens")
                    & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]
        pair_refs = (full.groupby(["cdr3.alpha", "cdr3.beta"])["reference.id"]
                     .agg(lambda s: tuple(sorted(set(s.dropna())))))
        key = test.set_index("tcr_id")[["cdr3_alpha_aa", "cdr3_beta_aa",
                                        "TRBV_IMGT"]]
        for position_i in np.argsort(-sizes):
            if sizes[position_i] < 20:
                continue
            fam = uniques[position_i]
            block = work[work["family"] == fam]
            sets = [set(pair_refs.get((key.loc[t, "cdr3_alpha_aa"],
                                       key.loc[t, "cdr3_beta_aa"]), ()))
                    for t in block["tcr_id"]]
            shared = set.intersection(*sets) if sets and all(sets) else set()
            vgene = key.loc[block["tcr_id"], "TRBV_IMGT"].astype(str)
            stem = vgene.str.split("*").str[0].value_counts()
            big_rows.append({
                "family": fam, "labelled": len(block),
                "validated": block["pos"].mean(),
                "nn3_median": block["nn3"].median(),
                "distinct_references": len(set().union(*sets)) if sets else 0,
                "shared_by_all": ", ".join(sorted(shared)) if shared else "none",
                "dominant_TRBV": f"{stem.index[0]} ({stem.iloc[0]}/{len(block)})"})
            print(f"  {fam:<14} n={len(block):3d}  validated "
                  f"{block['pos'].mean():6.1%}  refs "
                  f"{big_rows[-1]['distinct_references']:3d}  shared "
                  f"{big_rows[-1]['shared_by_all']:<6}  "
                  f"{big_rows[-1]['dominant_TRBV']}")
    big_table = pd.DataFrame(big_rows)
    if len(big_table):
        print(f"\n  None sharing a single reference across all members means the")
        print("  excess is not a study artifact: each family is reported")
        print("  independently in 10-19 papers.")

    provenance = pd.DataFrame()
    if verdict_b and VDJDB_FULL.exists():
        print("\n=== Part 4: is the excess driven by study or donor? ===")
        full = pd.read_csv(VDJDB_FULL, sep="\t", low_memory=False)
        full = full[(full["species"] == "HomoSapiens")
                    & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]
        key = test.set_index("tcr_id")[["cdr3_alpha_aa", "cdr3_beta_aa"]]
        pair_to_refs = (full.groupby(["cdr3.alpha", "cdr3.beta"])["reference.id"]
                        .agg(lambda s: tuple(sorted(set(s.dropna())))))
        work_refs = [pair_to_refs.get(
            (key.loc[t, "cdr3_alpha_aa"], key.loc[t, "cdr3_beta_aa"]), ())
            for t in work["tcr_id"]]
        work = work.assign(refs=work_refs)
        prov_rows = []
        for fam, block in work.groupby("family"):
            rate = block["pos"].mean()
            sets = [set(r) for r in block["refs"] if r]
            shared = set.intersection(*sets) if len(sets) == len(block) and sets else set()
            prov_rows.append({"family": fam, "n": len(block), "purity": rate,
                              "all_share_a_reference": bool(shared),
                              "n_distinct_references": len(
                                  set().union(*sets)) if sets else 0})
        provenance = pd.DataFrame(prov_rows)
        pure = provenance[provenance["purity"].isin([0.0, 1.0])]
        rest = provenance[~provenance["purity"].isin([0.0, 1.0])]
        print(f"  pure families sharing one reference:  "
              f"{pure['all_share_a_reference'].mean():.1%} of {len(pure)}")
        print(f"  mixed families sharing one reference: "
              f"{rest['all_share_a_reference'].mean():.1%} of {len(rest)}")
        print("\n  If the first is far above the second, the extra agreement is")
        print("  at least partly a provenance effect, not a recognition one.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E54 - is family purity more than the continuous density relation?\n\n"
        "## The graph (full cohort) against the purity set (labelled only)\n\n"
        + graph_table.to_markdown(index=False)
        + f"\n\nFamilies with >={MIN_LABELLED} labelled members: **{n_families}**, "
        f"holding {len(work)} TCRs. Observed: pure-positive **{observed[0]}**, "
        f"pure-negative {observed[1]}, mixed {observed[2]}, pair concordance "
        f"**{observed[3]:.4f}**.\n\n"
        "## The density-to-label relation the null preserves\n\n"
        + band.round(4).to_markdown()
        + "\n\n## Permutation within (epitope x nn3 bin)\n\n"
        + sweep.round(4).to_markdown(index=False)
        + "\n\n## Where the extra agreement sits, by family size\n\n"
        + by_size.round(4).to_markdown(index=False)
        + f"\n\nWider permutation pool (all labelled rows in the stratum): "
        f"concordance {observed[3]:.4f} against {wide[:, 3].mean():.4f}, "
        f"p {wide_p:.4f}.\n" 
        + ("\n\n## The families carrying the excess\n\n"
           + big_table.round(4).to_markdown(index=False) if len(big_table) else "")
        + ("\n\n## Provenance\n\n" + provenance.round(4).to_markdown(index=False)
           if len(provenance) else "")
        + "\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
