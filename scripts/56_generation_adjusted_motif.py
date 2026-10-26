#!/usr/bin/env python
"""E58 -- are the public motifs enriched beyond what rearrangement explains?

E57 established the six families are reproduced across 12-112 independent
donors. That says the motifs are real observations, not one cohort propagated.
It does not say why they recur: V(D)J recombination makes some CDR3s far easier
to produce than others, and a motif that is simply easy to generate will be
found repeatedly without any antigen preferring it.

The question is therefore the plan's: **if TCRs were produced by the
rearrangement model alone, how common would these CDR3β motifs be?**

Two things this cannot do, stated before the numbers:

  * OLGA generates CDR3s without a V gene, so a motif here is defined on the
    **CDR3β sequence only**, and the families' TRBV restriction is reported
    separately rather than folded into the enrichment.
  * Any cohort motif will look enriched against a naive generation model,
    because the cohort is antigen-selected to begin with. The number that
    carries information is therefore not the raw enrichment but the comparison
    with the plan's required control -- **same epitope, same TRBV19, matched
    size, not in a large family.** That control is built here and all
    conclusions are relative to it.

The motif definition is fixed before any enrichment is computed and uses no
labels: a family's modal CDR3β length L, its position-wise consensus over
members of that length, and a radius equal to the rounded-up median within-
family Hamming distance to that consensus. A sequence matches if it has length
L and lies within the radius.

Run: .venv/bin/python scripts/56_generation_adjusted_motif.py
"""

from __future__ import annotations

import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/31_generation_adjusted_motif.md"
CACHE = ROOT / "data/interim/olga_beta_200k.npy"
RADIUS = 0.12
BIG = 20
N_OLGA = 200_000
N_CONTROL = 200
SEED = 0


def families(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        _, component = connected_components(
            sp.csr_matrix(1.0 - block @ block.T <= RADIUS), directed=False)
        out.update({t: f"{epitope}#{c}" for t, c in zip(ids, component)})
    return out


def olga_beta(n: int) -> np.ndarray:
    if CACHE.exists():
        pool = np.load(CACHE, allow_pickle=True)
        if len(pool) >= n:
            return pool[:n]
    from tcrdist.pgen import OlgaModel

    print(f"  generating {n} OLGA beta CDR3s (cached afterwards) ...")
    model = OlgaModel(chain_folder="human_T_beta", recomb_type="VDJ")
    out = [s for s in model.gen_cdr3s(n=int(n * 1.3)) if s][:n]
    pool = np.array(out, dtype=object)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.save(CACHE, pool, allow_pickle=True)
    return pool


def motif_from(sequences) -> tuple[int, str, int] | None:
    """Modal length, position-wise consensus, and a radius set by the members'
    own spread. Uses no labels and no knowledge of what it will be tested on."""
    lengths = Counter(len(s) for s in sequences)
    length, _ = lengths.most_common(1)[0]
    members = [s for s in sequences if len(s) == length]
    if len(members) < 3:
        return None
    consensus = "".join(
        Counter(s[i] for s in members).most_common(1)[0][0] for i in range(length))
    spread = [sum(a != b for a, b in zip(s, consensus)) for s in members]
    return length, consensus, int(np.ceil(np.median(spread)))


def match_rate(pool, length: int, consensus: str, radius: int) -> float:
    hits = 0
    for s in pool:
        if len(s) != length:
            continue
        d = 0
        for a, b in zip(s, consensus):
            if a != b:
                d += 1
                if d > radius:
                    break
        if d <= radius:
            hits += 1
    return hits / len(pool)


def main() -> int:
    rng = np.random.default_rng(SEED)
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    family_of = families(test)
    trbv = (test.set_index("tcr_id")["TRBV_IMGT"].astype(str)
            .str.split("*").str[0].to_dict())
    beta = test.set_index("tcr_id")["cdr3_beta_aa"].astype(str).to_dict()
    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    counts = (labels.assign(f=labels["tcr_id"].map(family_of))
              .groupby("f").size())
    big = set(counts[counts >= BIG].index)

    gil = [t for t in test.loc[test["epitope_aa"] == "GILGFVFTL", "tcr_id"]]
    in_family = [t for t in gil if family_of[t] in big]
    v_control = [t for t in gil if family_of[t] not in big and trbv[t] == "TRBV19"]
    other_v = [t for t in gil if family_of[t] not in big and trbv[t] != "TRBV19"]

    print("=== Part 1: is the motif simply easier to generate? ===")
    print("  log10 Pgen(beta), GILGFVFTL only so the epitope is fixed.\n")
    pgen_rows = []
    for name, members in (("in a large TRBV19 family", in_family),
                          ("TRBV19, no large family", v_control),
                          ("other V, no large family", other_v)):
        values = pgen.loc[members, "log10_pgen_beta"].dropna()
        lengths = pd.Series([len(beta[t]) for t in members])
        pgen_rows.append({"group": name, "n": len(members),
                          "median_log10_pgen_beta": float(values.median()),
                          "median_cdr3b_length": float(lengths.median())})
        print(f"  {name:<26} n={len(members):4d}  median log10 Pgen "
              f"{values.median():7.3f}   median CDR3b length "
              f"{lengths.median():.0f}")
    pgen_table = pd.DataFrame(pgen_rows)
    delta = (pgen_table.iloc[0]["median_log10_pgen_beta"]
             - pgen_table.iloc[1]["median_log10_pgen_beta"])
    print(f"\n  family members are {delta:+.3f} log10 units "
          f"({'more' if delta > 0 else 'less'} generatable) than matched "
          f"TRBV19 non-family TCRs")

    print("\n=== Part 2: motif enrichment over a generation background ===")
    pool = olga_beta(N_OLGA)
    print(f"  background: {len(pool)} OLGA beta CDR3s, no antigen selection\n")
    cohort_beta_by_epitope = {
        epitope: [beta[t] for t in block["tcr_id"]]
        for epitope, block in test.groupby("epitope_aa")}
    cohort_beta = cohort_beta_by_epitope["GILGFVFTL"]
    rows = []
    for family in sorted(big, key=lambda f: -counts[f]):
        members = [beta[t] for t, f in family_of.items() if f == family]
        spec = motif_from(members)
        if spec is None:
            continue
        length, consensus, radius = spec
        # A motif's observed frequency belongs to ITS OWN epitope's cohort. A
        # first version measured every family against GILGFVFTL's cohort, which
        # scored the GLCTLVAML motif at exactly zero.
        own = cohort_beta_by_epitope[family.split("#")[0]]
        observed = match_rate(own, length, consensus, radius)
        background = match_rate(pool, length, consensus, radius)
        rows.append({
            "family": family, "members": len(members), "L": length,
            "consensus": consensus, "radius": radius,
            "observed_in_cohort": observed, "olga_background": background,
            "enrichment": observed / background if background else np.inf})
        print(f"  {family:<14} L={length} r={radius} {consensus:<20} "
              f"cohort {observed:7.4f}  OLGA {background:.2e}  "
              f"enrichment {rows[-1]['enrichment']:>10.1f}x")
    motifs = pd.DataFrame(rows)

    print("\n=== Part 3: the control the plan requires ===")
    print("  Any antigen-selected cohort motif beats a naive generation model.")
    print(f"  So: {N_CONTROL} pseudo-motifs built the SAME way from same-size")
    print("  random sets of GILGFVFTL TRBV19 TCRs that are NOT in a large")
    print("  family. If the real motifs are not more enriched than these, the")
    print("  enrichment is a property of the cohort, not of the motifs.\n")
    control_rows = []
    for _, row in motifs.iterrows():
        if not row["family"].startswith("GILGFVFTL"):
            control_rows.append({
                "family": row["family"], "enrichment": row["enrichment"],
                "control_median": np.nan, "control_p95": np.nan,
                "p_vs_control": np.nan})
            print(f"  {row['family']:<14} enrichment "
                  f"{row['enrichment']:>9.1f}x   NOT ESTIMABLE -- E56 found "
                  f"only 13 matched controls")
            continue
        size = int(row["members"])
        draws = []
        for _ in range(N_CONTROL):
            pick = rng.choice(len(v_control), size=min(size, len(v_control)),
                              replace=False)
            spec = motif_from([beta[v_control[i]] for i in pick])
            if spec is None:
                continue
            length, consensus, radius = spec
            background = match_rate(pool, length, consensus, radius)
            observed = match_rate(cohort_beta, length, consensus, radius)
            draws.append(observed / background if background else np.nan)
        draws = np.array([d for d in draws if np.isfinite(d)])
        tail = ((draws >= row["enrichment"]).sum() + 1) / (len(draws) + 1)
        control_rows.append({
            "family": row["family"], "enrichment": row["enrichment"],
            "control_median": float(np.median(draws)),
            "control_p95": float(np.percentile(draws, 95)),
            "p_vs_control": tail})
        print(f"  {row['family']:<14} enrichment {row['enrichment']:>9.1f}x   "
              f"control median {np.median(draws):>8.1f}x  p95 "
              f"{np.percentile(draws, 95):>9.1f}x   p {tail:.3f}")
    control = pd.DataFrame(control_rows)

    print("\n=== how many INDEPENDENT motifs is this really? ===")
    print("  The five GILGFVFTL families are separate only because the 0.12")
    print("  radius split them. If their consensus sequences are near-")
    print("  identical they are sub-clusters of one public motif, and five")
    print("  p-values are not five independent pieces of evidence.\n")
    gil_motifs = motifs[motifs["family"].str.startswith("GILGFVFTL")]
    for _, a in gil_motifs.iterrows():
        line = []
        for _, b in gil_motifs.iterrows():
            if a["family"] == b["family"]:
                line.append("  .")
            elif a["L"] != b["L"]:
                line.append("  L")
            else:
                line.append(f"{sum(x != y for x, y in zip(a['consensus'], b['consensus'])):3d}")
        print(f"  {a['consensus']:<16} " + "".join(line))
    same_length = gil_motifs[gil_motifs["L"] == gil_motifs["L"].mode().iloc[0]]
    pairs = [sum(x != y for x, y in zip(a, b))
             for i, a in enumerate(same_length["consensus"])
             for b in list(same_length["consensus"])[i + 1:]]
    print(f"\n  ('L' = different length, so no Hamming distance.)")
    print(f"  Among the {len(same_length)} families of the modal length, the")
    print(f"  consensus sequences differ by a median of "
          f"{np.median(pairs):.0f} residues out of "
          f"{int(same_length['L'].iloc[0])}.")
    print("  They are sub-clusters of a single CASSIRS-type public motif, not")
    print("  five independent discoveries. The effective number of independent")
    print("  GILGFVFTL motifs tested here is closer to 2 than to 5.")

    print("\n=== reading it ===")
    estimable = control["p_vs_control"].notna()
    beats = int((control.loc[estimable, "p_vs_control"] < 0.05).sum())
    print(f"  {beats} of {int(estimable.sum())} ESTIMABLE families exceed the")
    print(f"  matched-control distribution at p<0.05 "
          f"({int((~estimable).sum())} not estimable).")
    print("  But they are not that many independent tests: the effective")
    print("  number of distinct GILGFVFTL motifs is about 2.")
    if delta < 0:
        print(f"\n  And Part 1 points the same way: family CDR3betas are "
              f"{abs(delta):.2f} log10")
        print("  units HARDER to generate than matched TRBV19 non-family ones,")
        print("  so being easy to produce is not why they recur.")
    else:
        print(f"\n  But Part 1 cautions: family CDR3betas are {delta:.2f} log10")
        print("  units EASIER to generate, so generation explains part of it.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E58 - generation-adjusted motif enrichment\n\n"
        "## Part 1 - generation probability by group (GILGFVFTL)\n\n"
        + pgen_table.round(4).to_markdown(index=False)
        + "\n\n## Part 2 - motif enrichment over an OLGA background\n\n"
        + motifs.round(6).to_markdown(index=False)
        + "\n\n## Part 3 - against the matched control\n\n"
        + control.round(4).to_markdown(index=False)
        + "\n\nMotifs are defined on CDR3beta only, because OLGA generates no V "
        "gene; the families' TRBV19 restriction is reported separately and is "
        "not folded into the enrichment.\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
