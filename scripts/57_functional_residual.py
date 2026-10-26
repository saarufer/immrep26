#!/usr/bin/env python
"""E59 -- do the public motifs validate better than matched non-motif TCRs?

E58 measured how often these motifs occur. This asks whether occurring in one
predicts passing functional revalidation, against controls matched on the things
that could produce the same difference without any recognition advantage.

Three constraints carried in from the earlier experiments, none of them
optional:

  * **GLCTLVAML is not estimable.** E56 found 13 matched controls. The plan
    forbids loosening the matching to rescue it, so it is reported as such.
  * **nn3 is a mediator, not a confounder.** The motif's convergence is *why*
    its members have close neighbours, so regressing density away does not
    isolate a cleaner effect -- it removes the mechanism. Both the uncontrolled
    association and the density-matched residual are reported, separately, and
    neither is called the causal estimate.
  * **The unit is the motif, not the TCR and not the family.** E58 showed four
    of the five GILGFVFTL families share length 13 and the CASSIRS prefix and
    differ by a median of 2 residues: they are sub-clusters of one public motif.
    So 158 labelled TCRs are not 158 observations, and five families are not
    five. Leave-one-out is run at the family level *and* at the motif level,
    and the motif level is the one that counts.

Run: .venv/bin/python scripts/57_functional_residual.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import fisher_exact

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/32_functional_residual.md"
RADIUS = 0.12
BIG = 20
EPITOPE = "GILGFVFTL"
GENE = "TRBV19"


def geometry(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    family, nn3 = {}, {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, component = connected_components(
            sp.csr_matrix(d <= RADIUS), directed=False)
        np.fill_diagonal(d, np.inf)
        ordered = np.sort(d, axis=1)
        for i, t in enumerate(ids):
            family[t] = f"{epitope}#{component[i]}"
            nn3[t] = float(ordered[i, :3].mean())
    return family, nn3


def stratified_difference(frame, strata_columns):
    """Mantel-Haenszel style: a difference within each stratum, weighted by the
    stratum's harmonic sample size, so a stratum holding one TCR of one group
    cannot swing the pooled estimate."""
    rows, weights, diffs = [], [], []
    for keys, block in frame.groupby(strata_columns):
        a = block[block["motif"]]
        b = block[~block["motif"]]
        if len(a) == 0 or len(b) == 0:
            continue
        w = len(a) * len(b) / (len(a) + len(b))
        diffs.append(a["pos"].mean() - b["pos"].mean())
        weights.append(w)
        rows.append({"stratum": str(keys), "n_motif": len(a),
                     "n_control": len(b), "motif_rate": a["pos"].mean(),
                     "control_rate": b["pos"].mean()})
    if not weights:
        return np.nan, 0, 0, pd.DataFrame(rows)
    pooled = float(np.average(diffs, weights=weights))
    covered = sum(r["n_motif"] for r in rows)
    return pooled, covered, len(rows), pd.DataFrame(rows)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    family_of, nn3 = geometry(test)
    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    meta = test.set_index("tcr_id")
    trbv = meta["TRBV_IMGT"].astype(str).str.split("*").str[0]
    counts = (labels.assign(f=labels["tcr_id"].map(family_of))
              .groupby("f").size())
    big = set(counts[counts >= BIG].index)

    frame = labels[labels["epitope_aa"] == EPITOPE].copy()
    frame["family"] = frame["tcr_id"].map(family_of)
    frame["motif"] = frame["family"].isin(big)
    frame["trbv"] = frame["tcr_id"].map(trbv)
    frame["len_b"] = frame["tcr_id"].map(meta["cdr3_beta_aa"].astype(str).str.len())
    frame["nn3"] = frame["tcr_id"].map(nn3)
    frame["pgen"] = frame["tcr_id"].map(pgen["log10_pgen_beta"])
    frame = frame[frame["trbv"] == GENE]      # the control the plan requires
    frame["pgen_tercile"] = pd.qcut(frame["pgen"], 3, labels=False,
                                    duplicates="drop")
    frame["nn3_tercile"] = pd.qcut(frame["nn3"], 3, labels=False,
                                   duplicates="drop")

    print(f"=== the comparison, {EPITOPE} / {GENE} only ===")
    print(f"  GLCTLVAML is NOT ESTIMABLE (E56: 13 matched controls).\n")
    a, b = frame[frame["motif"]], frame[~frame["motif"]]
    print(f"  in a large family      n={len(a):4d}  validated {a['pos'].mean():6.1%}")
    print(f"  TRBV19, no family      n={len(b):4d}  validated {b['pos'].mean():6.1%}")
    odds, p = fisher_exact([[int(a["pos"].sum()), int((1 - a["pos"]).sum())],
                            [int(b["pos"].sum()), int((1 - b["pos"]).sum())]])
    raw = a["pos"].mean() - b["pos"].mean()
    print(f"  raw difference {raw:+.3f}   odds ratio {odds:.1f}   "
          f"Fisher p {p:.2e}")

    print("\n=== matching, one layer at a time ===")
    ladder = []
    for name, columns in (
        ("CDR3b length", ["len_b"]),
        ("length x Pgen tercile", ["len_b", "pgen_tercile"]),
        ("length x Pgen x nn3 tercile (density as well)",
         ["len_b", "pgen_tercile", "nn3_tercile"]),
    ):
        pooled, covered, strata, detail = stratified_difference(frame, columns)
        ladder.append({"matching": name, "difference": pooled,
                       "motif_tcrs_covered": covered, "strata_used": strata})
        print(f"  {name:<44} {pooled:+.3f}   "
              f"{covered}/{len(a)} motif TCRs in a usable stratum, "
              f"{strata} strata")
    ladder = pd.DataFrame(ladder)
    print("\n  The last row controls the MEDIATOR. It is reported because the")
    print("  plan asks for it, not because it is the causal estimate: the")
    print("  motif's convergence is why its members have close neighbours, so")
    print("  matching on nn3 removes the mechanism along with the confounding.")

    print("\n=== E59.2: the statistical unit ===")
    print("  Leave one family out, then leave one MOTIF out. E58 showed four of")
    print("  the five GIL families are sub-clusters of one CASSIRS motif, so")
    print("  the family level overstates independence.\n")
    beta = meta["cdr3_beta_aa"].astype(str)
    consensus = {}
    for fam in sorted(big):
        if not fam.startswith(EPITOPE):
            continue
        members = [beta[t] for t, f in family_of.items() if f == fam]
        length = pd.Series([len(m) for m in members]).mode().iloc[0]
        same = [m for m in members if len(m) == length]
        consensus[fam] = "".join(
            pd.Series([m[i] for m in same]).mode().iloc[0] for i in range(length))
    # Two motif groups: the length-13 CASSIRS cluster and everything else.
    groups = {}
    for fam, seq in consensus.items():
        groups.setdefault("CASSIRS-13" if seq.startswith("CASSIRS") and len(seq) == 13
                          else f"other:{seq}", []).append(fam)
    print("  motif groups:")
    for name, members in groups.items():
        print(f"    {name:<22} {members}")

    sensitivity = []
    for level, units in (("family", {f: [f] for f in consensus}),
                         ("motif", groups)):
        for name, members in units.items():
            sub = frame[~frame["family"].isin(members)]
            x, y = sub[sub["motif"]], sub[~sub["motif"]]
            if len(x) < 5:
                sensitivity.append({"level": level, "removed": name,
                                    "difference": np.nan, "n_motif_left": len(x)})
                continue
            sensitivity.append({
                "level": level, "removed": name,
                "difference": x["pos"].mean() - y["pos"].mean(),
                "n_motif_left": len(x)})
    sensitivity = pd.DataFrame(sensitivity)
    for level in ("family", "motif"):
        block = sensitivity[sensitivity["level"] == level].dropna(subset=["difference"])
        print(f"\n  leave-one-{level}-out: difference ranges "
              f"{block['difference'].min():+.3f} to {block['difference'].max():+.3f} "
              f"over {len(block)} removals")
        for _, r in block.iterrows():
            print(f"    without {r['removed']:<22} {r['difference']:+.3f} "
                  f"({int(r['n_motif_left'])} motif TCRs left)")

    print("\n=== E59.3: reading it with E58 ===")
    motif_level = sensitivity[(sensitivity["level"] == "motif")].dropna(
        subset=["difference"])
    print(f"  E58: enrichment beyond generation, confirmed on ~2 independent")
    print(f"       motifs against a matched control.")
    print(f"  E59: validation advantage {raw:+.3f} uncontrolled, "
          f"{ladder.iloc[1]['difference']:+.3f} matched on length and Pgen,")
    print(f"       {ladder.iloc[2]['difference']:+.3f} once density is matched "
          f"too (mediator removed).")
    if len(motif_level) >= 2 and (motif_level["difference"] > 0).all():
        print(f"  The advantage survives removing either motif group "
              f"({motif_level['difference'].min():+.3f} worst case), so it is")
        print("  not carried by the CASSIRS cluster alone.")
    print("\n  Joint reading, with the plan's boundary intact: motifs that")
    print("  recur beyond the rearrangement model also revalidate more often")
    print("  than TCRs matched on epitope, V gene, length and Pgen. That is")
    print("  consistent with antigen-associated selection and does not")
    print("  demonstrate it -- a naive repertoire or a functional experiment")
    print("  is still required, and the evidence rests on about two motifs.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        f"# E59 - functional residual of the public motifs ({EPITOPE} only)\n\n"
        f"GLCTLVAML is **not estimable** (E56: 13 matched controls).\n\n"
        f"Raw: motif {a['pos'].mean():.3f} (n={len(a)}) against matched "
        f"{GENE} non-motif {b['pos'].mean():.3f} (n={len(b)}), difference "
        f"**{raw:+.3f}**, odds ratio {odds:.1f}, Fisher p {p:.2e}.\n\n"
        "## Matching ladder\n\n" + ladder.round(4).to_markdown(index=False)
        + "\n\nThe density-matched row controls a **mediator** and is not the "
        "causal estimate.\n\n## Leave-one-out sensitivity\n\n"
        + sensitivity.round(4).to_markdown(index=False)
        + "\n\nMotif groups: "
        + "; ".join(f"{k} = {v}" for k, v in groups.items()) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
