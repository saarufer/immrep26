#!/usr/bin/env python
"""E129 -- Is alpha-beta pairing compatibility anchor-specific or generic?

The collaborator's E129, on Minervina's exo-TCR chain-replacement panels
(Zenodo 10.5281/zenodo.14885005). Their step 4 asks for three baselines --
single-chain distance, the sum of independent alphabeta distances, and a model
carrying a genuine alphabeta interaction term -- under anchor holdout.

**The identifiability problem that reshapes the experiment.** The anchor chain
is not a column in these files; it is implicit, one fixed exo-TCR per panel.
Within a panel the anchor is therefore *constant*, so an interaction term
between anchor and partner is unidentifiable there: nothing varies to interact
with. An interaction model needs several anchors, and it needs the anchor
*sequences*, which this release does not carry.

So the question is asked in the form the data can actually answer, which turns
out to be the paper's own title claim. If pairing compatibility is a real
interaction, then two different anchors should admit *different* partner sets
beyond what partner-chain genericness explains. That is testable from the
labels alone:

  lift(X -> Y) = P(partner is permissive for Y | permissive for X)
                 / P(partner is permissive for Y)

  lift ~ 1   pairing is anchor-specific; knowing X tells you nothing about Y,
             and the competition's natural alphabeta receptors need a genuine
             interaction model.
  lift >> 1  partner permissiveness is largely generic, and E126's
             single-chain margin already reads it -- there is no further
             channel here.

The decisive comparison is **GILb vs GILb_J12**: two different fixed-beta
anchors, the *same* epitope (GILGFVFTL), the *same* partner chain type (alpha),
the same donor pool. Any lift between those two cannot be epitope specificity.
The ten alpha-partner panels give the full matrix, splitting same-epitope pairs
from cross-epitope ones.

**Controls.** (1) Donor. A partner shared between two panels may be shared
because it is a public clonotype in one donor's repertoire, so every lift is
recomputed within donor. (2) The negative sets are the donors' repertoires
minus the enriched clones, 33-191 MB each, so the base rate P(permissive) is
small and a raw co-occurrence count would be meaningless -- hence lift rather
than concordance. (3) Expression: these partners are all sequenced from cells,
so "not expressed" and "expressed but not binding" are not confounded the way
they would be in a transfection panel -- the negatives are observed, expressed
chains that failed to enrich. Read depth is reported as the QC.

Run: .venv/bin/python scripts/112_chain_pairing.py
"""

from __future__ import annotations

import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
PANEL = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
             "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/pn/pos_neg")
REPORT = ROOT / "reports/90_chain_pairing.md"
UC = ["CDR3.amino.acid.sequence.x", "bestVGene.x", "donor.x",
      "Read.count.x", "FC"]
# the ten fixed-beta anchors (partner = alpha), with their epitope
ALPHA_PARTNER = {"GILb": "GILGFVFTL", "GILb_J12": "GILGFVFTL",
                 "LTDb": "LTDEMIAQY", "LLWb": "LLWNGPMAV",
                 "TTDb": "TTDPSFLGRY", "YLQb": "YLQPRTFLL",
                 "EIYb": "EIYKRWII", "FLCb": "FLCMKALL",
                 "GTSb": "GTSGSPIIDKK", "SFHb": "SFHSLHLLF"}


def load_panel(tag: str) -> pd.DataFrame:
    p = pd.read_csv(PANEL / f"{tag}_true_pos.tsv", sep="\t",
                    usecols=UC, low_memory=False)
    n = pd.read_csv(PANEL / f"{tag}_true_neg.tsv", sep="\t",
                    usecols=UC, low_memory=False)
    p["y"], n["y"] = 1, 0
    d = pd.concat([p, n], ignore_index=True)
    d = d.rename(columns={"CDR3.amino.acid.sequence.x": "cdr3",
                          "donor.x": "donor", "Read.count.x": "reads"})
    d = d.dropna(subset=["cdr3"])
    # a clone seen as both: trust the enriched call
    return (d.groupby(["cdr3", "donor"], as_index=False)
             .agg(y=("y", "max"), reads=("reads", "max")))


def lift(a: pd.DataFrame, b: pd.DataFrame) -> tuple:
    """lift of permissiveness in b given permissiveness in a, shared partners."""
    m = a.merge(b, on="cdr3", suffixes=("_a", "_b"))
    if len(m) < 50:
        return (np.nan,) * 5
    base = m["y_b"].mean()
    given = m.loc[m["y_a"] == 1, "y_b"]
    if len(given) < 10 or base == 0:
        return (np.nan,) * 5
    return (len(m), int(m["y_a"].sum()), base, given.mean(),
            given.mean() / base)


def main() -> int:
    print("=== loading the ten alpha-partner panels ===")
    P: dict[str, pd.DataFrame] = {}
    for tag in ALPHA_PARTNER:
        P[tag] = load_panel(tag)
        d = P[tag]
        print(f"  {tag:<10}{ALPHA_PARTNER[tag]:<13}{len(d):8d} clone-donor rows"
              f"  {int(d.y.sum()):6d} permissive ({d.y.mean()*100:.2f}%)"
              f"  {d.donor.nunique()} donors")

    print("\n=== 1. the decisive pair: two anchors, one epitope ===")
    n, k, base, cond, lf = lift(P["GILb"], P["GILb_J12"])
    print(f"  GILb -> GILb_J12   shared partners {n}  permissive in GILb {k}")
    print(f"    P(permissive in J12)            = {base*100:6.3f}%")
    print(f"    P(permissive in J12 | in GILb)  = {cond*100:6.3f}%")
    print(f"    LIFT = {lf:.2f}x")
    n2, k2, b2, c2, l2 = lift(P["GILb_J12"], P["GILb"])
    print(f"  GILb_J12 -> GILb   LIFT = {l2:.2f}x  "
          f"({b2*100:.3f}% -> {c2*100:.3f}%)")

    print("\n=== 2. full matrix of lifts (row -> column) ===")
    tags = list(ALPHA_PARTNER)
    M = pd.DataFrame(np.nan, index=tags, columns=tags)
    for x in tags:
        for y in tags:
            if x == y:
                continue
            M.loc[x, y] = lift(P[x], P[y])[4]
    print(M.round(2).to_string())

    print("\n=== 3. same epitope vs different epitope ===")
    same, diff = [], []
    for x, y in combinations(tags, 2):
        for a, b in ((x, y), (y, x)):
            v = M.loc[a, b]
            if np.isnan(v):
                continue
            (same if ALPHA_PARTNER[a] == ALPHA_PARTNER[b] else diff).append(v)
    print(f"  same epitope, different anchor : n={len(same):3d}  "
          f"median lift {np.median(same):.2f}  range "
          f"{min(same):.2f}-{max(same):.2f}")
    print(f"  different epitope              : n={len(diff):3d}  "
          f"median lift {np.median(diff):.2f}  range "
          f"{min(diff):.2f}-{max(diff):.2f}")

    print("\n=== 4. control: within donor ===")
    rows = []
    for a, b in [("GILb", "GILb_J12"), ("GILb_J12", "GILb")]:
        for dn in sorted(set(P[a].donor) & set(P[b].donor)):
            va = P[a][P[a].donor == dn]
            vb = P[b][P[b].donor == dn]
            n, k, bs, cd, lv = lift(va, vb)
            if not np.isnan(lv):
                print(f"  {a} -> {b}  donor {dn}: shared {n:6d}  "
                      f"base {bs*100:6.3f}%  cond {cd*100:6.3f}%  lift {lv:.2f}x")
                rows.append([a, b, dn, n, bs, cd, lv])

    print("\n=== 5. QC: read depth of permissive vs not ===")
    for tag in ("GILb", "GILb_J12", "LTDb"):
        d = P[tag]
        print(f"  {tag:<10}permissive median reads "
              f"{d.loc[d.y==1,'reads'].median():8.0f}   "
              f"non-permissive {d.loc[d.y==0,'reads'].median():8.0f}")

    lines = ["# E129 — Is αβ pairing compatibility anchor-specific or generic?",
             "", "Generated by `scripts/112_chain_pairing.py`.", "",
             "The anchor chain is implicit (one exo-TCR per panel) and its",
             "sequence is not in this release, so an anchor×partner interaction",
             "term is unidentifiable. The answerable form is the paper's own",
             "title claim, via lift(X→Y) = P(permissive for Y | permissive for X)",
             "/ P(permissive for Y).", "",
             "## Panels", "",
             "| anchor | epitope | clone-donor rows | permissive | % | donors |",
             "|---|---|---:|---:|---:|---:|"]
    for tag in ALPHA_PARTNER:
        d = P[tag]
        lines.append(f"| {tag} | {ALPHA_PARTNER[tag]} | {len(d)} "
                     f"| {int(d.y.sum())} | {d.y.mean()*100:.2f}% "
                     f"| {d.donor.nunique()} |")
    lines += ["", "## The decisive pair (same epitope, two anchors)", "",
              f"- `GILb → GILb_J12`: {n} shared partners, "
              f"P(permissive) {base*100:.3f}% → {cond*100:.3f}% given, "
              f"**lift {lf:.2f}×**",
              f"- `GILb_J12 → GILb`: **lift {l2:.2f}×** "
              f"({b2*100:.3f}% → {c2*100:.3f}%)", "",
              "## Full lift matrix", "", "```",
              M.round(2).to_string(), "```", "",
              "## Same epitope vs different epitope", "",
              f"- same epitope, different anchor: n={len(same)}, "
              f"median **{np.median(same):.2f}×**",
              f"- different epitope: n={len(diff)}, "
              f"median **{np.median(diff):.2f}×**", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
