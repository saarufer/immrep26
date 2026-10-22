#!/usr/bin/env python
"""Fold each TCR against anchor-matched decoy peptides (plan cards D08, C07).

Every structural read tried so far asks "how confident is this complex", and E14
showed that mostly reports how canonical the TCR is -- the MHC-peptide interface,
which contains no TCR at all, scored 0.647. A confidence number cannot
distinguish "this TCR recognises *this* peptide" from "this TCR folds into a
generic docking pose".

The relative version can. Folding the same TCR against its annotated peptide and
against decoys, then taking target minus decoy, cancels whatever is a property
of the TCR alone. If validated TCRs show a larger target-versus-decoy margin
than non-validating ones, that is peptide-specific structural evidence; if not,
the structural signal is generic and the whole line is closed.

**Decoys are anchor-matched on purpose.** A decoy that cannot occupy the groove
would make the margin measure peptide-MHC compatibility rather than TCR
specificity. Matching is on anchor *class* rather than exact residue, which is
what the allele actually constrains: HLA-A*02:01 takes an aliphatic residue at
P2 and at the C-terminus, while HLA-A*01:01 is anchored at P3 by an acidic
residue and at the C-terminus by an aromatic one. Requiring identical residues
would be stricter than biology and, as a first draft of this script discovered,
unsatisfiable -- no reported A*02:01 nonamer shares GILGFVFTL's P2 isoleucine
*and* its C-terminal leucine. Each decoy is a real, independently reported
epitope of the same allele and length.

Targets are not refolded -- the main 2,636-complex run already holds them -- so
this costs three folds per sampled TCR. Staged public-first: 60 TCRs per public
epitope, and the private extension only runs if the public result justifies it.

Run: python scripts/30_control_peptide_inputs.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.hla import load_hla_platforms  # noqa: E402

OUTDIR = ROOT / "data/interim/decoy/inputs"
MANIFEST = ROOT / "data/interim/decoy/manifest.csv"
PER_EPITOPE = 60
SEED = 0

# Real epitopes of the same HLA, length and anchor residues as each target.
DECOYS = {
    "GILGFVFTL": ["LLLDRLNQL", "FLYALALLL", "YLLEMLWRL"],
    "YLQPRTFLL": ["LLLDRLNQL", "LLYDANYFL", "YLLEMLWRL"],
    "GLCTLVAML": ["LLLDRLNQL", "FLYALALLL", "LLYDANYFL"],
    "NLVPMVATV": ["SLLMWITQV", "CLGGLLTMV", "LLFGYPVYV"],
    # A*01:01 decoys carry an acidic P3 and an aromatic C-terminus, like the
    # target; a first draft used FTSDYYQLY and LTSLTILQLY, both of which have a
    # serine at P3 and so fail the real anchor requirement.
    "TTDPSFLGRY": ["CTDDNALAYY", "ILDTAGKEEY", "QSDNGLDSDY"],
    # Private-epitope decoys, staged for phase 2 only.
    "CINGVCWTV": ["SLLMWITQV", "CLGGLLTMV", "LLFGYPVYV"],
    "LLWNGPMAV": ["SLLMWITQV", "CLGGLLTMV", "FLYNLLTRV"],
    "LTDEMIAQY": ["CTELKLSDY", "VTEHDTLLY", "PTDNYITTY"],
}

# What each allele actually constrains, used to verify the table above.
ANCHORS = {
    "A*02:01": {"p2": set("LMIVQAT"), "last": set("VLIAM")},
    "A*01:01": {"p3": set("DE"), "last": set("YF")},
}
HLA_OF = {
    "GILGFVFTL": "A*02:01", "YLQPRTFLL": "A*02:01", "GLCTLVAML": "A*02:01",
    "NLVPMVATV": "A*02:01", "CINGVCWTV": "A*02:01", "LLWNGPMAV": "A*02:01",
    "TTDPSFLGRY": "A*01:01", "LTDEMIAQY": "A*01:01",
}


def main() -> int:
    test = load_test()
    labels = build_training_frame()
    platforms = load_hla_platforms()

    # Verify the table against what the allele constrains, rather than trusting it.
    print("=== decoy anchor check (class, not identity) ===")
    for target, decoys in DECOYS.items():
        rule = ANCHORS[HLA_OF[target]]
        for decoy in decoys:
            checks = {"length": len(decoy) == len(target)}
            for site, allowed in rule.items():
                position = {"p2": 1, "p3": 2, "last": -1}[site]
                checks[f"{site}:{target[position]}/{decoy[position]}"] = (
                    target[position] in allowed and decoy[position] in allowed
                )
            ok = all(checks.values())
            detail = "  ".join(k for k, v in checks.items() if not v)
            print(f"  {target:<11} vs {decoy:<11} "
                  f"{'ok' if ok else 'FAIL: ' + detail}")
            if not ok:
                raise SystemExit(f"{decoy} is not a valid anchor-matched decoy for {target}")

    rng = np.random.default_rng(SEED)
    rows = []
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Phase 1: labelled epitopes only, balanced across the two classes.
    for epitope, cohort in labels.groupby("epitope_aa"):
        picks = []
        for label in (0, 1):
            pool = cohort[cohort["pos"] == label]
            take = min(PER_EPITOPE // 2, len(pool))
            picks.append(pool.sample(take, random_state=SEED))
        sample = pd.concat(picks)
        for row in sample.itertuples(index=False):
            for decoy in DECOYS[epitope]:
                name = f"{row.tcr_id}__{decoy}"
                target = OUTDIR / f"{name}.yaml"
                if not target.exists():
                    target.write_text(
                        "version: 1\nsequences:\n"
                        f"  - protein:\n      id: A\n      sequence: {row.full_seq_alpha_aa}\n      msa: empty\n"
                        f"  - protein:\n      id: B\n      sequence: {row.full_seq_beta_aa}\n      msa: empty\n"
                        f"  - protein:\n      id: M\n      sequence: {platforms[row.hla_long]}\n      msa: empty\n"
                        f"  - protein:\n      id: P\n      sequence: {decoy}\n      msa: empty\n"
                    )
                rows.append({"name": name, "tcr_id": row.tcr_id, "epitope": epitope,
                             "decoy": decoy, "pos": row.pos, "phase": 1})

    manifest = pd.DataFrame(rows)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, index=False)
    print(f"\nphase 1: {manifest.tcr_id.nunique()} TCRs x {len(DECOYS['GILGFVFTL'])} decoys "
          f"= {len(manifest)} folds")
    print(manifest.groupby(["epitope", "pos"]).tcr_id.nunique().unstack().to_string())
    print(f"\nat ~14 s per fold that is about {len(manifest) * 14 / 3600:.1f} GPU-hours")
    print(f"wrote inputs to {OUTDIR.relative_to(ROOT)}")
    print("next: ./scripts/31_fold_decoys.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
