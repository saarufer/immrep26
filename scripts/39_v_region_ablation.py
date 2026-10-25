#!/usr/bin/env python
"""Where in the V domain does E26's alignment signal sit: framework or CDR1/2?

E26 found that masking CDR3 out of the whole-V-domain alignment distance barely
changes its ranking (Spearman +0.990), and read that as germline confounding.
A reviewer pointed out that the reading is not forced: CDR1 and CDR2 are
germline-encoded *and* genuinely contact peptide-MHC, so a V-region signal can
be real recognition rather than a study/donor artefact. The two explanations
were not separated, and separating them is cheap.

The test splits the CDR3-masked region in two:

  framework   the scaffold, which contacts nothing in the interface
  loops       CDR1 + CDR2 + CDR2.5, which do contact peptide and MHC

Signal concentrated in framework is confounding, because framework similarity
is V-gene identity and nothing else. Signal in the loops is compatible with
real recognition, though it does not prove it -- the loops are still perfectly
correlated with V-gene identity, so this test can rule the confounding reading
*in* but cannot rule it out.

CDR1, CDR2 and CDR2.5 come from tcrdist3's IMGT-derived reference keyed on the
V gene, so no numbering scheme has to be invented. Coverage is 0 genes missing,
with CDR1 located in 99.5% of chains, CDR2 in 98.0% and CDR2.5 in 99.98%; a
chain whose reference loop does not appear verbatim keeps that loop unmasked
and is counted.

E26's own split was computed in a heredoc that was never saved -- the fifth
time that has cost something -- so the full and CDR3-masked variants are
recomputed here too, which also re-checks E26's numbers. The bandwidth E26 used
is lost with the heredoc, so every region is scored over a common sigma grid
and the conclusion rests on the ordering being stable, not on matching an
absolute value.

Run: python scripts/39_v_region_ablation.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

CACHE = ROOT / "data/interim/vregion"
REPORT = ROOT / "reports/16_v_region_ablation.md"
SIGMAS = (0.02, 0.05, 0.10, 0.20)
REGIONS = ("full", "no_cdr3", "loops", "framework")


def alignment_distance(sequences: list[str]) -> np.ndarray:
    """Reuse E26's distance exactly: Needleman-Wunsch, BLOSUM62, 11/1."""
    spec = importlib.util.spec_from_file_location(
        "fullseq", ROOT / "scripts/22_fullseq_alignment_distance.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("fullseq", module)
    spec.loader.exec_module(module)
    return module.cohort_distance(sequences)


def loop_table() -> dict[str, list[str]]:
    """CDR1, CDR2 and CDR2.5 per V gene, from tcrdist3's IMGT-derived reference."""
    from tcrdist.repertoire_db import RefGeneSet

    genes = RefGeneSet("alphabeta_gammadelta_db.tsv").all_genes["human"]
    # Only V genes carry CDR1/2/2.5; the reference also holds J and D genes.
    return {name: gene.cdrs_no_gaps[:3] for name, gene in genes.items()
            if getattr(gene, "region", None) == "V"
            and getattr(gene, "cdrs_no_gaps", None)}


def build_regions(test: pd.DataFrame, loops: dict) -> tuple[pd.DataFrame, dict]:
    """Four sequence variants per chain. Returns the frame and a miss tally."""
    out = {f"{region}_{chain}": [] for region in REGIONS for chain in "ab"}
    misses = {"cdr1": 0, "cdr2": 0, "cdr25": 0, "cdr3": 0, "gene": 0}

    for row in test.itertuples(index=False):
        for chain, gene_col, seq_col, cdr3_col in (
            ("a", "TRAV_IMGT", "full_seq_alpha_aa", "cdr3_alpha_aa"),
            ("b", "TRBV_IMGT", "full_seq_beta_aa", "cdr3_beta_aa"),
        ):
            full = getattr(row, seq_col)
            cdr3 = getattr(row, cdr3_col)
            gene = str(getattr(row, gene_col)).replace("_", "/")
            trio = loops.get(gene)
            if trio is None:
                misses["gene"] += 1
                trio = ["", "", ""]

            framework = full
            present = []
            for key, loop in zip(("cdr1", "cdr2", "cdr25"), trio):
                if loop and loop in framework:
                    framework = framework.replace(loop, "", 1)
                    present.append(loop)
                else:
                    misses[key] += 1
            if cdr3 in framework:
                framework = framework.replace(cdr3, "", 1)
            else:
                misses["cdr3"] += 1

            no_cdr3 = full.replace(cdr3, "", 1) if cdr3 in full else full

            out[f"full_{chain}"].append(full)
            out[f"no_cdr3_{chain}"].append(no_cdr3)
            out[f"loops_{chain}"].append("".join(present) or full[:1])
            out[f"framework_{chain}"].append(framework)

    frame = test[["tcr_id", "epitope_aa"]].copy().reset_index(drop=True)
    for key, values in out.items():
        frame[key] = values
    return frame, misses


def density(frame: pd.DataFrame, region: str, sigma: float) -> pd.Series:
    """Kernel density within each epitope cohort over alpha + beta distance."""
    scores = pd.Series(index=frame["tcr_id"], dtype=float)
    for epitope, cohort in frame.groupby("epitope_aa", sort=True):
        path = CACHE / f"{region}__{epitope}.npy"
        if path.exists():
            total = np.load(path)
        else:
            total = sum(alignment_distance(cohort[f"{region}_{c}"].tolist())
                        for c in "ab")
            CACHE.mkdir(parents=True, exist_ok=True)
            np.save(path, total)
        d = total.astype(float).copy()
        np.fill_diagonal(d, np.inf)
        scores.loc[cohort["tcr_id"]] = np.exp(-d / sigma).mean(axis=1)
    return scores


def score(labels: pd.DataFrame, series: pd.Series) -> float:
    per = []
    for _, block in labels.groupby("epitope_aa"):
        per.append(roc_auc_score(block["pos"], block["tcr_id"].map(series)))
    return float(np.mean(per))


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    frame, misses = build_regions(test, loop_table())

    n = len(test) * 2
    print("=== region construction ===")
    print(f"  chains: {n}; V genes absent from the reference: {misses['gene']}")
    for key in ("cdr1", "cdr2", "cdr25", "cdr3"):
        print(f"  {key:<6} not located verbatim: {misses[key]:4d}  "
              f"({misses[key] / n:.2%})")
    lengths = {r: int(np.mean([len(s) for s in frame[f"{r}_a"]])) for r in REGIONS}
    print(f"\n  mean alpha-chain length by region: {lengths}")

    print("\n=== mean per-epitope AUC by region and bandwidth ===")
    table = {}
    for region in REGIONS:
        row = {}
        for sigma in SIGMAS:
            row[sigma] = score(labels, density(frame, region, sigma))
        table[region] = row
        print(f"  {region:<10} " + "  ".join(
            f"s={s}: {row[s]:.4f}" for s in SIGMAS), flush=True)
    grid = pd.DataFrame(table).T
    grid.columns = [f"sigma={s}" for s in SIGMAS]
    grid["best"] = grid.max(axis=1)

    print("\n  E26 reported full 0.8379, CDR3-masked 0.8287, CDR3-only 0.7999.")
    print("  Its bandwidth is lost with the heredoc, so compare orderings.")

    best_sigma = SIGMAS[int(np.argmax([
        np.mean([table[r][s] for r in REGIONS]) for s in SIGMAS]))]
    print(f"\n=== within-epitope Spearman between regions (sigma={best_sigma}) ===")
    series = {r: density(frame, r, best_sigma) for r in REGIONS}
    correlations = pd.DataFrame(index=REGIONS, columns=REGIONS, dtype=float)
    for x in REGIONS:
        for y in REGIONS:
            per = []
            for _, cohort in frame.groupby("epitope_aa"):
                ids = cohort["tcr_id"]
                per.append(spearmanr(series[x][ids], series[y][ids]).statistic)
            correlations.loc[x, y] = float(np.mean(per))
    print(correlations.round(3).to_string())

    print("\n=== the question: framework or loops? ===")
    fw, lp, nc = (grid.loc["framework", "best"], grid.loc["loops", "best"],
                  grid.loc["no_cdr3", "best"])
    print(f"  CDR3-masked (framework + loops) {nc:.4f}")
    print(f"  framework only                  {fw:.4f}")
    print(f"  loops only (CDR1 + CDR2 + 2.5)  {lp:.4f}")
    rho_fw = correlations.loc["no_cdr3", "framework"]
    rho_lp = correlations.loc["no_cdr3", "loops"]
    print(f"\n  rank agreement with the CDR3-masked ranking:"
          f"  framework {rho_fw:+.3f}   loops {rho_lp:+.3f}")
    if fw >= lp and rho_fw > rho_lp:
        verdict = ("Framework carries it. Framework contacts nothing in the "
                   "interface, so this is germline confounding and E26's "
                   "original reading stands.")
    elif lp > fw and rho_lp > rho_fw:
        verdict = ("The loops carry it. CDR1/CDR2 do contact peptide and MHC, "
                   "so the signal is compatible with real recognition -- though "
                   "the loops are still a function of V gene, so confounding is "
                   "not excluded.")
    else:
        verdict = ("Split verdict: AUC and rank agreement disagree about which "
                   "region dominates, so this test does not separate them.")
    print(f"\n  {verdict}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Where in the V domain does the alignment signal sit?\n\n"
        "Framework contacts nothing in the interface; CDR1, CDR2 and CDR2.5 "
        "contact peptide and MHC. Splitting E26's CDR3-masked region between "
        "them tests whether its signal is germline confounding or compatible "
        "with recognition.\n\n## Mean per-epitope AUC\n\n"
        + grid.round(4).to_markdown()
        + f"\n\n## Within-epitope Spearman between regions (sigma={best_sigma})\n\n"
        + correlations.round(3).to_markdown()
        + f"\n\n## Verdict\n\n{verdict}\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
