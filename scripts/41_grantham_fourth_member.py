#!/usr/bin/env python
"""Grantham as a fourth density member, inside the rank average.

E43 tested Grantham as an extra component blended *on top of* the finished
model and got -0.0022. That is a different question from the one here: the
submission's density term is itself a rank average of three spaces, and the
proposal is to make Grantham the fourth *inside* that average, before the
0.7/0.3 structure blend.

The precedent is E25/E26, which did exactly this with the whole-V-domain
alignment distance: `density4_struct03.csv`, public 0.8525 against the
submitted 0.8498, selected in all five nested folds, and still not adopted --
first because +0.0026 is below the +0.02 threshold, then because E26 found its
gain is germline confounding, and E42 confirmed that by showing framework alone
reproduces its ranking at Spearman +0.994.

**That mechanistic objection does not transfer to Grantham.** Grantham scores
CDR3, which is the junction, so a gain would come from recognition-relevant
sequence rather than from V-gene identity. If it clears the bar there is no
standing argument against it.

Candidates, all scored identically: the current three spaces, plus Grantham,
plus the V-domain alignment, in every combination that answers the question.
Each is reported standalone, as a density average, and blended with structure,
under both the public score and strict nested leave-one-epitope-out -- the
latter because choosing a member is a selection, and E15 priced selection at
0.0048.

Run: python scripts/41_grantham_fourth_member.py
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

CACHE = ROOT / "data/interim/grantham"
REPORT = ROOT / "reports/18_grantham_fourth_member.md"
SIGMA_SCEPTR = 0.1
SIGMA_TCRDIST = 20.0
SIGMA_ALIGN = 0.1
STRUCTURE_WEIGHT = 0.3
THRESHOLD = 0.02


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def rank_within(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    return frame.assign(
        r=frame.groupby("epitope_aa")["v"].rank(pct=True)).set_index("tcr_id")["r"]


def per_epitope(labels: pd.DataFrame, series: pd.Series) -> dict:
    return {e: float(roc_auc_score(b["pos"], b["tcr_id"].map(series)))
            for e, b in labels.groupby("epitope_aa")}


def mean_auc(labels, series) -> float:
    return float(np.mean(list(per_epitope(labels, series).values())))


def density_from_matrices(test, matrices: dict, sigma: float) -> pd.Series:
    out = pd.Series(index=test["tcr_id"], dtype=float)
    for epitope, (distances, ids) in matrices.items():
        d = distances.astype(float).copy()
        np.fill_diagonal(d, np.inf)
        out.loc[ids] = np.exp(-d / sigma).mean(axis=1)
    return out


def sceptr_members(test) -> dict:
    out = {}
    for space in ("default", "small"):
        cached = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        vectors = cached["vectors"] / np.linalg.norm(
            cached["vectors"], axis=1, keepdims=True)
        matrices = {}
        for epitope, cohort in test.groupby("epitope_aa"):
            ids = cohort["tcr_id"].tolist()
            block = vectors[[position[t] for t in ids]]
            matrices[epitope] = (1.0 - block @ block.T, ids)
        out[f"SCEPTR {space}"] = density_from_matrices(test, matrices, SIGMA_SCEPTR)
    return out


def tcrdist_member(test) -> pd.Series:
    matrices = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/tcrdist/{epitope}.npz",
                         allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [position[t] for t in cohort["tcr_id"]]
        d = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
        matrices[epitope] = (d, cohort["tcr_id"].tolist())
    return density_from_matrices(test, matrices, SIGMA_TCRDIST)


def vdomain_member(test) -> pd.Series:
    matrices = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cached = np.load(ROOT / f"data/interim/fullseq/{epitope}.npz",
                         allow_pickle=True)
        position = {t: i for i, t in enumerate(cached["tcr_id"])}
        rows = [position[t] for t in cohort["tcr_id"]]
        d = (cached["alpha"] + cached["beta"]).astype(float)[np.ix_(rows, rows)]
        matrices[epitope] = (d, cohort["tcr_id"].tolist())
    return density_from_matrices(test, matrices, SIGMA_ALIGN)


def grantham_member(test, hypotheses, cards) -> pd.Series:
    """CDR3 alignment under Grantham cost, cached so this is cheap to re-run."""
    CACHE.mkdir(parents=True, exist_ok=True)
    matrix = hypotheses.grantham_matrix()
    matrices = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        path = CACHE / f"{epitope}.npy"
        if path.exists():
            d = np.load(path)
        else:
            d = cards.cdr3_distance(cohort, matrix, 11, 1)
            np.save(path, d)
        matrices[epitope] = (d.astype(float), cohort["tcr_id"].tolist())
    return density_from_matrices(test, matrices, SIGMA_ALIGN)


def combine(members: dict, names: tuple, test) -> pd.Series:
    return sum(rank_within(members[n], test) for n in names) / len(names)


def nested_loeo(test, labels, members: dict, candidates: list,
                structure_rank: pd.Series) -> dict:
    """Pick the member set on four epitopes, score the fifth. The cost of the
    choice is inside the number, which is the whole point of E15."""
    held = {}
    for epitope in sorted(labels["epitope_aa"].unique()):
        inner = labels[labels["epitope_aa"] != epitope]
        best, best_score = None, -np.inf
        for names in candidates:
            density = combine(members, names, test)
            blend = (1 - STRUCTURE_WEIGHT) * density + STRUCTURE_WEIGHT * structure_rank
            score = mean_auc(inner, blend)
            if score > best_score:
                best, best_score = names, score
        density = combine(members, best, test)
        blend = (1 - STRUCTURE_WEIGHT) * density + STRUCTURE_WEIGHT * structure_rank
        outer = labels[labels["epitope_aa"] == epitope]
        held[epitope] = (float(roc_auc_score(
            outer["pos"], outer["tcr_id"].map(blend))), best)
    return held


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    cards = load_module("cards", "scripts/27_remaining_sequence_cards.py")
    hypotheses = load_module("hyp", "scripts/40_new_hypotheses.py")

    print("building density members ...")
    members = sceptr_members(test)
    members["tcrdist3"] = tcrdist_member(test)
    members["Grantham CDR3"] = grantham_member(test, hypotheses, cards)
    members["V-domain align"] = vdomain_member(test)

    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure_rank = rank_within(structure["pae_tcr_pep_min"], test)

    print("\n=== members, standalone ===")
    for name, series in members.items():
        per = per_epitope(labels, series)
        print(f"  {name:<18} {np.mean(list(per.values())):.4f}   "
              f"GLC {per['GLCTLVAML']:.3f}")

    print("\n=== within-epitope Spearman between members ===")
    names = list(members)
    corr = pd.DataFrame(index=names, columns=names, dtype=float)
    for x in names:
        for y in names:
            corr.loc[x, y] = float(np.mean([
                spearmanr(c["tcr_id"].map(members[x]),
                          c["tcr_id"].map(members[y])).statistic
                for _, c in test.groupby("epitope_aa")]))
    print(corr.round(3).to_string())

    THREE = ("SCEPTR default", "SCEPTR small", "tcrdist3")
    sets = {
        "3 members (submitted)": THREE,
        "4 members + Grantham": THREE + ("Grantham CDR3",),
        "4 members + V-domain (E26)": THREE + ("V-domain align",),
        "5 members": THREE + ("Grantham CDR3", "V-domain align"),
        "Grantham replaces tcrdist3": ("SCEPTR default", "SCEPTR small",
                                       "Grantham CDR3"),
    }

    print("\n=== public score: density alone, and blended 0.7/0.3 with structure ===")
    rows = []
    for label, combo in sets.items():
        density = combine(members, combo, test)
        blend = (1 - STRUCTURE_WEIGHT) * density + STRUCTURE_WEIGHT * structure_rank
        per = per_epitope(labels, blend)
        rows.append({"configuration": label,
                     "density only": mean_auc(labels, density),
                     "with structure": mean_auc(labels, blend), **per})
        print(f"  {label:<28} density {rows[-1]['density only']:.4f}   "
              f"blended {rows[-1]['with structure']:.4f}")
    public = pd.DataFrame(rows).set_index("configuration")
    submitted = public.loc["3 members (submitted)", "with structure"]
    public["vs submitted"] = public["with structure"] - submitted

    print("\n  per-epitope, blended:")
    print(public[[c for c in public.columns
                  if c not in ("density only", "with structure", "vs submitted")]]
          .round(4).to_string())

    print("\n=== strict nested LOEO over the five configurations ===")
    print("  The member set is chosen on four epitopes and scored on the fifth,")
    print("  so the cost of choosing is inside the estimate (E15: 0.0048).\n")
    held = nested_loeo(test, labels, members, list(sets.values()), structure_rank)
    for epitope, (score, chosen) in held.items():
        picked = [n for n in chosen]
        print(f"  {epitope:<11} {score:.4f}   chose {len(picked)} members"
              f" ({'+Grantham' if 'Grantham CDR3' in picked else 'no Grantham'})")
    nested = float(np.mean([s for s, _ in held.values()]))
    print(f"\n  nested LOEO over member choice: {nested:.4f}")

    print("\n=== a priori, no selection: fix 4-member+Grantham in advance ===")
    fixed = combine(members, sets["4 members + Grantham"], test)
    fixed_blend = (1 - STRUCTURE_WEIGHT) * fixed + STRUCTURE_WEIGHT * structure_rank
    fixed_score = mean_auc(labels, fixed_blend)
    delta = fixed_score - submitted
    print(f"  4-member+Grantham {fixed_score:.4f} against submitted {submitted:.4f}"
          f"   delta {delta:+.4f}")

    print("\n=== private-side viability ===")
    priv = test[test["split"] == "private"]
    for epitope, cohort in priv.groupby("epitope_aa"):
        values = cohort["tcr_id"].map(members["Grantham CDR3"])
        print(f"  {epitope:<11} distinct values {values.nunique():4d}/{len(cohort)}"
              f"   sd {values.std():.4g}")

    verdict = ("CLEARS the +0.02 threshold -- rebuild the entry."
               if delta > THRESHOLD else
               f"does not clear +{THRESHOLD}; the entry stands.")
    print(f"\n  {delta:+.4f} {verdict}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Grantham as a fourth density member\n\n"
        "## Members, standalone\n\n"
        + pd.Series({n: mean_auc(labels, s) for n, s in members.items()},
                    name="mean AUC").round(4).to_markdown()
        + "\n\n## Within-epitope Spearman\n\n" + corr.round(3).to_markdown()
        + "\n\n## Configurations\n\n" + public.round(4).to_markdown()
        + f"\n\n## Nested LOEO over member choice\n\n{nested:.4f}\n\n"
        f"A priori 4-member+Grantham: {fixed_score:.4f}, "
        f"delta {delta:+.4f} against the submitted {submitted:.4f}. {verdict}\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
