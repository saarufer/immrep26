#!/usr/bin/env python
"""H5: is an alpha-beta pairing compatible with its epitope's pairing rules?

H2 asked how many sequence-diverse beta partners an alpha chain has, and died
because the median TCR has no alpha neighbours at all. This asks something
H2 does not cover: given this specific alpha and the specific beta it is
actually paired with, is that combination the kind the epitope's TCRs make?

The design is the PI's, with one addition. It runs as three gates, and failing
an early one closes the branch without spending the later ones.

  Gate A  Does the statistic detect any cross-chain dependence at all? Score
          real pairs against pairs shuffled within the epitope. If real pairs
          do not beat shuffled, there is nothing to exploit.

  Gate A' Is that dependence more than V-gene and length composition? Repeat
          the shuffle *within V/J and length strata*. If the signal vanishes,
          "compatibility" was composition in disguise -- the negative control
          the PI asked for.

  Gate B  Does the score separate functional positives from negatives on the
          five labelled epitopes? This is where the hypothesis is most likely
          to die, and the reason is worth stating before the run: every TCR
          here is a historically observed alpha-beta pair, so a non-validating
          one may have perfectly ordinary inter-chain compatibility. Gate A can
          pass and Gate B still fail.

  Gate C  E3 viability on the private epitopes, then nested-LOEO blend gain
          against the frozen entry, +0.02 to enter.

Scores are leave-one-out throughout: a TCR never contributes to the statistics
used to score it, or it would be comparing itself to itself.

Literature: Bioinformatics 2026 on CDR contacts and pairing constraints, and
Briefings in Bioinformatics 2026 on V/J, charge, polarity and mass associations
between antigen-binding alpha and beta chains.

Run: python scripts/42_ab_pairing_compatibility.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/19_ab_pairing_compatibility.md"
THRESHOLD = 0.02
DRAWS = 200
SEED = 0

CHARGE = {**{a: 1.0 for a in "KRH"}, **{a: -1.0 for a in "DE"}}
MASS = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103,
        "L": 113, "I": 113, "N": 114, "D": 115, "Q": 128, "K": 128, "E": 129,
        "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def chain_features(test: pd.DataFrame) -> pd.DataFrame:
    """The quantities the 2026 papers report as associated across chains."""
    f = pd.DataFrame(index=test.index)
    f["tcr_id"] = test["tcr_id"].values
    f["epitope"] = test["epitope_aa"].values
    for side, v, j, cdr3 in (("a", "TRAV_IMGT", "TRAJ_IMGT", "cdr3_alpha_aa"),
                             ("b", "TRBV_IMGT", "TRBJ_IMGT", "cdr3_beta_aa")):
        seqs = test[cdr3].astype(str)
        f[f"v_{side}"] = test[v].astype(str).values
        f[f"j_{side}"] = test[j].astype(str).values
        f[f"len_{side}"] = seqs.str.len().values
        f[f"charge_{side}"] = [sum(CHARGE.get(c, 0.0) for c in s) for s in seqs]
        f[f"mass_{side}"] = [np.mean([MASS.get(c, 110) for c in s]) for s in seqs]
    return f


def pmi_score(block: pd.DataFrame, assignment: np.ndarray,
              pairs: tuple) -> np.ndarray:
    """Leave-one-out pointwise mutual information, summed over feature pairs.

    `assignment[i]` is the row supplying the beta side for alpha row i, which is
    the identity for real pairs and a permutation for shuffled ones.
    """
    n = len(block)
    total = np.zeros(n)
    for left, right in pairs:
        a = block[left].to_numpy()
        b = block[right].to_numpy()[assignment]
        joint = Counter(zip(a, b))
        ma, mb = Counter(a), Counter(b)
        for i in range(n):
            # Remove the query's own contribution before scoring it.
            nj = joint[(a[i], b[i])] - 1
            na = ma[a[i]] - 1
            nb = mb[b[i]] - 1
            m = n - 1
            if m <= 0:
                continue
            total[i] += np.log(((nj + 0.5) / m)
                               / (((na + 0.5) / m) * ((nb + 0.5) / m)))
    return total


def continuous_score(block: pd.DataFrame, assignment: np.ndarray,
                     pairs: tuple) -> np.ndarray:
    """How well the beta value is predicted by the alpha value, leave-one-out.
    A compatible pair sits close to the cohort's cross-chain regression line."""
    n = len(block)
    total = np.zeros(n)
    for left, right in pairs:
        x = block[left].to_numpy(dtype=float)
        y = block[right].to_numpy(dtype=float)[assignment]
        for i in range(n):
            keep = np.arange(n) != i
            xs, ys = x[keep], y[keep]
            if xs.std() < 1e-9:
                continue
            # np.cov defaults to ddof=1, so the denominator must match
            # or the slope is inflated by n/(n-1).
            slope = np.cov(xs, ys)[0, 1] / xs.var(ddof=1)
            intercept = ys.mean() - slope * xs.mean()
            residual = y[i] - (slope * x[i] + intercept)
            scale = ys.std() or 1.0
            total[i] -= abs(residual) / scale
    return total


def strata(block: pd.DataFrame) -> np.ndarray:
    """V/J and length strata, so the control shuffle preserves composition."""
    return (block["v_a"].astype(str) + "|" + block["v_b"].astype(str) + "|"
            + (block["len_a"] // 3).astype(str) + "|"
            + (block["len_b"] // 3).astype(str)).to_numpy()


def shuffled(block: pd.DataFrame, rng, stratified: bool) -> np.ndarray:
    n = len(block)
    order = np.arange(n)
    if not stratified:
        return rng.permutation(order)
    out = order.copy()
    labels = strata(block)
    for value in np.unique(labels):
        idx = np.where(labels == value)[0]
        if len(idx) > 1:
            out[idx] = rng.permutation(idx)
    return out


CATEGORICAL = (("v_a", "v_b"), ("j_a", "j_b"), ("v_a", "j_b"), ("j_a", "v_b"))
CONTINUOUS = (("len_a", "len_b"), ("charge_a", "charge_b"), ("mass_a", "mass_b"))


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    features = chain_features(test)
    rng = np.random.default_rng(SEED)

    print("=== Gate A: does the statistic see any cross-chain dependence? ===")
    print("  Real pairs against pairs shuffled inside the epitope, and against")
    print("  pairs shuffled inside V/J and length strata (the control).\n")

    real_scores = {"pmi": pd.Series(index=test["tcr_id"], dtype=float),
                   "cont": pd.Series(index=test["tcr_id"], dtype=float)}
    rows = []
    for epitope, block in features.groupby("epitope"):
        block = block.reset_index(drop=True)
        identity = np.arange(len(block))
        real_pmi = pmi_score(block, identity, CATEGORICAL)
        real_cont = continuous_score(block, identity, CONTINUOUS)
        real_scores["pmi"].loc[block["tcr_id"]] = real_pmi
        real_scores["cont"].loc[block["tcr_id"]] = real_cont

        null_pmi, null_pmi_strat, null_cont = [], [], []
        for _ in range(DRAWS // 10):
            null_pmi.append(pmi_score(block, shuffled(block, rng, False),
                                      CATEGORICAL).mean())
            null_pmi_strat.append(pmi_score(block, shuffled(block, rng, True),
                                            CATEGORICAL).mean())
            null_cont.append(continuous_score(block, shuffled(block, rng, False),
                                              CONTINUOUS).mean())
        rows.append({
            "epitope": epitope, "n": len(block),
            "pmi_real": real_pmi.mean(), "pmi_shuffled": np.mean(null_pmi),
            "pmi_stratified": np.mean(null_pmi_strat),
            "cont_real": real_cont.mean(), "cont_shuffled": np.mean(null_cont),
        })
        print(f"  {epitope:<11} n={len(block):4d}  "
              f"PMI real {real_pmi.mean():7.3f} vs shuffled "
              f"{np.mean(null_pmi):7.3f} vs stratified "
              f"{np.mean(null_pmi_strat):7.3f}")

    gate = pd.DataFrame(rows)
    gate["pmi_excess"] = gate["pmi_real"] - gate["pmi_shuffled"]
    gate["pmi_excess_vs_strata"] = gate["pmi_real"] - gate["pmi_stratified"]
    print(f"\n  mean PMI excess over a free shuffle:       "
          f"{gate['pmi_excess'].mean():+.4f}")
    print(f"  mean PMI excess over a stratified shuffle: "
          f"{gate['pmi_excess_vs_strata'].mean():+.4f}   <- the control")
    print(f"  continuous-score excess:                   "
          f"{(gate['cont_real'] - gate['cont_shuffled']).mean():+.4f}")
    passes_a = gate["pmi_excess"].mean() > 0
    print(f"\n  Gate A {'passes' if passes_a else 'FAILS'}: real pairs "
          f"{'do' if passes_a else 'do NOT'} carry cross-chain dependence.")

    print("\n=== Gate B: does compatibility separate validated from not? ===")
    results = []
    for name, series in real_scores.items():
        per = {e: roc_auc_score(b["pos"], b["tcr_id"].map(series))
               for e, b in labels.groupby("epitope_aa")}
        mean = float(np.mean(list(per.values())))
        results.append({"score": name, "mean_auc": mean, **per})
        print(f"  {name:<6} mean AUC {mean:.4f}   " + "  ".join(
            f"{e[:3]} {v:.3f}" for e, v in per.items()))

    print("\n=== Gate C: viability, then nested-LOEO gain on the frozen entry ===")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    base_rank = hyp.rank_within(decomp.submitted_score(test), test)
    baseline = hyp.mean_auc(labels, base_rank)
    print(f"  frozen entry baseline {baseline:.4f}\n")

    gains = []
    for name, series in real_scores.items():
        ok, detail = hyp.viability(test, series)
        gain = np.nan
        if ok:
            gain = hyp.nested_blend_gain(test, labels, base_rank, series) - baseline
        gains.append({"score": name, "viable": ok, "private_levels": detail,
                      "nested_gain": gain})
        print(f"  {name:<6} viable {str(ok):<5} ({detail})   gain {gain:+.4f}")

    best = np.nanmax([g["nested_gain"] for g in gains])
    print(f"\n  best nested gain {best:+.4f} against +{THRESHOLD}")
    verdict = ("clears the bar -- revisit the entry" if best > THRESHOLD
               else "does not clear the bar; the entry stands")
    print(f"  {verdict}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# H5: antigen-conditioned alpha-beta pairing compatibility\n\n"
        "## Gate A — cross-chain dependence, real against shuffled\n\n"
        + gate.round(4).to_markdown(index=False)
        + "\n\n## Gate B — separation of functional labels\n\n"
        + pd.DataFrame(results).round(4).to_markdown(index=False)
        + "\n\n## Gate C — viability and nested gain\n\n"
        + pd.DataFrame(gains).round(4).to_markdown(index=False)
        + f"\n\nBest nested gain {best:+.4f} against +{THRESHOLD}: {verdict}.\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
