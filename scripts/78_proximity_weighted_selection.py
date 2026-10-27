#!/usr/bin/env python
"""Should the public epitopes be weighted by how much they resemble the private
ones when a configuration is chosen?

E86 measured that CINGVCWTV and LLWNGPMAV are structurally almost twins of
TTDPSFLGRY (distance 1.40 and 1.47) while GILGFVFTL and GLCTLVAML are the two
public cohorts least like the private split (5.33 and 5.18 mean distance). The
obvious move follows: when selecting a configuration, stop averaging the five
public epitopes equally and upweight the ones the private split actually looks
like. This is textbook covariate-shift correction.

E85 already rejected a *different* use of the same distances -- predicting the
private AUC *level* from the structurally nearest epitope's headroom fraction
failed its leave-one-epitope-out test at MAE 0.0261 against 0.0157. That says
nothing about this question. Predicting a level and predicting *which
configuration wins* are different, and only the second one is what selection
needs.

Two tests, in the order that makes the second interpretable:

  **1. The mechanism.** Over a bank of configurations, correlate each pair of
  epitopes' config-score vectors. If a configuration that wins on A also wins
  on B precisely when A and B are structurally alike, the weighting has a
  reason. If config transfer is unrelated to structural distance, weighting is
  just discarding epitopes.

  **2. The payoff.** Nested leave-one-epitope-out, but with the *inner*
  selection objective weighted by proximity to the held-out epitope instead of
  equal. This is exactly the procedure that would be used against the private
  split, simulated on public folds where the answer is known. Swept over a
  temperature, so the equal-weight limit and the nearest-neighbour-only limit
  are both on the same curve rather than being separate claims.

The prior is against: E15 measured that cross-validated selection over this same
bank costs **-0.0048** against fixing the configuration in advance, because with
five folds the selection step has worse signal-to-noise than the thing being
selected. Proximity weighting *reduces the effective number of folds further*,
so it has to beat equal weighting by enough to pay for that too.

Run: .venv/bin/python scripts/78_proximity_weighted_selection.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/53_proximity_weighted_selection.md"
PRIVATE = {"CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"}
# Temperatures for w_i ∝ exp(-d_i / tau). inf is equal weight, 0 is
# nearest-neighbour-only; both limits belong on the same curve.
TAUS = (np.inf, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.0)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def exact_spearman_p(x, y) -> float:
    """Exact two-sided permutation p. E49: scipy's asymptotic p is 60x wrong
    at n=5, and n=10 here is still far from asymptotic."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    observed = abs(spearmanr(x, y).statistic)
    n = len(x)
    if n > 9:                       # 10! = 3.6M, enumerate with a cap
        rng = np.random.default_rng(0)
        draws = 200_000
        count = sum(abs(spearmanr(x, rng.permutation(y)).statistic) >= observed - 1e-12
                    for _ in range(draws))
        return (count + 1) / (draws + 1)
    count = sum(abs(spearmanr(x, np.array(p)).statistic) >= observed - 1e-12
                for p in permutations(y))
    return count / np.math.factorial(n)


def weights_from_distance(d: np.ndarray, tau: float) -> np.ndarray:
    if np.isinf(tau):
        return np.ones_like(d) / len(d)
    if tau == 0.0:
        w = (d == d.min()).astype(float)
        return w / w.sum()
    w = np.exp(-(d - d.min()) / tau)
    return w / w.sum()


def main() -> int:
    test = load_test()
    labels = build_training_frame()
    strict = load_module("strict", "scripts/13_strict_loeo.py")
    simmap = load_module("simmap", "scripts/77_epitope_similarity_map.py")

    # ---------- the structural distance matrix, all eight ----------
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    struct = simmap.cohort_structure(test, vectors, where)
    all_names = list(struct.index)
    F = np.column_stack([
        struct["singleton"].astype(float).to_numpy(),
        np.log(struct["rows"].astype(float).to_numpy()),
        np.log1p(struct["max_family"].astype(float).to_numpy()),
        struct["close_pairs"].astype(float).to_numpy(),
        struct["nn3_iqr"].astype(float).to_numpy(),
        struct["nn3_sd"].astype(float).to_numpy()])
    Z = (F - F.mean(0)) / (F.std(0) + 1e-12)
    D = pd.DataFrame(np.linalg.norm(Z[:, None, :] - Z[None, :, :], axis=2),
                     index=all_names, columns=all_names)

    # ---------- the configuration bank ----------
    print("building the 323-configuration bank (E15's, reused) ...")
    ranks = strict.build_component_ranks(test).set_index("tcr_id")
    lab = labels[["tcr_id", "epitope_aa", "pos"]].set_index("tcr_id")
    frame = lab.join(ranks.drop(columns=["epitope_aa"]))
    frame = frame[frame["pos"].notna()].copy()
    frame["pos"] = frame["pos"].astype(int)
    configs = strict.candidate_configs([c for c in ranks.columns
                                        if c != "epitope_aa"])
    epitopes = sorted(frame["epitope_aa"].unique())
    print(f"  {len(configs)} configurations x {len(epitopes)} public epitopes")

    # per-config, per-epitope AUC -- the whole experiment runs off this matrix
    from sklearn.metrics import roc_auc_score
    A = np.zeros((len(configs), len(epitopes)))
    for j, e in enumerate(epitopes):
        rows = (frame["epitope_aa"] == e).to_numpy()
        y = frame.loc[rows, "pos"].to_numpy()
        idx = frame.index[rows]
        for i, c in enumerate(configs):
            A[i, j] = roc_auc_score(y, strict.score(frame, c, idx))
    print(f"  config x epitope AUC matrix: {A.shape}, "
          f"range {A.min():.4f}-{A.max():.4f}")

    lines = ["# E87 — Should the public epitopes be weighted by resemblance to "
             "the private ones?", "",
             "Generated by `scripts/78_proximity_weighted_selection.py`.", ""]

    # ---------- test 1: does config transfer track structural distance? ----------
    print("\n=== 1. does configuration transfer track structural distance? ===")
    print("  Spearman between the two epitopes' 323-vectors of config AUC.")
    print(f"  {'pair':<24}{'struct dist':>12}{'config transfer ρ':>20}")
    pair_d, pair_r, pair_lab = [], [], []
    for a in range(len(epitopes)):
        for b in range(a + 1, len(epitopes)):
            ea, eb = epitopes[a], epitopes[b]
            r = spearmanr(A[:, a], A[:, b]).statistic
            d = D.loc[ea, eb]
            pair_d.append(d)
            pair_r.append(r)
            pair_lab.append(f"{ea[:4]}-{eb[:4]}")
            print(f"  {ea[:4] + '-' + eb[:4]:<24}{d:12.2f}{r:20.3f}")
    rho = spearmanr(pair_d, pair_r).statistic
    p = exact_spearman_p(pair_d, pair_r)
    print(f"\n  Spearman(structural distance, config transfer) = {rho:+.3f}, "
          f"exact p = {p:.3f}  (n = {len(pair_d)} pairs)")
    print("  a NEGATIVE correlation is what would justify proximity weighting")

    lines += ["## 1. Does configuration transfer track structural distance?", "",
              "Each cell is the Spearman correlation between two epitopes' "
              f"{len(configs)}-vectors",
              "of configuration AUC: do the same configurations win on both?", "",
              "| pair | structural distance | config transfer ρ |",
              "|---|---:|---:|"]
    for lab, d_, r_ in sorted(zip(pair_lab, pair_d, pair_r), key=lambda t: t[1]):
        lines.append(f"| {lab} | {d_:.2f} | {r_:+.3f} |")
    lines += ["",
              f"Spearman(structural distance, config transfer) = **{rho:+.3f}**, "
              f"exact p = **{p:.3f}** over {len(pair_d)} pairs.",
              "A *negative* correlation is what would justify proximity weighting.",
              ""]

    # ---------- test 1b: does ANY label-free similarity predict transfer? ----------
    # Rejecting one distance metric is not the same as closing the question.
    # Every similarity below is computable on the private rows, which is the
    # constraint that makes this hard: the quantity that would obviously
    # predict configuration transfer is configuration transfer itself, and
    # that needs labels.
    print("\n=== 1b. does any label-free similarity predict config transfer? ===")
    cham = simmap.chamfer(test, vectors, where, all_names)
    cham = pd.DataFrame(cham, index=all_names, columns=all_names)
    hla = struct["hla"].to_dict()
    sing = struct["singleton"].astype(float).to_dict()
    nrow = struct["rows"].astype(float).to_dict()
    iqr = struct["nn3_iqr"].astype(float).to_dict()
    maxf = struct["max_family"].astype(float).to_dict()
    cands = {
        "structural distance (E86)": lambda a, b: D.loc[a, b],
        "repertoire Chamfer (E86)": lambda a, b: cham.loc[a, b],
        "|log cohort size ratio|": lambda a, b: abs(np.log(nrow[a] / nrow[b])),
        "|singleton % difference|": lambda a, b: abs(sing[a] - sing[b]),
        "|nn3 IQR difference|": lambda a, b: abs(iqr[a] - iqr[b]),
        "|log largest-family ratio|": lambda a, b: abs(np.log(maxf[a] / maxf[b])),
        "HLA mismatch (0/1)": lambda a, b: float(hla[a] != hla[b]),
    }
    print(f"  {'label-free similarity':<28}{'ρ with transfer':>17}{'exact p':>10}")
    cand_rows = []
    for name, fn in cands.items():
        x = [fn(epitopes[a], epitopes[b])
             for a in range(len(epitopes)) for b in range(a + 1, len(epitopes))]
        if len(set(x)) < 2:
            continue
        r_ = spearmanr(x, pair_r).statistic
        p_ = exact_spearman_p(x, pair_r)
        print(f"  {name:<28}{r_:+17.3f}{p_:10.3f}")
        cand_rows.append([name, r_, p_])
    print("  (negative ρ = more alike means better transfer, which is what")
    print("   proximity weighting assumes)")

    lines += ["### 1b. Does any label-free similarity predict transfer?", "",
              "Rejecting one metric is not closing the question. Every candidate",
              "below is computable on the private rows -- which is the binding",
              "constraint, because the quantity that would obviously predict",
              "configuration transfer is configuration transfer itself, and that",
              "needs labels.", "",
              "| label-free similarity | ρ with config transfer | exact p |",
              "|---|---:|---:|"]
    for r in cand_rows:
        lines.append(f"| {r[0]} | {r[1]:+.3f} | {r[2]:.3f} |")
    lines += ["",
              "A *negative* ρ -- more alike means better transfer -- is what",
              "proximity weighting assumes.", ""]

    # ---------- test 2: the payoff, nested LOEO with a weighted inner objective ----
    print("\n=== 2. nested LOEO with a proximity-weighted inner objective ===")
    print("  For each held-out epitope the configuration is chosen on the other")
    print("  four, weighted by exp(-d/tau) toward the held-out one. tau=inf is")
    print("  E15's equal-weight nested LOEO; tau=0 selects on the single nearest.")
    oracle = float(A.mean(axis=1).max())
    best_apriori_idx = int(A.mean(axis=1).argmax())
    print(f"\n  oracle (best config on all five)            {oracle:.4f}")
    print(f"  {'tau':>6}{'nested LOEO':>14}{'vs equal':>11}   "
          f"distinct configs chosen")
    curve = []
    equal_score = None
    for tau in TAUS:
        held = []
        chosen = []
        for j, e in enumerate(epitopes):
            others = [k for k in range(len(epitopes)) if k != j]
            d = np.array([D.loc[e, epitopes[k]] for k in others])
            w = weights_from_distance(d, tau)
            obj = A[:, others] @ w
            pick = int(obj.argmax())
            chosen.append(pick)
            held.append(A[pick, j])
        s = float(np.mean(held))
        if np.isinf(tau):
            equal_score = s
        n_distinct = len(set(chosen))
        tau_lab = "inf" if np.isinf(tau) else f"{tau:g}"
        print(f"  {tau_lab:>6}{s:14.4f}{s - equal_score:+11.4f}   {n_distinct}")
        curve.append([tau_lab, s, s - equal_score, n_distinct, list(held)])

    lines += ["## 2. The payoff: nested LOEO with a weighted inner objective", "",
              "For each held-out epitope the configuration is chosen on the other",
              "four, weighted by `exp(-d/tau)` toward the held-out one. `tau=inf`",
              "reproduces E15's equal-weight nested LOEO; `tau=0` selects on the",
              "single structurally nearest epitope.", "",
              f"Oracle (best configuration on all five): **{oracle:.4f}**", "",
              "| tau | nested LOEO | vs equal weight | distinct configs chosen |",
              "|---|---:|---:|---:|"]
    for r in curve:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {r[2]:+.4f} | {r[3]} |")
    lines.append("")

    # per-fold detail for the two limits
    print("\n  per-fold, equal weight vs nearest-neighbour-only:")
    print(f"  {'held out':<12}{'equal':>9}{'tau=0':>9}{'oracle fold':>13}")
    eq = [r for r in curve if r[0] == "inf"][0][4]
    nn = [r for r in curve if r[0] == "0"][0][4]
    for j, e in enumerate(epitopes):
        print(f"  {e:<12}{eq[j]:9.4f}{nn[j]:9.4f}{A[:, j].max():13.4f}")
    lines += ["Per fold, the two limits against that fold's own oracle:", "",
              "| held out | equal weight | tau=0 (nearest only) | fold oracle |",
              "|---|---:|---:|---:|"]
    for j, e in enumerate(epitopes):
        lines.append(f"| {e} | {eq[j]:.4f} | {nn[j]:.4f} | {A[:, j].max():.4f} |")
    lines.append("")

    # ---------- what the weights would be on the private split ----------
    print("\n=== 3. what the weights would actually be, for the private split ===")
    pub = [e for e in all_names if e not in PRIVATE]
    for tau in (np.inf, 2.0, 1.0, 0.0):
        tau_lab = "inf" if np.isinf(tau) else f"{tau:g}"
        agg = np.zeros(len(pub))
        for pe in sorted(PRIVATE):
            d = np.array([D.loc[pe, p] for p in pub])
            agg += weights_from_distance(d, tau) / len(PRIVATE)
        print(f"  tau={tau_lab:<5} " + "  ".join(
            f"{p[:4]} {w:.3f}" for p, w in zip(pub, agg)))
    lines += ["## 3. The weights this would put on the public epitopes", "",
              "Averaged over the three private epitopes.", "",
              "| tau | " + " | ".join(p[:4] for p in pub) + " |",
              "|---|" + "---:|" * len(pub)]
    for tau in (np.inf, 2.0, 1.0, 0.0):
        tau_lab = "inf" if np.isinf(tau) else f"{tau:g}"
        agg = np.zeros(len(pub))
        for pe in sorted(PRIVATE):
            d = np.array([D.loc[pe, p] for p in pub])
            agg += weights_from_distance(d, tau) / len(PRIVATE)
        lines.append(f"| {tau_lab} | " + " | ".join(f"{w:.3f}" for w in agg) + " |")
    lines.append("")

    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
