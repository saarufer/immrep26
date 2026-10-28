#!/usr/bin/env python
"""A per-epitope kernel bandwidth, chosen without labels.

D runs three kernel densities at fixed bandwidths — SCEPTR `default` and
`small` at σ=0.1, tcrdist at σ=20 — the same for every cohort. E96 showed the
cohorts are not alike enough for that to be obviously right: the share of the
within-epitope distance explained by CDR3 length alone varies 15× across them.

E16 already learned a bandwidth by gradient descent and failed: it converged to
0.20–0.26 against a true optimum near 0.07 and scored 0.8199 against plain
cosine's 0.8396. But E16 learned **one global** bandwidth from labels. This is a
different question — a **per-epitope** bandwidth chosen from a **label-free**
property of the cohort, which is the only kind that can be applied to the
private split at all.

The scouting measurement that motivates it:

  per-epitope oracle σ    0.8459
  global σ = 0.1 (D's)    0.8396      headroom +0.0064

and, unexpectedly, the obvious scale statistics track the oracle **with the
wrong sign**: Spearman(median pairwise distance, oracle σ) = −0.821 over the
five public epitopes. A *looser* cohort wants a *smaller* bandwidth, which is
the reverse of the standard adaptive-KDE rule. The mechanism is visible once
stated: in a sparse cohort the only real signal is the handful of genuinely
close pairs, and a wide kernel washes them out; in a dense cohort a narrow
kernel throws away the family. Both are the same statement — **hold the
effective number of contributing neighbours constant** — and that is the rule
this tests, against two cruder inverse-scale rules and the oracle.

Every rule carries at most one global constant, fitted on four epitopes and
applied to the fifth, so the +0.0064 is a ceiling that includes no selection
cost and the nested figure is what would actually be realised.

Run: .venv/bin/python scripts/87_per_epitope_bandwidth.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/64_per_epitope_bandwidth.md"
SPACES = ("default", "small")
SIG_SCEPTR = np.array([0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5])
SIG_TCRDIST = np.array([4.0, 6.0, 10.0, 14.0, 20.0, 30.0, 45.0, 70.0])
EFF_GRID = np.array([1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 35.0])
# p = 0 reproduces D's fixed bandwidths exactly, so the rule may decline.
COEF_GRID = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0])


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def exact_p(x, y) -> float:
    obs = abs(stats.spearmanr(x, y).statistic)
    n = len(x)
    hits = sum(abs(stats.spearmanr(x, list(p)).statistic) >= obs - 1e-12
               for p in permutations(y))
    return hits / np.math.factorial(n)


def effective_neighbours(D: np.ndarray, sigma: float) -> float:
    """Mean over TCRs of the kernel mass each one's neighbourhood carries.

    Scale-free and label-free: it counts how many neighbours actually
    contribute at this bandwidth, which is the quantity the rule holds fixed.
    """
    K = np.exp(-D / sigma)
    K[~np.isfinite(K)] = 0.0
    return float(K.sum(axis=1).mean())


def sigma_for_eff(D: np.ndarray, target: float, grid: np.ndarray) -> float:
    """The grid bandwidth whose effective neighbour count is closest to target."""
    eff = np.array([effective_neighbours(D, s) for s in grid])
    return float(grid[int(np.argmin(np.abs(np.log(eff + 1e-12)
                                           - np.log(target))))])


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    y = scored.set_index("tcr_id")["pos"]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    epitopes = sorted(test["epitope_aa"].unique())

    # ---------- distance matrices per (space, epitope) ----------
    Dmats: dict[str, dict] = {}
    for space in SPACES:
        c = np.load(ROOT / f"data/interim/sceptr/{space}.npz", allow_pickle=True)
        W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
        po = {t: i for i, t in enumerate(c["tcr_id"])}
        Dmats[space] = {}
        for e, g in test.groupby("epitope_aa"):
            X = W[[po[t] for t in g["tcr_id"]]]
            D = 1.0 - X @ X.T
            np.fill_diagonal(D, np.inf)
            Dmats[space][e] = (D, g["tcr_id"].to_numpy())
    Dmats["tcrdist"] = {}
    for e, g in test.groupby("epitope_aa"):
        c = np.load(ROOT / f"data/interim/tcrdist/{e}.npz", allow_pickle=True)
        po = {t: i for i, t in enumerate(c["tcr_id"])}
        rows = [po[t] for t in g["tcr_id"]]
        D = (c["alpha"] + c["beta"]).astype(float)[np.ix_(rows, rows)]
        np.fill_diagonal(D, np.inf)
        Dmats["tcrdist"][e] = (D, g["tcr_id"].to_numpy())
    GRID = {"default": SIG_SCEPTR, "small": SIG_SCEPTR, "tcrdist": SIG_TCRDIST}
    FIXED = {"default": 0.1, "small": 0.1, "tcrdist": 20.0}
    pub = [e for e in epitopes if y.reindex(Dmats["default"][e][1]).notna().sum() >= 10]

    # ---------- the three σ rules, all label-free ----------
    def sigmas_fixed(space):
        return {e: FIXED[space] for e in epitopes}

    def sigmas_eff(space, target):
        return {e: sigma_for_eff(Dmats[space][e][0], target, GRID[space])
                for e in epitopes}

    scale_nn1, scale_med = {}, {}
    for sp in ("default", "small", "tcrdist"):
        scale_nn1[sp], scale_med[sp] = {}, {}
        for e in epitopes:
            D = Dmats[sp][e][0]
            iu = np.triu_indices(len(D), 1)
            scale_nn1[sp][e] = float(np.median(np.sort(D, axis=1)[:, 0]))
            scale_med[sp][e] = float(np.median(D[iu]))

    def sigmas_inverse(space, power, scale):
        """sigma_e = sigma_base[space] * (geometric-mean scale / this cohort's)^p.

        Dimensionally correct -- the exponent is the only free number and the
        component's own fixed bandwidth is the anchor, so **p = 0 reproduces D
        exactly**. An earlier draft used a bare coefficient as the bandwidth,
        which is dimensionless nonsense for tcrdist (distances in the hundreds)
        and silently saturated at the grid floor for every cohort.
        """
        sc = scale[space]
        gm = float(np.exp(np.mean(np.log(list(sc.values())))))
        out = {}
        for e in epitopes:
            s = FIXED[space] * (gm / sc[e]) ** power
            out[e] = float(GRID[space][int(np.argmin(np.abs(
                np.log(GRID[space]) - np.log(max(s, 1e-12)))))])
        return out

    def sigmas_invmed(space, power):
        return sigmas_inverse(space, power, scale_med)

    def sigmas_invnn1(space, power):
        return sigmas_inverse(space, power, scale_nn1)

    def density(space, sig_by_ep) -> pd.Series:
        out = pd.Series(index=test["tcr_id"], dtype=float)
        for e in epitopes:
            D, ids = Dmats[space][e]
            out.loc[ids] = np.exp(-D / sig_by_ep[e]).mean(axis=1)
        return out

    def macro(sig_by_space, eps=None) -> float:
        parts = [hyp.rank_within(density(sp, sig_by_space[sp]), test)
                 for sp in ("default", "small", "tcrdist")]
        s = sum(parts) / 3
        per = hyp.per_epitope(labels, s)
        use = eps if eps is not None else pub
        return float(np.mean([per[e] for e in use if e in per]))

    print("=== the three density components at D's fixed bandwidths ===")
    base_sig = {sp: sigmas_fixed(sp) for sp in ("default", "small", "tcrdist")}
    base = macro(base_sig)
    print(f"  density3 macro at fixed sigma: {base:.4f}")

    # ---------- oracle per-epitope, for the ceiling ----------
    print("\n=== per-epitope oracle sigma, per component (ceiling, uses labels) ===")
    oracle_sig = {}
    for sp in ("default", "small", "tcrdist"):
        oracle_sig[sp] = {}
        for e in epitopes:
            D, ids = Dmats[sp][e]
            yy = y.reindex(ids)
            m = yy.notna().to_numpy()
            if m.sum() < 10:
                oracle_sig[sp][e] = FIXED[sp]
                continue
            aucs = [roc_auc_score(yy[m], np.exp(-D / s).mean(axis=1)[m])
                    for s in GRID[sp]]
            oracle_sig[sp][e] = float(GRID[sp][int(np.argmax(aucs))])
        print(f"  {sp:<9} " + "  ".join(f"{e[:3]} {oracle_sig[sp][e]:g}"
                                        for e in pub))
    ceiling = macro(oracle_sig)
    print(f"  density3 macro at oracle sigma: {ceiling:.4f}  "
          f"(headroom {ceiling-base:+.4f})")

    # ---------- does the effective-neighbour rule reproduce the oracle? ----------
    print("\n=== does a label-free rule pick the oracle sigma? (component 'default') ===")
    D0 = Dmats["default"]
    med = {e: float(np.median(D0[e][0][np.triu_indices(len(D0[e][0]), 1)]))
           for e in epitopes}
    orc = [oracle_sig["default"][e] for e in pub]
    mp = [med[e] for e in pub]
    rho = stats.spearmanr(mp, orc).statistic
    print(f"  Spearman(median pairwise distance, oracle sigma) = {rho:+.3f}, "
          f"exact p = {exact_p(mp, orc):.3f}  (n={len(pub)})")
    for target in EFF_GRID[:4]:
        picked = sigmas_eff("default", target)
        r = stats.spearmanr([picked[e] for e in pub], orc).statistic
        print(f"  effective-neighbour target {target:5.1f}: picks "
              + "  ".join(f"{e[:3]} {picked[e]:g}" for e in pub)
              + f"   rho with oracle {r:+.2f}")

    # ---------- the crux: can a label-free statistic even span the oracle? ----------
    # A rule can only reproduce the oracle if its own spread across cohorts is
    # at least as large as the oracle's. Correlation is not enough.
    print("\n=== dynamic range: how far does each quantity move across cohorts? ===")
    print(f"  {'quantity':<28}{'min':>9}{'max':>9}{'ratio':>8}")
    for nm, d in (("oracle sigma (default)",
                   {e: oracle_sig["default"][e] for e in pub}),
                  ("median pairwise distance", {e: scale_med["default"][e]
                                                for e in pub}),
                  ("median nn1 distance", {e: scale_nn1["default"][e]
                                           for e in pub}),
                  ("effective neighbours @0.1",
                   {e: effective_neighbours(Dmats["default"][e][0], 0.1)
                    for e in pub})):
        v = np.array(list(d.values()))
        print(f"  {nm:<28}{v.min():9.3f}{v.max():9.3f}{v.max()/v.min():8.2f}x")

    # ---------- nested LOEO over every rule ----------
    print("\n=== nested LOEO: fit the rule's one constant on 4 epitopes, apply to the 5th ===")
    rules = {
        "fixed sigma (D's own)": ("fixed", [None]),
        "constant effective neighbours": ("eff", EFF_GRID),
        "inverse median distance": ("inv", COEF_GRID),
        "inverse median nn1": ("invnn", COEF_GRID),
    }
    print(f"  {'rule':<34}{'nested macro':>14}{'vs D-fixed':>12}"
          f"   constants chosen")
    results = []
    for name, (kind, grid) in rules.items():
        held, chosen = [], []
        for e in pub:
            tr = [x for x in pub if x != e]
            if kind == "fixed":
                sig = {sp: sigmas_fixed(sp) for sp in ("default", "small", "tcrdist")}
                best = None
            else:
                scores = []
                for cst in grid:
                    f = {"eff": sigmas_eff, "inv": sigmas_invmed,
                         "invnn": sigmas_invnn1}[kind]
                    sig_c = {sp: f(sp, cst) for sp in ("default", "small", "tcrdist")}
                    scores.append(macro(sig_c, eps=tr))
                bi = int(np.argmax(scores))
                best = grid[bi]
                f = {"eff": sigmas_eff, "inv": sigmas_invmed,
                     "invnn": sigmas_invnn1}[kind]
                sig = {sp: f(sp, best) for sp in ("default", "small", "tcrdist")}
            held.append(macro(sig, eps=[e]))
            chosen.append(best)
        nested = float(np.mean(held))
        cs = "fixed" if kind == "fixed" else \
            "/".join(f"{c:g}" for c in chosen)
        print(f"  {name:<34}{nested:14.4f}{nested-base:+12.4f}   {cs}")
        results.append([name, nested, nested - base, cs])

    # ---------- what the surviving rule picks, including on the private split ----------
    print("\n=== what the inverse-nn1 rule picks, all eight cohorts ===")
    # the constant chosen by a majority of folds above
    from collections import Counter
    coef = float(sorted(Counter(
        [float(x) for x in results[-1][3].split("/")]).items(),
        key=lambda kv: (-kv[1], kv[0]))[0][0])
    picks = {sp: sigmas_invnn1(sp, coef) for sp in ("default", "small", "tcrdist")}
    print(f"  exponent p = {coef:g} (modal across folds; p=0 would mean 'decline')")
    print(f"  {'epitope':<12}{'split':<9}{'median nn1':>12}"
          + "".join(f"{sp:>10}" for sp in ("default", "small", "tcrdist"))
          + f"{'D uses':>22}")
    pick_rows = []
    sp_of = test.drop_duplicates("epitope_aa").set_index("epitope_aa")["split"]
    for e in epitopes:
        print(f"  {e:<12}{sp_of[e]:<9}{scale_nn1['default'][e]:12.3f}"
              + "".join(f"{picks[sp][e]:10g}" for sp in
                        ("default", "small", "tcrdist"))
              + f"{'0.1 / 0.1 / 20':>22}")
        pick_rows.append([e, sp_of[e], scale_nn1["default"][e],
                          picks["default"][e], picks["small"][e],
                          picks["tcrdist"][e]])

    # ---------- put it back inside the full D ----------
    print("\n=== the rule inside the full submitted D ===")
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)

    def density3(sig_by_space) -> pd.Series:
        return sum(hyp.rank_within(density(sp, sig_by_space[sp]), test)
                   for sp in ("default", "small", "tcrdist")) / 3

    d3_old, d3_new = density3(base_sig), density3(picks)
    d_new = dr + (d3_new - d3_old) / 4
    per_old, per_new = hyp.per_epitope(labels, dr), hyp.per_epitope(labels, d_new)
    print(f"  {'epitope':<12}{'D':>9}{'D with rule':>13}{'diff':>9}")
    for e in pub:
        print(f"  {e:<12}{per_old[e]:9.4f}{per_new[e]:13.4f}"
              f"{per_new[e]-per_old[e]:+9.4f}")
    mo = float(np.mean([per_old[e] for e in pub]))
    mn = float(np.mean([per_new[e] for e in pub]))
    print(f"  {'MACRO':<12}{mo:9.4f}{mn:13.4f}{mn-mo:+9.4f}")
    # nested over {D, D with rule}
    held = []
    for e in pub:
        tr = [x for x in pub if x != e]
        pick = max((("D", per_old), ("D+rule", per_new)),
                   key=lambda kv: np.mean([kv[1][t] for t in tr]))
        held.append(pick[1][e])
    print(f"  nested LOEO over the two: {np.mean(held):.4f} "
          f"({np.mean(held)-mo:+.4f})")

    lines = ["# E97 — A per-epitope kernel bandwidth, chosen without labels", "",
             "Generated by `scripts/87_per_epitope_bandwidth.py`.", "",
             "## The ceiling", "",
             f"- density3 at D's fixed bandwidths: **{base:.4f}**",
             f"- density3 at the per-epitope oracle: **{ceiling:.4f}** "
             f"(headroom **{ceiling-base:+.4f}**)", "",
             "Per-epitope oracle bandwidths:", "",
             "| component | " + " | ".join(e for e in pub) + " |",
             "|---|" + "---:|" * len(pub)]
    for sp in ("default", "small", "tcrdist"):
        lines.append(f"| {sp} | " + " | ".join(f"{oracle_sig[sp][e]:g}"
                                               for e in pub) + " |")
    lines += ["",
              "## The sign is backwards from standard adaptive KDE", "",
              f"Spearman(median pairwise distance, oracle σ) = **{rho:+.3f}**, "
              f"exact p **{exact_p(mp, orc):.3f}** over {len(pub)} epitopes. A",
              "*looser* cohort wants a *narrower* kernel, i.e. what is really being",
              "held constant is the number of neighbours that contribute.", "",
              "## Nested LOEO", "",
              "Each rule carries at most one global constant, fitted on four",
              "epitopes and applied to the fifth.", "",
              "| rule | nested macro | vs D's fixed σ | constants chosen |",
              "|---|---:|---:|---|"]
    for r in results:
        lines.append(f"| {r[0]} | **{r[1]:.4f}** | {r[2]:+.4f} | {r[3]} |")
    lines += ["", "## Dynamic range is the whole story", "",
              "A rule can only reproduce an oracle it can span, and the three rules",
              "order by exactly that:", "",
              "| quantity | min | max | ratio |", "|---|---:|---:|---:|"]
    for nm, d in (("**oracle σ (target)**",
                   {e: oracle_sig["default"][e] for e in pub}),
                  ("median pairwise distance", {e: scale_med["default"][e]
                                                for e in pub}),
                  ("median nn1 distance", {e: scale_nn1["default"][e]
                                           for e in pub}),
                  ("effective neighbours @0.1",
                   {e: effective_neighbours(Dmats["default"][e][0], 0.1)
                    for e in pub})):
        v = np.array(list(d.values()))
        lines.append(f"| {nm} | {v.min():.3f} | {v.max():.3f} "
                     f"| **{v.max()/v.min():.2f}×** |")
    lines += ["", "## What the surviving rule picks, all eight cohorts", "",
              f"Exponent p = {coef:g}, modal across folds (p=0 would reproduce D). "
              "D uses 0.1 / 0.1 / 20.", "",
              "| epitope | split | median nn1 | default | small | tcrdist |",
              "|---|---|---:|---:|---:|---:|"]
    for r in pick_rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:.3f} | {r[3]:g} | {r[4]:g} "
                     f"| {r[5]:g} |")
    lines += ["", "## Inside the full submitted D", "",
              "| epitope | D | D with rule | diff |", "|---|---:|---:|---:|"]
    for e in pub:
        lines.append(f"| {e} | {per_old[e]:.4f} | {per_new[e]:.4f} "
                     f"| {per_new[e]-per_old[e]:+.4f} |")
    lines += [f"| **MACRO** | {mo:.4f} | {mn:.4f} | **{mn-mo:+.4f}** |", "",
              f"Nested LOEO over the two: **{np.mean(held):.4f}** "
              f"({np.mean(held)-mo:+.4f}).", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
