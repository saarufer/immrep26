#!/usr/bin/env python
"""E151 -- what the IMMREP26 private leaderboard can and cannot resolve.

The private split is three peptide cohorts of 243 / 229 / 130 rows. After the
board was published, the honest field ran from 0.846 down to 0.467 with no gap,
which invites reading rank differences that the evaluation cannot support. This
quantifies that directly, on the public rows, by bootstrapping a split shaped
like the private one.

Two quantities, and the distinction between them is the whole point:

  **marginal sd**   the spread of one team's macro-AUC under row resampling.
                    Every team is scored on the *same* 602 rows, so this wobble
                    is common to the field and does not by itself reorder it.
  **paired sd**     the spread of the *difference* between two teams. This is
                    what reorders the board, and it shrinks with how correlated
                    the two predictions are: Var(A-B) = 2 sigma^2 (1 - rho).

For two unrelated entries rho(bootstrap AUC) ~ 0, so sd(diff) = sqrt(2)*sigma
and the 95% resolvable gap is 1.96*sqrt(2)*sigma.

Also reported, since it is the quantity the C1 post-mortem turns on: the
expected margin of C1 over D on a private-shaped split, against the -0.008
actually realised.

Run: .venv/bin/python scripts/129_leaderboard_resolution.py
"""

from __future__ import annotations

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

from immrep.data import build_training_frame  # noqa: E402

REPORT = ROOT / "reports/108_leaderboard_resolution.md"

# the private split's shape, and three public cohorts standing in for it
PRIVATE_SIZES = [243, 229, 130]
COHORTS = [("YLQPRTFLL", 243), ("TTDPSFLGRY", 229), ("GLCTLVAML", 130)]
B = 4000

SUBS = {
    "D": "d_density3_m30pool_m50_struct.csv",
    "C1": "C1_consensus_S4_FHB.csv",
    "S4": "S4_D_stack9.csv",
    "F": "f_density3_m30_m50_struct_nopool.csv",
    "density": "density_sceptr_default_sigma0.1.csv",
}
PAIRS = [("C1", "S4"), ("D", "F"), ("D", "C1"), ("D", "density")]


def main() -> int:
    lab = build_training_frame()[["tcr_id", "epitope_aa", "pos"]].dropna(subset=["pos"]).copy()
    lab["pos"] = lab["pos"].astype(int)

    missing = [f for f in SUBS.values() if not (ROOT / "submissions" / f).exists()]
    if missing:
        print("submissions/ is gitignored; rebuild with `make submission` first.")
        print("missing:", ", ".join(missing))
        return 1

    for name, fn in SUBS.items():
        s = pd.read_csv(ROOT / "submissions" / fn).set_index("tcr_id")["prediction"]
        lab[name] = lab["tcr_id"].map(s)
        lab[name] = lab.groupby("epitope_aa")[name].rank(pct=True)

    y = lab["pos"].to_numpy()
    ep = lab["epitope_aa"].to_numpy()
    rng = np.random.default_rng(1)
    pools = {e: np.where(ep == e)[0] for e, _ in COHORTS}

    boot: dict[str, list[float]] = {k: [] for k in SUBS}
    for _ in range(B):
        idx = [rng.choice(pools[e], size=n, replace=True) for e, n in COHORTS]
        if any(len(set(y[i])) < 2 for i in idx):
            continue
        for k in SUBS:
            x = lab[k].to_numpy()
            boot[k].append(float(np.mean([roc_auc_score(y[i], x[i]) for i in idx])))

    A = {k: np.array(v) for k, v in boot.items()}
    n = len(A["D"])
    sigma = float(np.mean([A[k].std(ddof=1) for k in SUBS]))

    print(f"=== one entry, fixed cohorts {PRIVATE_SIZES}, row bootstrap, B={n} ===")
    for k in SUBS:
        print(f"  {k:8s} mean {A[k].mean():.4f}   sd {A[k].std(ddof=1):.4f}")
    print(f"\n  marginal sd (mean over entries)  sigma = {sigma:.4f}")

    print("\n=== difference between two entries, scored on the SAME rows ===")
    print(f"  {'pair':<18}{'rho(pred)':>10}{'rho(AUC)':>10}{'sd(diff)':>10}{'95% gap':>10}")
    rows = []
    for a, b in PAIRS:
        rp = float(np.mean([spearmanr(lab.loc[lab.epitope_aa == e, a],
                                      lab.loc[lab.epitope_aa == e, b]).statistic
                            for e, _ in COHORTS]))
        rb = float(np.corrcoef(A[a], A[b])[0, 1])
        sd = float((A[a] - A[b]).std(ddof=1))
        print(f"  {a + '-' + b:<18}{rp:+10.3f}{rb:+10.3f}{sd:10.4f}{1.96 * sd:10.3f}")
        rows.append([f"{a} vs {b}", rp, rb, sd, 1.96 * sd])

    unrel_sd = float(np.sqrt(2) * sigma)
    print(f"  {'unrelated (rho=0)':<18}{0.0:+10.3f}{0.0:+10.3f}{unrel_sd:10.4f}"
          f"{1.96 * unrel_sd:10.3f}")
    rows.append(["two unrelated entries", 0.0, 0.0, unrel_sd, 1.96 * unrel_sd])

    dc = A["D"] - A["C1"]
    print(f"\n=== C1 over D on a private-shaped split ===")
    print(f"  expected margin  {-dc.mean():+.4f} +/- {dc.std(ddof=1):.4f}")
    print(f"  realised margin  -0.0080   ({(-0.0080 - -dc.mean()) / dc.std(ddof=1):+.1f} sd)")

    lines = ["# E151 — what the private leaderboard can resolve", "",
             "Generated by `scripts/129_leaderboard_resolution.py`.", "",
             f"Public rows bootstrapped into a private-shaped split "
             f"({' / '.join(map(str, PRIVATE_SIZES))} rows over three cohorts), B={n}.", "",
             "## One entry", "", "| entry | mean | sd |", "|---|---:|---:|"]
    for k in SUBS:
        lines.append(f"| {k} | {A[k].mean():.4f} | {A[k].std(ddof=1):.4f} |")
    lines += ["", f"Marginal sd **{sigma:.4f}** — common to every team, since all are",
              "scored on the same rows, so it does not by itself reorder the board.", "",
              "## Difference between two entries", "",
              "| pair | ρ(prediction) | ρ(bootstrap AUC) | sd(difference) | gap for 95% |",
              "|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:+.3f} | {r[2]:+.3f} | {r[3]:.4f} | **{r[4]:.3f}** |")
    lines += ["", "## C1 over D", "",
              f"Expected on a private-shaped split: **{-dc.mean():+.4f} ± {dc.std(ddof=1):.4f}**. ",
              "Realised on the private split: **−0.0080**, about "
              f"{abs((-0.0080 - -dc.mean()) / dc.std(ddof=1)):.0f} sd away. The relationship",
              "between the two changed between the public and private peptides; the",
              "private split alone cannot establish the weaker claim that D beats C1",
              "there, since 0.008 is 1.3 paired sd.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
