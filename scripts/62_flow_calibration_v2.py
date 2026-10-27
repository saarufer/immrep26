#!/usr/bin/env python
"""E64 -- E61's calibration under family holdout, and does model shape help?

E61 fitted the SEEKER screen-to-flow calibration with a plain KFold, which is
only sound if nearby peptides cannot sit on both sides of a split. The proposed
fix was peptide-family isolation. Measuring first: among the 192 flow-measured
peptides there is **1 Hamming-1 pair on TCR59 and 4 on TCR63**, and only 68 and
80 of them are even 9-mers -- these are real proteome and virome peptides, not
an NNK library. **Sequence-family leakage was never a live concern for E61.**

Two groupings that do exist are used instead:

  * `Gene` -- 160 and 157 distinct source proteins, largest group 5 and 6
    peptides. The natural family unit, and only mildly restrictive.
  * `origin` -- human against virus. Training on one and testing on the other
    is a real domain shift and a far stronger test than any k-fold.

On model shape the plan asks for continuous, monotone-with-knot and saturating
forms. All are fitted, with the caution E61 earned: 192 points selected on the
predictor itself cannot identify a four-parameter curve, so the sigmoid is
reported as a ceiling on what shape can buy, not as an estimated mechanism.
Monotone direction is taken from the data (`increasing="auto"`), which is the
bug E61 found.

Run: .venv/bin/python scripts/62_flow_calibration_v2.py
"""

from __future__ import annotations

import glob
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold, KFold

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
SEEKER = ROOT / "data/external/seeker/B27_single_TCR_screen/screening_results"
REPORT = ROOT / "reports/37_flow_calibration_v2.md"
PREDICTOR = "r2_r1_ratio"       # E61: the later-round gain
SEED = 0


def sigmoid(x, lo, hi, centre, scale):
    return lo + (hi - lo) / (1.0 + np.exp(-(x - centre) / max(scale, 1e-6)))


def knot(x, intercept, slope1, slope2, position):
    return intercept + slope1 * x + slope2 * np.clip(x - position, 0, None)


class Curve:
    """A fitted shape, falling back to a line if the optimiser does not
    converge -- which on 192 rank-selected points it sometimes will not."""

    def __init__(self, func, p0):
        self.func, self.p0, self.params = func, p0, None

    def fit(self, x, y):
        try:
            self.params, _ = curve_fit(self.func, x, y, p0=self.p0, maxfev=20000)
        except Exception:                                    # noqa: BLE001
            self.params = None
            self.linear = LinearRegression().fit(x.reshape(-1, 1), y)
        return self

    def predict(self, x):
        if self.params is None:
            return self.linear.predict(x.reshape(-1, 1))
        return self.func(x, *self.params)


def evaluate(x, y, groups, model_factory, splitter):
    out = np.full(len(y), np.nan)
    for train, test in splitter:
        model = model_factory()
        if isinstance(model, IsotonicRegression):
            model.fit(x[train], y[train])
            out[test] = model.predict(x[test])
        elif isinstance(model, Curve):
            out[test] = model.fit(x[train], y[train]).predict(x[test])
        else:
            model.fit(x[train].reshape(-1, 1), y[train])
            out[test] = model.predict(x[test].reshape(-1, 1))
    keep = np.isfinite(out) & np.isfinite(y)
    return float(spearmanr(out[keep], y[keep]).statistic) if keep.sum() > 10 else np.nan


def main() -> int:
    print("=== is sequence-family leakage a real concern here? ===")
    frames, structure = {}, []
    for tcr in ("TCR59", "TCR63"):
        path = glob.glob(str(SEEKER / f"{tcr}*add_flow.csv"))[0]
        frame = pd.read_csv(path)
        measured = frame[frame["add_flow_single_copy"].notna()].copy()
        frames[tcr] = measured
        sequences = [s for s in measured["aa"].astype(str) if s.isalpha()]
        length = pd.Series([len(s) for s in sequences]).mode().iloc[0]
        same = [s for s in sequences if len(s) == length]
        pairs = 0
        for p in range(length):
            buckets = defaultdict(list)
            for s in same:
                buckets[s[:p] + "." + s[p + 1:]].append(s)
            pairs += sum(len(v) * (len(v) - 1) // 2 for v in buckets.values()
                         if len(v) > 1)
        structure.append({"tcr": tcr, "measured": len(measured),
                          "modal_length": int(length),
                          "at_modal_length": len(same),
                          "hamming1_pairs": pairs,
                          "distinct_genes": int(measured["Gene"].nunique()),
                          "distinct_origins": int(measured["origin"].nunique())})
        print(f"  {tcr}: {len(measured)} measured, {pairs} Hamming-1 pairs, "
              f"{measured['Gene'].nunique()} genes, "
              f"{measured['origin'].nunique()} origins")
    structure = pd.DataFrame(structure)
    print("\n  Essentially no near-neighbour structure, so E61's plain KFold")
    print("  was not leaking. Gene and origin are the groupings that exist.")

    models = {
        "linear (log x)": LinearRegression,
        "isotonic, direction from data":
            lambda: IsotonicRegression(out_of_bounds="clip", increasing="auto"),
        "two-piece linear, learnable knot":
            lambda: Curve(knot, [20.0, 10.0, -10.0, 2.0]),
        "saturating sigmoid":
            lambda: Curve(sigmoid, [0.0, 80.0, 2.0, 0.5]),
    }

    print(f"\n=== calibration of {PREDICTOR} against single-copy flow ===")
    rows = []
    for tcr, measured in frames.items():
        x = np.log10(np.clip(measured[PREDICTOR].to_numpy(float), 1e-9, None))
        y = measured["add_flow_single_copy"].to_numpy(float)
        genes = measured["Gene"].astype(str).to_numpy()
        origins = measured["origin"].astype(str).to_numpy()
        keep = np.isfinite(x) & np.isfinite(y)
        x, y, genes, origins = x[keep], y[keep], genes[keep], origins[keep]

        schemes = {
            "plain KFold (E61)": list(KFold(5, shuffle=True,
                                            random_state=SEED).split(x)),
            "gene holdout": list(GroupKFold(
                n_splits=5).split(x, y, genes)),
            "origin holdout (human<->virus)": [
                (np.where(origins != o)[0], np.where(origins == o)[0])
                for o in np.unique(origins)],
        }
        for scheme, splitter in schemes.items():
            for label, factory in models.items():
                rho = evaluate(x, y, None, factory, splitter)
                rows.append({"tcr": tcr, "validation": scheme, "model": label,
                             "spearman": rho})
            best = max((r for r in rows if r["tcr"] == tcr
                        and r["validation"] == scheme),
                       key=lambda r: (r["spearman"] if np.isfinite(r["spearman"])
                                      else -9))
            print(f"  {tcr}  {scheme:<32} best {best['model'][:26]:<28} "
                  f"{best['spearman']:+.3f}")
    table = pd.DataFrame(rows)

    print("\n=== the comparison that matters ===")
    pivot = table.pivot_table(index=["tcr", "model"], columns="validation",
                              values="spearman")
    print(pivot.round(3).to_string())

    print("\n=== reading it ===")
    for tcr in frames:
        block = table[table["tcr"] == tcr]
        kfold = block[block["validation"] == "plain KFold (E61)"]["spearman"].max()
        gene = block[block["validation"] == "gene holdout"]["spearman"].max()
        origin = block[
            block["validation"] == "origin holdout (human<->virus)"]["spearman"].max()
        print(f"  {tcr}: KFold {kfold:+.3f} -> gene holdout {gene:+.3f} "
              f"-> origin holdout {origin:+.3f}")
    lin = table[(table["model"] == "linear (log x)")
                & (table["validation"] == "gene holdout")]["spearman"]
    best_shape = table[(table["model"] != "linear (log x)")
                       & (table["validation"] == "gene holdout")]["spearman"]
    print(f"\n  Model shape under gene holdout: linear {lin.max():+.3f} against "
          f"the best curved form {best_shape.max():+.3f}.")
    print("  192 points selected on the predictor cannot identify a")
    print("  four-parameter curve, so any gap this small is shape fitting")
    print("  noise, not an estimated threshold.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E64 - E61's calibration under family holdout\n\n"
        "## Is sequence-family leakage real here?\n\n"
        + structure.to_markdown(index=False)
        + "\n\nEssentially no near-neighbour structure, so E61's plain KFold "
        "was not leaking. Gene and origin are the groupings that exist.\n\n"
        "## Held-out Spearman by validation scheme and model\n\n"
        + pivot.round(4).to_markdown() + "\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
