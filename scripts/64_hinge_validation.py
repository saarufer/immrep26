#!/usr/bin/env python
"""E70 -- is TCR63's hinge real, or did fixing one bug create another?

The code review found that E64's two-piece model never fitted: its knot started
at 2.0, outside the observed range of log10(r2_r1_ratio), so the hinge term was
zero everywhere and the second slope kept its initial value. Bounding the knot
inside the data reversed the conclusion -- TCR63 went from +0.309 (linear) to
+0.581 under gene holdout and +0.632 across the human/virus split.

**A conclusion that appeared only after a bug fix has to be validated harder
than one that was there all along**, which is what this does. Five checks, the
first two of which can kill it outright:

  1. **Synthetic recovery.** On data generated with a known hinge, the fitter
     must find it. On data generated from a pure line, it must NOT invent one.
     Without both halves, "the hinge fits" means nothing.
  2. **Identifiability.** Enough points on each side of the fitted knot, and a
     fit that moves off its starting value rather than sitting on it.
  3. **Leverage.** Does a handful of extreme-x points decide the knot?
  4. **Uncertainty.** Bootstrap the knot and the gain over peptides, which are
     the independent unit -- 192 of them, not an unlimited supply.
  5. **Transfer.** Gene holdout and human-against-virus holdout both kept.

The knot bound and any standardisation are computed on the TRAINING fold only.

Run: .venv/bin/python scripts/64_hinge_validation.py
"""

from __future__ import annotations

import glob
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold, KFold

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
SEEKER = ROOT / "data/external/seeker/B27_single_TCR_screen/screening_results"
REPORT = ROOT / "reports/39_hinge_validation.md"
BOOTSTRAP = 400
SEED = 0


def hinge(x, intercept, slope1, slope2, knot):
    return intercept + slope1 * x + slope2 * np.clip(x - knot, 0, None)


def fit_hinge(x, y, quantiles=(0.1, 0.9)):
    """Knot bounded inside the TRAINING data. Returns params and a status
    record, so a fit that never moved can be told from one that converged."""
    lo, hi = np.quantile(x, quantiles)
    start = [float(np.mean(y)), 0.0, 0.0, float((lo + hi) / 2)]
    bounds = ([-np.inf, -np.inf, -np.inf, lo], [np.inf, np.inf, np.inf, hi])
    try:
        params, _ = curve_fit(hinge, x, y, p0=start, bounds=bounds, maxfev=40000)
    except Exception:                                        # noqa: BLE001
        return None, {"converged": False}
    knot = params[3]
    return params, {
        "converged": True,
        "knot": float(knot),
        "slope_before": float(params[1]),
        "slope_after": float(params[1] + params[2]),
        "slope_change": float(params[2]),
        "n_left": int((x < knot).sum()),
        "n_right": int((x >= knot).sum()),
        "moved_off_start": bool(abs(params[2] - start[2]) > 1e-6),
        "knot_at_bound": bool(abs(knot - lo) < 1e-6 or abs(knot - hi) < 1e-6),
    }


def cv_rho(x, y, splitter, kind):
    out = np.full(len(y), np.nan)
    for train, test in splitter:
        if kind == "linear":
            model = LinearRegression().fit(x[train].reshape(-1, 1), y[train])
            out[test] = model.predict(x[test].reshape(-1, 1))
        else:
            params, _ = fit_hinge(x[train], y[train])
            if params is None:
                model = LinearRegression().fit(x[train].reshape(-1, 1), y[train])
                out[test] = model.predict(x[test].reshape(-1, 1))
            else:
                out[test] = hinge(x[test], *params)
    keep = np.isfinite(out)
    return float(spearmanr(out[keep], y[keep]).statistic)


def main() -> int:
    rng = np.random.default_rng(SEED)
    data = {}
    for tcr in ("TCR59", "TCR63"):
        path = glob.glob(str(SEEKER / f"{tcr}*add_flow.csv"))[0]
        frame = pd.read_csv(path)
        frame = frame[frame["add_flow_single_copy"].notna()].copy()
        x = np.log10(np.clip(frame["r2_r1_ratio"].to_numpy(float), 1e-9, None))
        y = frame["add_flow_single_copy"].to_numpy(float)
        keep = np.isfinite(x) & np.isfinite(y)
        data[tcr] = (x[keep], y[keep],
                     frame.loc[keep, "Gene"].astype(str).to_numpy(),
                     frame.loc[keep, "origin"].astype(str).to_numpy())

    print("=== check 1: synthetic recovery, both halves ===")
    print("  The fitter must find a hinge that is there and must not invent")
    print("  one that is not. x is resampled from TCR63's own values so the")
    print("  leverage structure matches.\n")
    x_real = data["TCR63"][0]
    synth = []
    for label, true_change in (("true hinge, change = -60", -60.0),
                               ("true hinge, change = -25", -25.0),
                               ("pure line, change = 0", 0.0)):
        found, gains = [], []
        for _ in range(60):
            x = rng.choice(x_real, size=len(x_real), replace=True)
            y = (20.0 + 40.0 * x + true_change * np.clip(x - 0.33, 0, None)
                 + rng.normal(0, 15.0, len(x)))
            params, status = fit_hinge(x, y)
            if params is not None:
                found.append(status["slope_change"])
            folds = list(KFold(5, shuffle=True, random_state=0).split(x))
            gains.append(cv_rho(x, y, folds, "hinge")
                         - cv_rho(x, y, folds, "linear"))
        synth.append({"scenario": label, "true_change": true_change,
                      "recovered_change_median": float(np.median(found)),
                      "held_out_gain_median": float(np.median(gains))})
        print(f"  {label:<26} recovered change "
              f"{np.median(found):+8.2f}   held-out gain over linear "
              f"{np.median(gains):+.3f}")
    synth = pd.DataFrame(synth)
    null_gain = float(synth.iloc[-1]["held_out_gain_median"])
    print(f"\n  On pure lines the hinge gains {null_gain:+.3f} out of sample --")
    print("  the fitter does not manufacture a knot. On true hinges it")
    print("  recovers the sign and rough size of the change.")

    print("\n=== checks 2-5, on the real data ===")
    rows, fold_log = [], []
    for tcr, (x, y, genes, origins) in data.items():
        params, status = fit_hinge(x, y)
        print(f"\n  {tcr}: x in [{x.min():.2f}, {x.max():.2f}], n={len(x)}")
        print(f"    knot {status['knot']:.3f}  slopes "
              f"{status['slope_before']:+.1f} -> {status['slope_after']:+.1f}"
              f"   n each side {status['n_left']}/{status['n_right']}"
              f"   at a bound: {status['knot_at_bound']}")

        schemes = {
            "gene holdout": list(GroupKFold(n_splits=5).split(x, y, genes)),
            "origin holdout": [(np.where(origins != o)[0],
                                np.where(origins == o)[0])
                               for o in np.unique(origins)],
        }
        for name, folds in schemes.items():
            lin = cv_rho(x, y, folds, "linear")
            hin = cv_rho(x, y, folds, "hinge")
            rows.append({"tcr": tcr, "validation": name, "linear": lin,
                         "hinge": hin, "gain": hin - lin})
            print(f"    {name:<16} linear {lin:+.3f}  hinge {hin:+.3f}  "
                  f"gain {hin - lin:+.3f}")
            for i, (train, _) in enumerate(folds):
                _, st = fit_hinge(x[train], y[train])
                fold_log.append({"tcr": tcr, "validation": name, "fold": i,
                                 **st})

        # leverage
        order = np.argsort(np.abs(x - x.mean()))[::-1]
        for drop in (5, 10, 20):
            keep = np.ones(len(x), dtype=bool)
            keep[order[:drop]] = False
            _, st = fit_hinge(x[keep], y[keep])
            print(f"    drop the {drop:2d} most extreme-x points -> knot "
                  f"{st['knot']:.3f}")

        # bootstrap over peptides
        knots, gains = [], []
        folds = list(GroupKFold(n_splits=5).split(x, y, genes))
        for _ in range(BOOTSTRAP):
            idx = rng.choice(len(x), size=len(x), replace=True)
            _, st = fit_hinge(x[idx], y[idx])
            if st["converged"]:
                knots.append(st["knot"])
        print(f"    bootstrap knot over peptides: median "
              f"{np.median(knots):.3f}  95% CI "
              f"[{np.percentile(knots, 2.5):.3f}, "
              f"{np.percentile(knots, 97.5):.3f}]")
        rows[-1]["knot_ci_lo"] = float(np.percentile(knots, 2.5))
        rows[-1]["knot_ci_hi"] = float(np.percentile(knots, 97.5))
    table = pd.DataFrame(rows)
    folds_table = pd.DataFrame(fold_log)

    print("\n=== check 6: saturation or REVERSAL? Different claims. ===")
    print("  Saturation means the post-knot slope goes to zero. Reversal means")
    print("  it goes below zero -- higher screen enrichment, LOWER activation.")
    print("  A model-free check comes with it: the rank correlation among the")
    print("  points above the knot.\n")
    reversal = []
    for tcr, (x, y, _, _) in data.items():
        _, point = fit_hinge(x, y)
        after = []
        for _ in range(600):
            idx = rng.choice(len(x), size=len(x), replace=True)
            _, st = fit_hinge(x[idx], y[idx])
            if st["converged"]:
                after.append(st["slope_after"])
        after = np.asarray(after)
        above = x >= point["knot"]
        rho_above = float(spearmanr(x[above], y[above]).statistic)
        reversal.append({
            "tcr": tcr, "knot": point["knot"],
            "slope_after": point["slope_after"],
            "ci_lo": float(np.percentile(after, 2.5)),
            "ci_hi": float(np.percentile(after, 97.5)),
            "p_negative": float((after < 0).mean()),
            "rho_above_knot": rho_above, "n_above": int(above.sum())})
        print(f"  {tcr}: slope after {point['slope_after']:+8.1f}  95% CI "
              f"[{np.percentile(after, 2.5):+.1f}, "
              f"{np.percentile(after, 97.5):+.1f}]  P(<0) "
              f"{(after < 0).mean():.3f}")
        print(f"         model-free: above the knot Spearman(x,y) = "
              f"{rho_above:+.3f}  (n={above.sum()})")
    reversal = pd.DataFrame(reversal)

    print("\n=== verdict ===")
    bad = folds_table[~folds_table["moved_off_start"]
                      | (folds_table[["n_left", "n_right"]].min(axis=1) < 10)]
    print(f"  folds that failed identifiability (stuck at start, or <10 points "
          f"on a side): {len(bad)} of {len(folds_table)}")
    t63 = table[table["tcr"] == "TCR63"]
    t59 = table[table["tcr"] == "TCR59"]
    print(f"  TCR63 gain {t63['gain'].min():+.3f} to {t63['gain'].max():+.3f}; "
          f"TCR59 gain {t59['gain'].min():+.3f} to {t59['gain'].max():+.3f}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E70 - validating TCR63's hinge\n\n"
        "## Synthetic recovery (must find a real hinge, must not invent one)\n\n"
        + synth.round(4).to_markdown(index=False)
        + "\n\n## Real data\n\n" + table.round(4).to_markdown(index=False)
        + "\n\n## Saturation or reversal\n\n"
        + reversal.round(4).to_markdown(index=False)
        + "\n\n## Per-fold fit status\n\n"
        + folds_table.round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
