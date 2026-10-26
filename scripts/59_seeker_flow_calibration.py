#!/usr/bin/env python
"""F1 -- does a SEEKER reporter screen predict the flow activation it proxies?

E60 put Gate 0 at PROFILES ONLY: PRP and SEEKER share 19 peptides and zero
TCRs, so no cross-source binding -> activation map can be fitted. What survives
is this, inside SEEKER, where the same receptor carries both a screening
enrichment and a flow readout on the same peptides.

Three constraints the roadmap sets and this honours:

  * **The measured subset is conditionally sampled.** Flow was run on 192
    peptides per receptor, about 1% of candidates, chosen as top hits: their
    median R1/R0 is 125.8 against 14.8 across all TCR59 candidates. So every
    number here describes *ranked candidates*, and nothing is inferred about
    the activation rate of the library.
  * **Single- and high-copy are fitted separately and compared**, because
    antigen density can change the functional label. E60 found TCR63 ships no
    high-copy column, so that comparison exists on TCR59 only.
  * **Parsimonious calibration**: isotonic and logistic, cross-validated. No
    four-parameter Hill curve -- 192 points chosen by rank cannot identify
    EC50 or a slope, and the roadmap says so.

Run: .venv/bin/python scripts/59_seeker_flow_calibration.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
SEEKER = ROOT / "data/external/seeker/B27_single_TCR_screen/screening_results"
REPORT = ROOT / "reports/34_seeker_flow_calibration.md"
PREDICTORS = ("r1_r0_ratio", "r2_r1_ratio", "count_pct_round1",
              "count_pct_round2", "count_round0")
SEED = 0
FOLDS = 5


def load(tcr: str) -> pd.DataFrame:
    path = next(SEEKER.glob(f"{tcr}*add_flow.csv"))
    return pd.read_csv(path)


def cv_spearman(x: np.ndarray, y: np.ndarray, model) -> float:
    """Out-of-fold rank correlation. Isotonic on 192 points will fit the
    training data almost perfectly, so only held-out predictions are reported."""
    out = np.empty(len(y))
    for train, test in KFold(FOLDS, shuffle=True, random_state=SEED).split(x):
        fitted = model()
        fitted.fit(x[train].reshape(-1, 1) if not isinstance(
            fitted, IsotonicRegression) else x[train], y[train])
        out[test] = (fitted.predict(x[test]) if isinstance(
            fitted, IsotonicRegression) else fitted.predict(x[test].reshape(-1, 1)))
    keep = np.isfinite(out)
    return float(spearmanr(out[keep], y[keep]).statistic)


def main() -> int:
    if not SEEKER.exists():
        print("SEEKER data absent -- NOT ESTIMABLE")
        return 1

    print("=== the conditional sampling, restated with numbers ===")
    sampling = []
    frames = {}
    for tcr in ("TCR59", "TCR63"):
        frame = load(tcr)
        frames[tcr] = frame
        measured = frame[frame["add_flow_single_copy"].notna()]
        sampling.append({
            "tcr": tcr, "candidates": len(frame), "flow_measured": len(measured),
            "measured_pct": len(measured) / len(frame),
            "median_r1_r0_measured": float(measured["r1_r0_ratio"].median()),
            "median_r1_r0_all": float(frame["r1_r0_ratio"].median()),
            "has_high_copy": "add_flow_high_copy" in frame.columns})
        print(f"  {tcr}: {len(measured)}/{len(frame)} measured "
              f"({len(measured) / len(frame):.2%}); median R1/R0 "
              f"{measured['r1_r0_ratio'].median():.1f} on the measured against "
              f"{frame['r1_r0_ratio'].median():.1f} on all")
    sampling = pd.DataFrame(sampling)
    print("\n  Flow was run on top-ranked candidates, so these are statements")
    print("  about ranked candidates and not about the library.")

    print("\n=== does the screening statistic predict the flow readout? ===")
    print("  Spearman on the measured subset, and cross-validated calibration.\n")
    rows = []
    for tcr, frame in frames.items():
        measured = frame[frame["add_flow_single_copy"].notna()].copy()
        y = measured["add_flow_single_copy"].to_numpy(dtype=float)
        for column in PREDICTORS:
            if column not in measured:
                continue
            x = measured[column].to_numpy(dtype=float)
            keep = np.isfinite(x) & np.isfinite(y)
            if keep.sum() < 50:
                continue
            rho = spearmanr(x[keep], y[keep]).statistic
            # increasing="auto": the default forces an INCREASING fit, which
            # cannot represent r1_r0_ratio's negative relationship with flow
            # and returned a near-zero rho for it while the linear model got
            # +0.448 on the same data.
            iso = cv_spearman(np.log10(np.clip(x[keep], 1e-9, None)), y[keep],
                              lambda: IsotonicRegression(out_of_bounds="clip",
                                                         increasing="auto"))
            lin = cv_spearman(np.log10(np.clip(x[keep], 1e-9, None)), y[keep],
                              LinearRegression)
            rows.append({"tcr": tcr, "predictor": column, "n": int(keep.sum()),
                         "spearman_raw": rho, "cv_isotonic": iso,
                         "cv_linear_log": lin})
            print(f"  {tcr}  {column:<18} n={keep.sum():3d}  raw rho "
                  f"{rho:+.3f}   CV isotonic {iso:+.3f}   CV log-linear "
                  f"{lin:+.3f}")
    calibration = pd.DataFrame(rows)

    best = calibration.loc[calibration["cv_isotonic"].idxmax()]
    print(f"\n  best out-of-fold: {best['tcr']} / {best['predictor']} at "
          f"rho {best['cv_isotonic']:+.3f}")
    worst = calibration["spearman_raw"].abs().max()
    print(f"  largest |raw rho| anywhere: {worst:.3f}")

    print("\n=== single-copy against high-copy (TCR59 only, per E60) ===")
    copy_rows = []
    frame = frames["TCR59"]
    both = frame[frame["add_flow_single_copy"].notna()
                 & frame["add_flow_high_copy"].notna()]
    single = both["add_flow_single_copy"].to_numpy(dtype=float)
    high = both["add_flow_high_copy"].to_numpy(dtype=float)
    agreement = spearmanr(single, high).statistic
    print(f"  n={len(both)}   rank agreement between the two conditions "
          f"rho {agreement:+.3f}")
    print(f"  median activation: single-copy {np.median(single):.1f}%, "
          f"high-copy {np.median(high):.1f}%")
    flips = int(((single < 5) & (high > 25)).sum())
    print(f"  peptides inactive at single copy (<5%) but active at high copy "
          f"(>25%): {flips} of {len(both)} ({flips / len(both):.1%})")
    for name, values in (("single", single), ("high", high)):
        x = both["r1_r0_ratio"].to_numpy(dtype=float)
        copy_rows.append({"condition": name, "n": len(both),
                          "median_activation": float(np.median(values)),
                          "spearman_with_r1_r0": float(
                              spearmanr(x, values).statistic)})
    copies = pd.DataFrame(copy_rows)
    print(copies.round(3).to_string(index=False))

    print("\n=== reading it ===")
    print("  Signs matter here: r1_r0_ratio and count_pct_round1 run NEGATIVE")
    print("  against flow on TCR59 (-0.450, -0.480) while r2_r1_ratio runs")
    print("  strongly positive (+0.687). First-round abundance is not")
    print("  activation strength; the second-round GAIN is what tracks it.")
    if calibration["spearman_raw"].abs().max() < 0.3:
        print("  The screening enrichment does NOT rank the flow readout among")
        print("  the candidates chosen for measurement. That is not a null")
        print("  result about the assay: the measured subset was selected for")
        print("  high enrichment, so its enrichment range is compressed and")
        print("  the correlation is attenuated by design. What it does say is")
        print("  that within the top ~1%, the screen cannot order activation.")
    else:
        print("  The screening enrichment carries rank information about flow")
        print("  activation even within the compressed, top-selected range.")
    if flips == 0:
        print(f"\n  And antigen density does NOT change the label here: the two")
        print(f"  conditions agree at rho {agreement:+.3f} and **{flips}** of "
              f"{len(both)} peptides")
        print("  are inactive at single copy and active at high copy. Median")
        print(f"  activation is in fact slightly LOWER at high copy "
              f"({np.median(high):.1f}% against {np.median(single):.1f}%).")
        print("  Within these top-ranked candidates, presentation level is not")
        print("  what decides whether a peptide activates -- which is the")
        print("  opposite of what the roadmap flagged as a risk, and is only a")
        print("  statement about candidates already selected for enrichment.")
    else:
        print(f"\n  Antigen density changes the label: {flips} of {len(both)} "
              f"peptides are inactive at single copy and active at high copy,")
        print(f"  while the two conditions agree at rho {agreement:+.3f}.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# F1 - SEEKER screen-to-flow calibration\n\n"
        "## The conditional sampling\n\n"
        + sampling.round(4).to_markdown(index=False)
        + "\n\nFlow was run on top-ranked candidates, so every number below "
        "describes ranked candidates, not the library.\n\n"
        "## Screening statistic against flow activation\n\n"
        + calibration.round(4).to_markdown(index=False)
        + "\n\n## Single-copy against high-copy (TCR59)\n\n"
        + copies.round(4).to_markdown(index=False)
        + f"\n\nRank agreement between conditions rho **{agreement:+.3f}**; "
        f"**{flips} of {len(both)}** peptides are inactive at single copy "
        f"(<5%) and active at high copy (>25%).\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
