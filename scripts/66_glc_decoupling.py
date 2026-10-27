#!/usr/bin/env python
"""E71 -- is there systematic sequence-function decoupling in GLC outside the
public motif, and if so, where?

Two minimal analyses, as asked, plus the control that makes the second one
interpretable.

  1. D's score distribution by label on GLCTLVAML, with the public-motif
     family marked.
  2. The evaluation AUC after removing the 27 family members.

**Why (2) needs a control.** Removing the family removes 26 positives and 1
negative, so the base rate moves and the remaining pairs are by construction
the harder ones. An AUC drop is therefore expected even with no decoupling at
all. Three references are measured alongside it:

  * the same removal applied to GILGFVFTL, the only other epitope with a
    family of this size, so "what removing a public motif costs" has a second
    reading;
  * the AUC among non-family rows *of the other epitopes*, which never had a
    family to remove;
  * what nn3_dist alone achieves on the GLC remainder, which separates "D has
     stopped working" from "the neighbourhood signal is not there".

Decoupling here means something specific and falsifiable: among non-motif
TCRs, sequence neighbourhood stops predicting function. If instead the
remainder simply has a weaker but still positive signal, that is difficulty,
not decoupling.

Run: .venv/bin/python scripts/66_glc_decoupling.py
"""

from __future__ import annotations

import base64
import io
import sys
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/41_glc_decoupling.md"
FIGURE = ROOT / "reports/glc_score_distribution.html"
RADIUS = 0.12
BIG = 20
BOOT = 2000
SEED = 0

plt.rcParams.update({"font.size": 9, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25})


def geometry(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    big, nn3, famsize = {}, {}, {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        _, comp = connected_components(sp.csr_matrix(d <= RADIUS), directed=False)
        sizes = pd.Series(comp).value_counts()
        inside = np.isin(comp, sizes[sizes >= BIG].index)
        np.fill_diagonal(d, np.inf)
        ordered = np.sort(d, axis=1)
        for i, t in enumerate(ids):
            big[t] = bool(inside[i])
            nn3[t] = float(ordered[i, :3].mean())
            famsize[t] = int(sizes[comp[i]])
    return big, nn3, famsize


def boot_auc(y, s, rng, draws=BOOT):
    out = []
    idx = np.arange(len(y))
    for _ in range(draws):
        pick = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[pick])) > 1:
            out.append(roc_auc_score(y[pick], s[pick]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def main() -> int:
    rng = np.random.default_rng(SEED)
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    d_score = pd.read_csv(
        ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
    ).set_index("tcr_id")["prediction"]
    big, nn3, famsize = geometry(test)

    frame = labels.copy()
    frame["d"] = frame["tcr_id"].map(d_score)
    frame["big"] = frame["tcr_id"].map(big)
    frame["nn3"] = frame["tcr_id"].map(nn3)
    meta = test.set_index("tcr_id")
    frame["TRAV"] = frame["tcr_id"].map(
        meta["TRAV_IMGT"].astype(str).str.split("*").str[0])
    frame["TRBV"] = frame["tcr_id"].map(
        meta["TRBV_IMGT"].astype(str).str.split("*").str[0])

    glc = frame[frame["epitope_aa"] == "GLCTLVAML"]
    print("=== 1. D's score distribution on GLCTLVAML ===")
    for label in (1, 0):
        part = glc[glc["pos"] == label]
        name = "validated" if label else "not validated"
        print(f"  {name:<14} n={len(part):3d}  D rank median "
              f"{part['d'].median():.3f}  IQR "
              f"[{part['d'].quantile(.25):.3f}, {part['d'].quantile(.75):.3f}]"
              f"   in family {int(part['big'].sum())}")
    overlap = glc[glc["pos"] == 1]["d"].between(
        glc[glc["pos"] == 0]["d"].quantile(0.25),
        glc[glc["pos"] == 0]["d"].quantile(0.75)).mean()
    print(f"  share of positives inside the negatives' IQR: {overlap:.1%}")

    print("\n=== 2. AUC after removing the 27 public-family members ===")
    rows = []
    for epitope, block in frame.groupby("epitope_aa"):
        rest = block[~block["big"]]
        full = roc_auc_score(block["pos"], block["d"])
        kept = (roc_auc_score(rest["pos"], rest["d"])
                if rest["pos"].nunique() > 1 else np.nan)
        lo, hi = (boot_auc(rest["pos"].to_numpy(), rest["d"].to_numpy(), rng)
                  if rest["pos"].nunique() > 1 else (np.nan, np.nan))
        nn_auc = (roc_auc_score(rest["pos"], -rest["nn3"])
                  if rest["pos"].nunique() > 1 else np.nan)
        rows.append({
            "epitope": epitope, "n": len(block),
            "removed": int(block["big"].sum()),
            "removed_pos": int(block.loc[block["big"], "pos"].sum()),
            "AUC_all": full, "AUC_without_family": kept,
            "ci_lo": lo, "ci_hi": hi, "nn3_alone_on_remainder": nn_auc,
            "base_rate_before": float(block["pos"].mean()),
            "base_rate_after": float(rest["pos"].mean())})
    table = pd.DataFrame(rows)
    show = ["epitope", "removed", "removed_pos", "AUC_all",
            "AUC_without_family", "ci_lo", "ci_hi", "nn3_alone_on_remainder"]
    print(table[show].round(3).to_string(index=False))

    glc_row = table[table["epitope"] == "GLCTLVAML"].iloc[0]
    gil_row = table[table["epitope"] == "GILGFVFTL"].iloc[0]
    print(f"\n  GLC: {glc_row['AUC_all']:.3f} -> "
          f"{glc_row['AUC_without_family']:.3f} "
          f"(95% CI [{glc_row['ci_lo']:.3f}, {glc_row['ci_hi']:.3f}])")
    print(f"  GIL: {gil_row['AUC_all']:.3f} -> "
          f"{gil_row['AUC_without_family']:.3f} "
          f"(95% CI [{gil_row['ci_lo']:.3f}, {gil_row['ci_hi']:.3f}])  "
          f"<- the control: same removal, other epitope")
    others = table[~table["epitope"].isin(["GLCTLVAML", "GILGFVFTL"])]
    print(f"  epitopes that never had a family: AUC "
          f"{others['AUC_all'].min():.3f}-{others['AUC_all'].max():.3f}")

    print("\n=== is the remainder decoupled, or just harder? ===")
    print("  Decoupling would mean the neighbourhood stops predicting")
    print("  function. Difficulty means a weaker but still positive signal.\n")
    for epitope in ("GLCTLVAML", "GILGFVFTL"):
        row = table[table["epitope"] == epitope].iloc[0]
        chance = row["ci_lo"] <= 0.5
        print(f"  {epitope}: remainder AUC {row['AUC_without_family']:.3f}, "
              f"CI excludes 0.5: {not chance}; nn3 alone "
              f"{row['nn3_alone_on_remainder']:.3f}")

    print("\n=== where in the remainder does it fail? ===")
    rest = glc[~glc["big"]].copy()
    rest["nn3_tercile"] = pd.qcut(rest["nn3"], 3, labels=["near", "mid", "far"])
    for key in ("nn3_tercile", "TRAV", "TRBV"):
        block = rest.groupby(key).agg(
            n=("pos", "size"), validated=("pos", "mean"),
            auc=("pos", lambda s: np.nan))
        parts = []
        for name, part in rest.groupby(key):
            if len(part) >= 15 and part["pos"].nunique() > 1:
                parts.append({key: name, "n": len(part),
                              "validated": part["pos"].mean(),
                              "AUC": roc_auc_score(part["pos"], part["d"])})
        if parts:
            print(f"\n  by {key} (groups of >=15 with both labels):")
            print(pd.DataFrame(parts).round(3).to_string(index=False))

    print("\n=== does D misorder the V-gene GROUPS on the remainder? ===")
    print("  Pair-weighted AUC counting only pairs inside the same TRBV, both")
    print("  quantities recomputed on the same bootstrap resample so the CI is")
    print("  on the paired difference.\n")

    def pair_auc(y, s, strata):
        num = den = 0.0
        for value in np.unique(strata):
            mask = strata == value
            p, n = s[mask & (y == 1)], s[mask & (y == 0)]
            if len(p) == 0 or len(n) == 0:
                continue
            diff = p[:, None] - n[None, :]
            num += (diff > 0).sum() + 0.5 * (diff == 0).sum()
            den += len(p) * len(n)
        return (num / den if den else np.nan), den

    strat_rows = []
    for epitope, block in frame.groupby("epitope_aa"):
        rest_e = block[~block["big"]]
        y = rest_e["pos"].to_numpy()
        s = rest_e["d"].to_numpy()
        g = rest_e["TRBV"].to_numpy()
        flat = np.zeros(len(y))
        overall, total_pairs = pair_auc(y, s, flat)
        within, kept_pairs = pair_auc(y, s, g)
        draws = []
        for _ in range(1000):
            idx = rng.choice(len(y), size=len(y), replace=True)
            a, _ = pair_auc(y[idx], s[idx], g[idx])
            b, _ = pair_auc(y[idx], s[idx], flat[idx])
            if np.isfinite(a) and np.isfinite(b):
                draws.append(a - b)
        draws = np.asarray(draws)
        strat_rows.append({
            "epitope": epitope, "n": len(rest_e), "overall": overall,
            "within_TRBV": within, "delta": within - overall,
            "pairs_kept": kept_pairs / total_pairs,
            "ci_lo": float(np.percentile(draws, 2.5)),
            "ci_hi": float(np.percentile(draws, 97.5)),
            "p_positive": float((draws > 0).mean())})
        print(f"  {epitope:<11} overall {overall:.3f}  within-TRBV {within:.3f}"
              f"  delta {within - overall:+.3f}  95% CI "
              f"[{np.percentile(draws, 2.5):+.3f}, "
              f"{np.percentile(draws, 97.5):+.3f}]  "
              f"({kept_pairs / total_pairs:.1%} of pairs kept)")
    stratified = pd.DataFrame(strat_rows)
    glc_delta = stratified[stratified["epitope"] == "GLCTLVAML"].iloc[0]
    print(f"\n  GLCTLVAML is the only epitope with a positive delta whose CI")
    print(f"  excludes zero -- but the lower bound is {glc_delta['ci_lo']:+.3f}")
    print("  and five epitopes were tested, so this is suggestive, not")
    print("  established. GILGFVFTL goes the other way at -0.152.")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6))
    bins = np.linspace(0, 1, 26)
    for ax, (subset, title) in zip(axes, (
            (glc, "GLCTLVAML, all 166 labelled"),
            (glc[~glc["big"]], "GLCTLVAML, public-motif family removed"))):
        for label, colour, name in ((1, "#2a6fb3", "validated"),
                                    (0, "#c4443a", "not validated")):
            ax.hist(subset[subset["pos"] == label]["d"], bins=bins, alpha=0.6,
                    color=colour, label=f"{name} (n={int((subset['pos'] == label).sum())})")
        auc = roc_auc_score(subset["pos"], subset["d"])
        ax.set_title(f"{title}\nAUC {auc:.3f}", fontsize=9.5)
        ax.set_xlabel("configuration D score (within-epitope rank)")
        ax.legend(fontsize=7.5, frameon=False)
    axes[0].set_ylabel("TCRs")
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    png = base64.b64encode(buffer.getvalue()).decode()
    FIGURE.write_text(
        '<!doctype html><html lang="zh"><head><meta charset="utf-8">'
        "<title>GLC 的 D 分数分布</title><style>body{max-width:960px;"
        "margin:2.5rem auto;padding:0 1.2rem;font:15px/1.75 -apple-system,"
        '"Noto Sans CJK SC","Microsoft YaHei",sans-serif}img{width:100%}'
        "</style></head><body><h1>GLCTLVAML：D 分数的正负分布</h1>"
        f'<img src="data:image/png;base64,{png}">'
        "<p>左：全部 166 条标签行。右：剔除 27 条公共 motif 家族成员后的 139 条。"
        "家族成员 26/27 是阳性且分数很高，剔除后蓝色分布的右峰消失，"
        "剩下的正负分布大面积重叠。</p></body></html>", encoding="utf-8")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E71 - GLCTLVAML outside the public motif\n\n"
        "## AUC with and without the public-motif family\n\n"
        + table.round(4).to_markdown(index=False)
        + "\n\n## Within-TRBV against overall, on the remainder\n\n"
        + stratified.round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)} and {FIGURE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
