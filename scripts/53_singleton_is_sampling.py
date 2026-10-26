#!/usr/bin/env python
"""E55 -- how much of "singleton" is just sample size?

E53 showed singleton status is substantially a property of the TCR: it is
predictable from CDR3 length, Pgen and interface confidence. That answers "which
TCRs lose their neighbours first", but not "would they have neighbours in a
larger sample", and those are different claims. The AIRR literature is explicit
that publicness tracks sample size, donor count and sequencing depth, so not
observing a neighbour is not evidence that none exists.

E21 already thinned the public cohorts, but it measured what happened to the
AUC. The statistic never extracted is the one that answers this directly:

    as the cohort shrinks, how many TCRs that HAD neighbours become singletons,
    and does that differ between validated and non-validated TCRs?

The decisive comparison is at matched size. Each public cohort is thinned to the
three private cohort sizes -- 130, 229, 243 -- and its singleton rate compared
with what the private cohorts actually show at those sizes (97.5%, 95.6%,
85.4%). If a thinned public cohort reaches the private rate, private sparsity is
sample size. If it does not, the private cohorts are sparser than sampling alone
explains.

The graph is rebuilt on the thinned **full** cohort, including unlabelled rows,
because a neighbour that happens to carry no label is still a neighbour; purity
and the positive/negative split are then read off the labelled survivors only.
This is the same discipline E54 used.

What it cannot do, stated plainly: this subsamples a finite candidate library.
It cannot reconstruct the antigen-specific repertoire of a real immune system,
so it bounds the sampling explanation from below, not from above.

Run: .venv/bin/python scripts/53_singleton_is_sampling.py
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

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/28_singleton_is_sampling.md"
FIGURE = ROOT / "reports/singleton_vs_cohort_size.html"
RADIUS = 0.12
REPEATS = 200
PRIVATE_SIZES = {"LTDEMIAQY": 130, "LLWNGPMAV": 229, "CINGVCWTV": 243}
SWEEP = (80, 100, 130, 160, 192, 229, 243, 300, 366, 387, 419, 500, 670)
SEED = 0

plt.rcParams.update({"font.size": 9, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25})


def distances(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        out[epitope] = (ids, 1.0 - block @ block.T)
    return out


def singletons(adjacency) -> np.ndarray:
    """True where a node has no edge within the radius. Single-linkage families
    of size one, which is the project's definition."""
    _, component = connected_components(sp.csr_matrix(adjacency), directed=False)
    sizes = np.bincount(component)
    return sizes[component] == 1


def main() -> int:
    rng = np.random.default_rng(SEED)
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    label_of = dict(zip(labels["tcr_id"], labels["pos"]))
    dist = distances(test)
    public = sorted(labels["epitope_aa"].unique())

    print("=== observed, at each cohort's own size ===")
    observed = []
    for epitope, (ids, d) in dist.items():
        single = singletons(d <= RADIUS)
        y = np.array([label_of.get(t, -1) for t in ids])
        row = {"epitope": epitope, "n": len(ids),
               "split": test.loc[test["epitope_aa"] == epitope, "split"].iloc[0],
               "singleton_all": float(single.mean())}
        if (y >= 0).any():
            row["singleton_pos"] = float(single[y == 1].mean())
            row["singleton_neg"] = float(single[y == 0].mean())
        observed.append(row)
    observed = pd.DataFrame(observed)
    print(observed.round(3).to_string(index=False))

    print(f"\n=== thinning each public cohort, {REPEATS} draws per size ===")
    print("  Graph rebuilt on the thinned FULL cohort; the positive/negative")
    print("  split is read off the labelled survivors.\n")
    rows = []
    for epitope in public:
        ids, d = dist[epitope]
        adjacency_full = d <= RADIUS
        base_single = singletons(adjacency_full)
        y = np.array([label_of.get(t, -1) for t in ids])
        had_neighbour = ~base_single
        for n in SWEEP:
            if n > len(ids):
                continue
            acc = {k: [] for k in ("all", "pos", "neg", "lost_pos", "lost_neg")}
            for _ in range(REPEATS):
                pick = rng.choice(len(ids), size=n, replace=False)
                single = singletons(adjacency_full[np.ix_(pick, pick)])
                sub_y = y[pick]
                sub_had = had_neighbour[pick]
                acc["all"].append(single.mean())
                if (sub_y == 1).any():
                    acc["pos"].append(single[sub_y == 1].mean())
                if (sub_y == 0).any():
                    acc["neg"].append(single[sub_y == 0].mean())
                # Of those that had a neighbour in the full cohort, how many
                # lost it purely because the cohort got smaller?
                for key, want in (("lost_pos", 1), ("lost_neg", 0)):
                    mask = sub_had & (sub_y == want)
                    if mask.any():
                        acc[key].append(single[mask].mean())
            rows.append({
                "epitope": epitope, "n": n,
                **{f"singleton_{k}": float(np.mean(v)) if v else np.nan
                   for k, v in acc.items() if k in ("all", "pos", "neg")},
                "lost_neighbour_pos": float(np.mean(acc["lost_pos"]))
                if acc["lost_pos"] else np.nan,
                "lost_neighbour_neg": float(np.mean(acc["lost_neg"]))
                if acc["lost_neg"] else np.nan})
            print(f"  {epitope:<11} n={n:4d}  singleton all "
                  f"{rows[-1]['singleton_all']:.3f}  pos "
                  f"{rows[-1]['singleton_pos']:.3f}  neg "
                  f"{rows[-1]['singleton_neg']:.3f}   of those that had a "
                  f"neighbour, lost it: pos "
                  f"{rows[-1]['lost_neighbour_pos']:.3f} / neg "
                  f"{rows[-1]['lost_neighbour_neg']:.3f}")
    sweep = pd.DataFrame(rows)

    print("\n=== the decisive comparison: public thinned to private sizes ===")
    matched = []
    for private_epitope, n in sorted(PRIVATE_SIZES.items(), key=lambda kv: kv[1]):
        actual = float(observed.loc[observed["epitope"] == private_epitope,
                                    "singleton_all"].iloc[0])
        at_n = sweep[sweep["n"] == n]
        row = {"private epitope": private_epitope, "n": n,
               "actual singleton rate": actual,
               "public thinned, mean": float(at_n["singleton_all"].mean()),
               "public thinned, max": float(at_n["singleton_all"].max()),
               "public epitopes available": len(at_n)}
        matched.append(row)
        print(f"\n  {private_epitope} at n={n}: actually {actual:.1%} singleton")
        for _, r in at_n.sort_values("singleton_all").iterrows():
            print(f"     {r['epitope']:<11} thinned to {n}: "
                  f"{r['singleton_all']:.1%}")
        print(f"     -> public range {at_n['singleton_all'].min():.1%}"
              f"-{at_n['singleton_all'].max():.1%}, "
              f"private is {'inside' if at_n['singleton_all'].min() <= actual <= at_n['singleton_all'].max() else 'OUTSIDE'}")
    matched = pd.DataFrame(matched)

    print("\n=== positives against negatives ===")
    full = sweep[sweep["n"] == sweep.groupby("epitope")["n"].transform("max")]
    small = sweep[sweep["n"] == 130]
    merged = full.merge(small, on="epitope", suffixes=("_full", "_130"))
    merged["pos_increase"] = (merged["singleton_pos_130"]
                              - merged["singleton_pos_full"])
    merged["neg_increase"] = (merged["singleton_neg_130"]
                              - merged["singleton_neg_full"])
    print(merged[["epitope", "n_full", "singleton_pos_full", "singleton_pos_130",
                  "pos_increase", "singleton_neg_full", "singleton_neg_130",
                  "neg_increase"]].round(3).to_string(index=False))
    print(f"\n  mean increase going to n=130: positives "
          f"{merged['pos_increase'].mean():+.3f}, negatives "
          f"{merged['neg_increase'].mean():+.3f}")
    at130 = sweep[sweep["n"] == 130]
    print(f"\n  But conditional on HAVING had a neighbour, negatives lose it")
    print(f"  more readily: at n=130 the loss rate is "
          f"{at130['lost_neighbour_pos'].mean():.3f} for positives against "
          f"{at130['lost_neighbour_neg'].mean():.3f} for negatives, a ratio of "
          f"{at130['lost_neighbour_neg'].mean() / at130['lost_neighbour_pos'].mean():.2f}.")
    print("  A positive sits in a family with redundant membership; a negative's")
    print("  neighbour is more often a single chance edge that subsampling cuts.")

    print("\n=== how much does size matter against which epitope it is? ===")
    spread_at_130 = (at130["singleton_all"].max() - at130["singleton_all"].min())
    within = []
    for epitope in public:
        block = sweep[sweep["epitope"] == epitope]
        within.append(block["singleton_all"].max() - block["singleton_all"].min())
    print(f"  spread ACROSS epitopes at a fixed n=130:      "
          f"{spread_at_130:.1%}")
    print(f"  spread WITHIN an epitope across the whole sweep: "
          f"{np.mean(within):.1%} on average")
    print(f"  Sample size moves the singleton rate, but which epitope it is")
    print(f"  moves it {spread_at_130 / np.mean(within):.1f}x as far.")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    for epitope in public:
        block = sweep[sweep["epitope"] == epitope]
        ax.plot(block["n"], block["singleton_all"], marker="o", ms=3.5, lw=1.3,
                label=epitope)
    offsets = {"CINGVCWTV": (10, 10), "LLWNGPMAV": (-78, -6),
               "LTDEMIAQY": (8, -20)}
    for epitope, n in PRIVATE_SIZES.items():
        actual = float(observed.loc[observed["epitope"] == epitope,
                                    "singleton_all"].iloc[0])
        ax.scatter([n], [actual], marker="*", s=190, zorder=6,
                   edgecolors="black", linewidths=0.6, color="#d62728")
        ax.annotate(f"{epitope} {actual:.0%}", (n, actual),
                    textcoords="offset points",
                    xytext=offsets.get(epitope, (6, -16)), fontsize=7.5,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8,
                              pad=1))
    ax.set_xlabel("cohort size (subsampled from the full cohort)")
    ax.set_ylabel("fraction of TCRs with no neighbour at cosine 0.12")
    ax.set_ylim(0, 1.03)
    ax.legend(fontsize=7.5, frameon=False, loc="lower left", ncol=2,
              title="public, thinned", title_fontsize=7.5)
    ax.set_title("Singleton rate against cohort size — red stars are the "
                 "private cohorts at their own size", fontsize=10)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    png = base64.b64encode(buffer.getvalue()).decode()
    FIGURE.write_text(
        '<!doctype html><html lang="zh"><head><meta charset="utf-8">'
        "<title>singleton 率与队列大小</title><style>body{max-width:900px;"
        "margin:2.5rem auto;padding:0 1.2rem;font:15px/1.75 -apple-system,"
        '"Noto Sans CJK SC","Microsoft YaHei",sans-serif}img{width:100%}</style>'
        "</head><body><h1>singleton 率随队列大小怎么变</h1>"
        f'<img src="data:image/png;base64,{png}">'
        "<p>每条曲线是一个公开表位被随机稀疏化到不同规模后的 singleton 比例"
        f"（{REPEATS} 次重复）。红星是三个私有表位在它们自己的规模上的实际值。"
        "红星落在曲线族之内，说明该私有队列的稀疏程度可以用采样规模解释；"
        "落在之外则说明它比单纯的采样更稀疏。</p></body></html>",
        encoding="utf-8")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E55 - how much of \"singleton\" is sample size?\n\n"
        "## Observed, at each cohort's own size\n\n"
        + observed.round(4).to_markdown(index=False)
        + "\n\n## Public cohorts thinned to the private sizes\n\n"
        + matched.round(4).to_markdown(index=False)
        + "\n\n## Full sweep\n\n" + sweep.round(4).to_markdown(index=False)
        + "\n\n## Positives against negatives, full size against n=130\n\n"
        + merged[["epitope", "n_full", "singleton_pos_full",
                  "singleton_pos_130", "pos_increase", "singleton_neg_full",
                  "singleton_neg_130", "neg_increase"]].round(4).to_markdown(
                      index=False)
        + "\n\nCaveat: this subsamples a finite candidate library and cannot "
        "reconstruct a real antigen-specific repertoire, so it bounds the "
        "sampling explanation from below.\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)} and {FIGURE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
