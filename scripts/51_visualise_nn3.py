#!/usr/bin/env python
"""Visualise the nn3_dist label distribution, for intuition rather than a number.

E53 established two things numerically that are hard to feel from a table: the
validation rate declines smoothly in nn3_dist with no knee, and the private
cohorts sit far to the right of the public ones. Four panels, each answering one
question a reader would ask.

Axis labels are English on purpose: a CJK font is not guaranteed in this
environment and a missing glyph renders as a box, which would be worse than
English in a figure whose surrounding prose is Chinese.

Output: reports/nn3_distribution.html, self-contained (figures base64-inlined).

Run: .venv/bin/python scripts/51_visualise_nn3.py
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

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

OUT = ROOT / "reports/nn3_distribution.html"
RADIUS = 0.12
POS_COLOUR = "#2a6fb3"
NEG_COLOUR = "#c4443a"
PRIVATE = ("CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY")

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "figure.dpi": 150,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
})


def neighbours(test):
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    rows = []
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        d = 1.0 - block @ block.T
        np.fill_diagonal(d, np.inf)
        ordered = np.sort(d, axis=1)
        _, fam = connected_components(
            sp.csr_matrix(d <= RADIUS), directed=False)
        sizes = pd.Series(fam).value_counts()
        clustered = np.isin(fam, sizes[sizes >= 2].index)
        for i, t in enumerate(ids):
            rows.append({"tcr_id": t, "epitope_aa": epitope,
                         "nn1": ordered[i, 0], "nn3": ordered[i, :3].mean(),
                         "clustered": bool(clustered[i])})
    return pd.DataFrame(rows).set_index("tcr_id")


def embed(fig) -> str:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode()


def figure_strips(frame, labels, order):
    """Every labelled TCR as a dot: positives above the line, negatives below.
    The most direct view of how much the two groups overlap."""
    fig, axes = plt.subplots(len(order), 1, figsize=(8.2, 1.35 * len(order)),
                             sharex=True)
    rng = np.random.default_rng(0)
    for ax, epitope in zip(axes, order):
        block = labels[labels["epitope_aa"] == epitope]
        x = block["tcr_id"].map(frame["nn3"]).to_numpy()
        y = block["pos"].to_numpy()
        auc = roc_auc_score(y, -x)
        for mask, sign, colour, name in ((y == 1, 1, POS_COLOUR, "validated"),
                                         (y == 0, -1, NEG_COLOUR, "not")):
            ax.scatter(x[mask], sign * (0.25 + rng.random(mask.sum()) * 0.55),
                       s=5, c=colour, alpha=0.5, linewidths=0, label=name)
        ax.axvline(RADIUS, color="#333", lw=1, ls="--")
        ax.axhline(0, color="#999", lw=0.6)
        ax.set_yticks([])
        ax.set_ylim(-1.0, 1.0)
        ax.set_ylabel(epitope, rotation=0, ha="right", va="center", fontsize=8)
        ax.text(0.998, 0.97, f"nn3 alone AUC {auc:.3f}   n={len(block)}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                color="#444", bbox=dict(facecolor="white", edgecolor="none",
                                        alpha=0.85, pad=1.5))
    axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.08), fontsize=8,
                   frameon=False, markerscale=2.5, ncol=2)
    axes[0].text(RADIUS + 0.008, 0.80, "radius 0.12", fontsize=7, color="#333",
                 bbox=dict(facecolor="white", edgecolor="none", alpha=0.85,
                           pad=1))
    axes[-1].set_xlabel("nn3_dist  (mean cosine distance to the 3 nearest "
                        "cohort members; left = has close relatives)")
    fig.suptitle("Every labelled TCR, positioned by how close its nearest "
                 "relatives are", y=1.06, fontsize=10.5)
    return embed(fig)


def figure_densities(frame, labels, order):
    fig, axes = plt.subplots(1, len(order), figsize=(2.55 * len(order), 2.5),
                             sharey=True)
    bins = np.linspace(0, 0.8, 33)
    for ax, epitope in zip(axes, order):
        block = labels[labels["epitope_aa"] == epitope]
        x = block["tcr_id"].map(frame["nn3"]).to_numpy()
        y = block["pos"].to_numpy()
        for mask, colour, name in ((y == 1, POS_COLOUR, "validated"),
                                   (y == 0, NEG_COLOUR, "not validated")):
            ax.hist(x[mask], bins=bins, density=True, alpha=0.55,
                    color=colour, label=name)
        ax.axvline(RADIUS, color="#333", lw=1, ls="--")
        ax.set_title(f"{epitope}\n{int(y.sum())}+ / {int((1 - y).sum())}-",
                     fontsize=8.5)
        ax.set_xlabel("nn3_dist")
    axes[0].set_ylabel("density")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("The same thing as overlapping histograms, per epitope",
                 y=1.04, fontsize=10.5)
    return embed(fig)


def figure_rate(frame, labels, order):
    fig, ax = plt.subplots(figsize=(7.4, 3.5))
    edges = np.linspace(0, 0.75, 16)
    centres = (edges[:-1] + edges[1:]) / 2
    for epitope in order:
        block = labels[labels["epitope_aa"] == epitope].copy()
        block["nn3"] = block["tcr_id"].map(frame["nn3"])
        rate, keep = [], []
        for lo, hi, centre in zip(edges[:-1], edges[1:], centres):
            band = block[(block["nn3"] >= lo) & (block["nn3"] < hi)]
            if len(band) >= 8:
                rate.append(band["pos"].mean())
                keep.append(centre)
        ax.plot(keep, rate, marker="o", ms=3.5, lw=1.1, alpha=0.75,
                label=epitope)
    everything = labels.copy()
    everything["nn3"] = everything["tcr_id"].map(frame["nn3"])
    pooled, keep = [], []
    for lo, hi, centre in zip(edges[:-1], edges[1:], centres):
        band = everything[(everything["nn3"] >= lo) & (everything["nn3"] < hi)]
        if len(band) >= 20:
            pooled.append(band["pos"].mean())
            keep.append(centre)
    ax.plot(keep, pooled, color="#111", lw=2.4, label="all five pooled",
            zorder=5)
    ax.axvline(RADIUS, color="#333", lw=1, ls="--")
    ax.text(RADIUS + 0.006, 0.03, "radius 0.12", fontsize=7.5, color="#333")
    ax.axhline(0.5, color="#999", lw=0.6, ls=":")
    ax.set_xlabel("nn3_dist")
    ax.set_ylabel("fraction validated")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=7.5, frameon=False, ncol=2)
    ax.set_title("Validation rate against nn3_dist — the slope is smooth, "
                 "and 0.12 is not where anything happens", fontsize=10)
    return embed(fig)


def figure_transfer(frame, test):
    fig, ax = plt.subplots(figsize=(7.4, 3.5))
    for epitope, cohort in test.groupby("epitope_aa"):
        values = np.sort(cohort["tcr_id"].map(frame["nn3"]).to_numpy())
        y = np.arange(1, len(values) + 1) / len(values)
        private = epitope in PRIVATE
        ax.plot(values, y, lw=2.2 if private else 1.2,
                ls="--" if private else "-", alpha=0.95 if private else 0.7,
                label=f"{epitope}{' (private)' if private else ''}")
    ax.axvline(RADIUS, color="#333", lw=1, ls="--")
    ax.text(RADIUS + 0.006, 0.95, "radius 0.12", fontsize=7.5, color="#333")
    ax.set_xlabel("nn3_dist")
    ax.set_ylabel("cumulative fraction of the cohort")
    ax.set_title("Sparseness is not a private-only problem: TTDPSFLGRY and "
                 "NLVPMVATV are\npublic and sparser than two of the three "
                 "private cohorts", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    return embed(fig)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    frame = neighbours(test)
    order = sorted(labels["epitope_aa"].unique())

    stats = []
    for epitope, cohort in test.groupby("epitope_aa"):
        values = cohort["tcr_id"].map(frame["nn3"])
        block = labels[labels["epitope_aa"] == epitope]
        auc = (roc_auc_score(block["pos"],
                             -block["tcr_id"].map(frame["nn3"]).to_numpy())
               if len(block) else np.nan)
        stats.append({
            "epitope": epitope, "split": cohort["split"].iloc[0],
            "n": len(cohort), "median nn3": values.median(),
            "below 0.12": f"{(values <= RADIUS).mean():.1%}",
            "nn3 alone AUC": f"{auc:.3f}" if not np.isnan(auc) else "—",
        })
    table = pd.DataFrame(stats).sort_values("median nn3")

    panels = [
        ("每条有标签的 TCR 一个点", figure_strips(frame, labels, order),
         "正样本在上、负样本在下，横轴是它到最近三个同队列成员的平均余弦距离。"
         "越靠左代表亲戚越近。两群确实分开，但<b>是连续地分开，不是被某条线分开</b>——"
         "虚线是项目里用的 0.12 半径，它落在左边一小撮里，右边绝大多数点都在它之外。"),
        ("同一件事画成直方图", figure_densities(frame, labels, order),
         "红蓝重叠的面积就是 nn3 单独排序会犯的错。注意 GLCTLVAML 的两条几乎完全重合"
         "（它的 nn3 单独 AUC 最低），而 YLQPRTFLL 分得最开。"),
        ("验证率随 nn3 怎么变", figure_rate(frame, labels, order),
         "黑粗线是五个表位合并。它从 0.9 一路平滑滑到 0.1，<b>没有任何一处出现台阶</b>。"
         "E53 量过：最大单步 −0.216，平均步长 0.107，比值只有 2.0。"
         "真有阈值的话会有一步远远甩开其他。所以「是不是 singleton」是一个"
         "<b>有损的二值摘要</b>，不是一个真实的类别边界。"),
        ("稀疏程度：私有 vs 公开", figure_transfer(frame, test),
         "这张图<b>完全不需要标签</b>，所以三个私有表位也画得出来（粗虚线）。"
         "它纠正了一个容易得出的悲观结论：按 nn3 中位数排序是 "
         "GIL 0.190 &lt; GLC 0.214 &lt; <b>LLW 0.414</b> &lt; <b>CIN 0.447</b> &lt; "
         "YLQ 0.452 &lt; TTD 0.532 &lt; NLV 0.535 &lt; <b>LTD 0.595</b>——"
         "<b>TTDPSFLGRY 和 NLVPMVATV 比两个私有表位还稀疏</b>，而它们的 nn3 单独 AUC "
         "是 0.860 和 0.845。只有 LTDEMIAQY 比所有公开表位都稀疏。"
         "所以私有队列的稀疏本身<b>不是</b>低分的充分理由——这和 E21 的结论一致："
         "稀疏度不能预测 AUC。"),
    ]

    html = f"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8">
<title>nn3_dist 的标签分布</title>
<style>
 body{{max-width:980px;margin:2.5rem auto;padding:0 1.2rem;
  font:15px/1.75 -apple-system,"Segoe UI","Noto Sans CJK SC","Microsoft YaHei",sans-serif;
  color:#1a1a1a}}
 h1{{font-size:1.5rem;margin:0 0 .3rem}}
 h2{{font-size:1.1rem;margin:2.6rem 0 .5rem;padding-bottom:.3rem;
  border-bottom:1px solid #e3e3e3}}
 .sub{{color:#666;margin:0 0 1.8rem;font-size:.92rem}}
 img{{width:100%;height:auto;margin:.5rem 0}}
 p{{margin:.5rem 0 0}}
 table{{border-collapse:collapse;font-size:.88rem;margin:1rem 0}}
 th,td{{padding:.3rem .7rem;border-bottom:1px solid #eaeaea;text-align:right}}
 th:first-child,td:first-child{{text-align:left}}
 tr.private td{{background:#fbf6ee}}
 code{{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px;font-size:.9em}}
</style></head><body>
<h1>nn3_dist 的标签分布</h1>
<p class="sub">每条 TCR 到它<b>同表位队列内</b>最近三个成员的平均余弦距离，
对照验证标签。由 <code>scripts/51_visualise_nn3.py</code> 生成；
图内标签用英文以避免中文字体缺失。</p>
<table><tr><th>表位</th><th>划分</th><th>n</th><th>nn3 中位数</th>
<th>≤0.12 占比</th><th>nn3 单独 AUC</th></tr>
{"".join(
    f'<tr class="{"private" if r["split"] == "private" else ""}">'
    f'<td>{r["epitope"]}</td><td>{r["split"]}</td><td>{r["n"]}</td>'
    f'<td>{r["median nn3"]:.3f}</td><td>{r["below 0.12"]}</td>'
    f'<td>{r["nn3 alone AUC"]}</td></tr>'
    for _, r in table.iterrows())}
</table>
{"".join(f'<h2>{title}</h2><img src="data:image/png;base64,{png}" alt="{title}">'
         f'<p>{note}</p>' for title, png, note in panels)}
<h2>一句不能省的说明</h2>
<p>0.12 这个半径在项目里<b>只用于分析记账</b>——定义克隆型家族和配对类别——
<b>不是预测器</b>。所以第三张图不代表模型白扔了 0.16 AUC。它代表的是：
E22、E50、E52 这些基于类别的分析都建在一个有损摘要之上，这是局限，不是可以去捡的分数。
而连续变量本身已经被利用了：提交的候选 D 用 top-8 聚合直接读它，这就是它在双孤立配对上
拿到 +0.0143 的原因。</p>
</body></html>"""
    OUT.write_text(html, encoding="utf-8")
    print(table.to_string(index=False))
    print(f"\nwrote {OUT.relative_to(ROOT)} "
          f"({OUT.stat().st_size / 1024:.0f} KB, self-contained)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
