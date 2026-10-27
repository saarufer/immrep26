#!/usr/bin/env python
"""How alike are the eight epitope cohorts, and which public one does each
private one resemble?

E85 needed this for one purpose -- finding each private epitope's nearest public
analogue -- and found something worth drawing properly: the private cohorts are
almost pure singletons (CINGVCWTV's largest family is two TCRs) while GILGFVFTL
and GLCTLVAML, the two epitopes nearly every intuition in this project was built
on, are the family-rich outliers. That makes TTDPSFLGRY the right public
analogue for the private split, and it is worth being able to see it rather than
read it off a table.

"Alike" is not one thing, so three axes are computed and drawn side by side,
because they order the epitopes differently and each answers a different
question:

  1. **Cohort structure** -- singleton share, family sizes, the spread of nn3.
     This is what decides whether the density mechanism has any dynamic range,
     so it is the axis that matters for transfer. Label-free, hence legitimate
     on the private rows.

  2. **TCR repertoire proximity** -- symmetrised Chamfer distance between the
     cohorts in SCEPTR space: for every TCR in A, the cosine distance to its
     nearest neighbour in B, averaged, and the same the other way. This asks
     whether the actual receptors look alike, which is a different question from
     whether they are organised alike.

  3. **Antigen** -- BLOSUM62 similarity between the peptides themselves,
     normalised by self-score. E75 established that six of the eight are
     HLA-A\\*02:01 nonamers, so this axis is mostly flat by construction and is
     included to show that it carries no grouping, not because it is expected to.

Purely descriptive. Nothing is fitted, and E85 already established that
structural matching *fails* its own leave-one-epitope-out test (MAE 0.0261
against the headroom model's 0.0157), so this map is a risk picture and must not
be used to re-weight the margin model.

Standardisation note: E85 standardised the structure features over the five
public epitopes only, because it was predicting public held-out folds. Here all
eight are standardised together, since the private cohorts are part of the
picture being drawn. Both versions of the public-only distances are printed so
the change is visible.

Run: .venv/bin/python scripts/77_epitope_similarity_map.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

REPORT = ROOT / "reports/52_epitope_similarity_map.md"
FIGURE = ROOT / "reports/figures/epitope_similarity_map.png"
FAMILY_RADIUS = 0.12
PRIVATE = {"CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"}
# E47, for the annotation only -- not an input to any distance here.
PUBLISHED = {"GILGFVFTL": 0.871, "YLQPRTFLL": 0.928, "NLVPMVATV": 0.809,
             "TTDPSFLGRY": 0.826, "GLCTLVAML": 0.690}


def cohort_structure(test, vectors, where):
    """The label-free family diagnostics, per epitope. Same as E85."""
    rows = {}
    for e, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vectors[[where[t] for t in ids]]
        dm = 1.0 - block @ block.T
        _, fam = connected_components(
            sp.csr_matrix(dm <= FAMILY_RADIUS), directed=False)
        sizes = pd.Series(fam).value_counts()
        iu = np.triu_indices(len(ids), 1)
        close = float((dm[iu] <= FAMILY_RADIUS).mean())
        np.fill_diagonal(dm, np.inf)
        nn3 = np.sort(dm, axis=1)[:, :3].mean(axis=1)
        q75, q25 = np.percentile(nn3, [75, 25])
        rows[e] = {
            "rows": len(ids),
            "singleton": float((sizes[fam] == 1).mean()),
            "max_family": int(sizes.max()),
            "close_pairs": close,
            "nn3_iqr": float(q75 - q25),
            "nn3_sd": float(nn3.std()),
            "hla": cohort["hla_short"].mode().iloc[0],
        }
    return pd.DataFrame(rows).T


def chamfer(test, vectors, where, names):
    """Symmetrised nearest-neighbour distance between two cohorts in SCEPTR space."""
    blocks = {}
    for e in names:
        ids = test[test["epitope_aa"] == e]["tcr_id"].to_numpy()
        blocks[e] = vectors[[where[t] for t in ids]]
    n = len(names)
    out = np.zeros((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i >= j:
                continue
            d = 1.0 - blocks[a] @ blocks[b].T
            out[i, j] = out[j, i] = 0.5 * (d.min(axis=1).mean()
                                           + d.min(axis=0).mean())
    return out


def peptide_similarity(names):
    """BLOSUM62 global alignment, normalised by self-score."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    n = len(names)
    out = np.zeros((n, n))
    self_score = {p: aligner.score(p, p) for p in names}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            out[i, j] = 2.0 * aligner.score(a, b) / (self_score[a] + self_score[b])
    return out


def heatmap(ax, mat, names, title, cmap, private, fmt="{:.2f}",
            mask_diagonal=False):
    """Diagonals are masked where they would otherwise dominate the scale, and
    the label colour is chosen from the cell's actual luminance rather than
    guessed from the colormap's direction."""
    m = np.array(mat, dtype=float)
    if mask_diagonal:
        np.fill_diagonal(m, np.nan)
    shown = np.ma.masked_invalid(m)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("#e9e9ec")
    im = ax.imshow(shown, cmap=cm)
    lo, hi = float(np.nanmin(m)), float(np.nanmax(m))
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    label = lambda e: f"{e}{'  ●' if e in private else ''}"  # noqa: E731
    ax.set_xticklabels([label(e) for e in names], rotation=90, fontsize=7.5)
    ax.set_yticklabels([label(e) for e in names], fontsize=7.5)
    for t, e in zip(ax.get_xticklabels(), names):
        t.set_color("#b8003c" if e in private else "#1a1a1a")
    for t, e in zip(ax.get_yticklabels(), names):
        t.set_color("#b8003c" if e in private else "#1a1a1a")
    for i in range(len(names)):
        for j in range(len(names)):
            if not np.isfinite(m[i, j]):
                continue
            r, g, b, _ = cm((m[i, j] - lo) / (hi - lo + 1e-12))
            lum = 0.299 * r + 0.587 * g + 0.114 * b
            ax.text(j, i, fmt.format(m[i, j]), ha="center", va="center",
                    fontsize=6.6, color="white" if lum < 0.55 else "#1a1a1a")
    for k, e in enumerate(names):
        if e not in private:
            continue
        ax.add_patch(Rectangle((k - .5, -.5), 1, len(names), fill=False,
                               ec="#b8003c", lw=1.1, alpha=.6))
        ax.add_patch(Rectangle((-.5, k - .5), len(names), 1, fill=False,
                               ec="#b8003c", lw=1.1, alpha=.6))
    ax.set_title(title, fontsize=9, pad=8)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=6.5)


def main() -> int:
    test = load_test()
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}

    struct = cohort_structure(test, vectors, where)
    names = list(struct.index)
    print("=== cohort structure, all eight ===")
    print(struct.to_string())

    # ---------- axis 1: cohort structure ----------
    feat_cols = ["singleton", "max_family", "close_pairs", "nn3_iqr", "nn3_sd"]
    F = struct[feat_cols].astype(float).to_numpy().copy()
    F = np.column_stack([
        F[:, 0], np.log(struct["rows"].astype(float).to_numpy()),
        np.log1p(F[:, 1]), F[:, 2], F[:, 3], F[:, 4]])
    Z = (F - F.mean(0)) / (F.std(0) + 1e-12)
    D_struct = np.linalg.norm(Z[:, None, :] - Z[None, :, :], axis=2)

    # the same thing standardised over the public five only, as E85 did
    pub_idx = [i for i, e in enumerate(names) if e not in PRIVATE]
    mu_p, sd_p = F[pub_idx].mean(0), F[pub_idx].std(0) + 1e-12
    Zp = (F - mu_p) / sd_p
    D_struct_pub = np.linalg.norm(Zp[:, None, :] - Zp[None, :, :], axis=2)
    print("\n  effect of standardising over all eight rather than the public five:")
    print("  (public-public distances only, E85's basis vs this one)")
    for i in pub_idx:
        for j in pub_idx:
            if i < j:
                print(f"    {names[i][:4]}-{names[j][:4]}  E85 {D_struct_pub[i,j]:5.2f}"
                      f"   all-eight {D_struct[i,j]:5.2f}")

    # ---------- axis 2: TCR repertoire proximity ----------
    print("\nmeasuring cross-cohort Chamfer distance in SCEPTR space ...")
    D_rep = chamfer(test, vectors, where, names)

    # ---------- axis 3: the antigen ----------
    S_pep = peptide_similarity(names)

    # ---------- the ranking the question asks for ----------
    lines = ["# E86 — A similarity map of the eight epitope cohorts", "",
             "Three axes, because they order the epitopes differently.",
             "Descriptive only: E85 established that matching on structure *fails*",
             "its own leave-one-epitope-out test, so this is a risk picture and not",
             "an input to the margin model.",
             "Generated by `scripts/77_epitope_similarity_map.py`.", ""]

    lines += ["## The cohorts", "",
              "| epitope | split | HLA | rows | singleton % | largest family "
              "| close pairs % | nn3 IQR | nn3 sd |", "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for e in names:
        r = struct.loc[e]
        lines.append(f"| {e} | {'private' if e in PRIVATE else 'public'} | {r['hla']} "
                     f"| {int(r['rows'])} | {float(r['singleton'])*100:.1f} "
                     f"| {int(r['max_family'])} | {float(r['close_pairs'])*100:.2f} "
                     f"| {float(r['nn3_iqr']):.4f} | {float(r['nn3_sd']):.4f} |")
    lines.append("")

    rank_tables = {}
    for axis, mat, higher_is_closer in (
            ("cohort structure", D_struct, False),
            ("TCR repertoire", D_rep, False),
            ("peptide (BLOSUM62)", S_pep, True)):
        print(f"\n=== nearest public epitope to each private one: {axis} ===")
        tbl = []
        for e in sorted(PRIVATE):
            i = names.index(e)
            others = [(names[j], mat[i, j]) for j in pub_idx]
            others.sort(key=lambda kv: -kv[1] if higher_is_closer else kv[1])
            print(f"  {e:<11} " + "  ".join(
                f"{n[:4]} {v:.2f}" for n, v in others))
            tbl.append([e] + [f"{n} ({v:.2f})" for n, v in others])
        rank_tables[axis] = tbl

    lines += ["## Nearest public epitope to each private one", ""]
    for axis, tbl in rank_tables.items():
        lines += [f"### {axis}", "",
                  "| private | 1st | 2nd | 3rd | 4th | 5th |",
                  "|---|---|---|---|---|---|"]
        for r in tbl:
            lines.append("| " + " | ".join(r) + " |")
        lines.append("")

    # how much do the three axes agree?
    from scipy.stats import spearmanr
    iu = np.triu_indices(len(names), 1)
    print("\n=== do the three axes agree? Spearman over the 28 epitope pairs ===")
    pairs = [("structure", D_struct), ("repertoire", D_rep),
             ("peptide", -S_pep)]   # negate so all three are distances
    agree = []
    for a in range(3):
        for b in range(a + 1, 3):
            rho = spearmanr(pairs[a][1][iu], pairs[b][1][iu]).statistic
            print(f"  {pairs[a][0]:<11} vs {pairs[b][0]:<11} ρ = {rho:+.3f}")
            agree.append([pairs[a][0], pairs[b][0], rho])
    lines += ["## Do the three axes agree?", "",
              "Spearman over the 28 epitope pairs, all three as distances.", "",
              "| axis | axis | ρ |", "|---|---|---:|"]
    for r in agree:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:+.3f} |")
    lines.append("")

    # The bar panel puts tcrdist3's published AUROC under each public epitope,
    # which invites reading "closer to the private cohorts" as "easier". Check
    # it before letting the figure imply it.
    print("\n=== does proximity to the private cohorts track difficulty? ===")
    priv_idx = [i for i, e in enumerate(names) if e in PRIVATE]
    mean_to_priv = {names[j]: float(np.mean([D_struct[i, j] for i in priv_idx]))
                    for j in pub_idx}
    xs = [mean_to_priv[names[j]] for j in pub_idx]
    ys = [PUBLISHED[names[j]] for j in pub_idx]
    rho_diff = spearmanr(xs, ys).statistic
    for j in pub_idx:
        print(f"  {names[j]:<11} mean distance to the three private "
              f"{mean_to_priv[names[j]]:5.2f}   tcrdist3 {PUBLISHED[names[j]]:.3f}")
    print(f"  Spearman = {rho_diff:+.3f} over five points -- "
          f"{'no' if abs(rho_diff) < 0.9 else 'a'} usable relationship")
    lines += ["## Does proximity to the private cohorts track difficulty?", "",
              "The bar panel labels each public epitope with tcrdist3's published",
              "AUROC, which invites reading *closer* as *easier*. It does not hold.",
              "", "| public epitope | mean distance to the three private | tcrdist3 |",
              "|---|---:|---:|"]
    for j in pub_idx:
        lines.append(f"| {names[j]} | {mean_to_priv[names[j]]:.2f} "
                     f"| {PUBLISHED[names[j]]:.3f} |")
    lines += ["", f"Spearman **{rho_diff:+.3f}** over five points.", ""]

    # ---------- the figure ----------
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titlesize": 9})
    fig = plt.figure(figsize=(15.5, 10.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0],
                          hspace=0.58, wspace=0.34,
                          left=0.065, right=0.975, top=0.895, bottom=0.135)

    heatmap(fig.add_subplot(gs[0, 0]), D_struct, names,
            "1  Cohort structure — standardised distance\n"
            "(singleton share, family sizes, nn3 spread; label-free)",
            "viridis_r", PRIVATE, mask_diagonal=True)
    heatmap(fig.add_subplot(gs[0, 1]), D_rep, names,
            "2  TCR repertoire — symmetrised Chamfer distance\n"
            "(nearest-neighbour cosine between cohorts, SCEPTR)",
            "magma_r", PRIVATE, fmt="{:.3f}", mask_diagonal=True)
    heatmap(fig.add_subplot(gs[0, 2]), S_pep, names,
            "3  Antigen — BLOSUM62 peptide similarity\n"
            "(six of eight are HLA-A*02:01 nonamers)",
            "Blues", PRIVATE, mask_diagonal=True)

    # dendrograms on axes 1 and 2
    for col, (mat, label) in enumerate((
            (D_struct, "cohort structure"), (D_rep, "TCR repertoire"))):
        ax = fig.add_subplot(gs[1, col])
        link = linkage(squareform(mat, checks=False), method="average")
        dn = dendrogram(link, labels=names, ax=ax, color_threshold=0,
                        above_threshold_color="#555")
        ax.set_title(f"Hierarchical clustering — {label}", fontsize=9, pad=6)
        ax.tick_params(axis="x", labelsize=7.5, rotation=90)
        ax.tick_params(axis="y", labelsize=6.5)
        for t in ax.get_xticklabels():
            if t.get_text() in PRIVATE:
                t.set_color("#b8003c")
                t.set_fontweight("bold")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    # the ranking panel: each private epitope's distance to the five public ones
    ax = fig.add_subplot(gs[1, 2])
    pub_names = [names[j] for j in pub_idx]
    order = np.argsort([D_struct[names.index("CINGVCWTV"), j] for j in pub_idx])
    pub_sorted = [pub_names[k] for k in order]
    width, colours = 0.26, ("#b8003c", "#e07a3f", "#7a4fa3")
    for k, e in enumerate(sorted(PRIVATE)):
        i = names.index(e)
        vals = [D_struct[i, names.index(p)] for p in pub_sorted]
        ax.bar(np.arange(len(pub_sorted)) + (k - 1) * width, vals, width,
               label=e, color=colours[k], edgecolor="white", linewidth=.5)
    ax.set_xticks(range(len(pub_sorted)))
    ax.set_xticklabels([f"{p}\n(tcrdist3 {PUBLISHED[p]:.2f})" for p in pub_sorted],
                       fontsize=7, rotation=0)
    ax.set_ylabel("standardised structural distance", fontsize=8)
    ax.set_title("Each private cohort's distance to the five public ones\n"
                 "(lower = more alike; cohort-structure axis)", fontsize=9, pad=6)
    ax.legend(fontsize=7.5, frameon=False)
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(axis="y", alpha=.25, lw=.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle("IMMREP26 — how alike are the eight epitope cohorts?    "
                 "● = private (602 rows, decides the ranking)",
                 fontsize=12.5, y=0.965)
    fig.text(0.065, 0.016,
             "E86. Descriptive. E85 showed that matching the margin model on "
             "structural similarity fails its own leave-one-epitope-out test "
             "(MAE 0.0261 vs 0.0157), so this is a risk picture, not a re-weighting.",
             fontsize=7.5, color="#555")
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=165)
    print(f"\nwrote {FIGURE.relative_to(ROOT)}")

    lines += ["## Figure", "", f"![similarity map](figures/{FIGURE.name})", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
