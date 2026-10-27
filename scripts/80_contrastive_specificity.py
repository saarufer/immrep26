#!/usr/bin/env python
"""E88 found peptide-specific signal. Two things follow, and both are tested here.

E88 established that the SCEPTR distance space places a validated TCR nearer
its own epitope's cohort than the other seven at **2.6x the permutation null**,
and that this survives removing the query's family, donor, publication and
forty nearest neighbours from every cohort, and survives dropping the V and J
genes from the representation. The mechanism does carry antigen information.

That raises exactly two questions, and they pull in opposite directions.

**1. Is it antigen selection, or is it sequence-space geometry?** E81's sharpest
control was a SCEPTR trained purely on OLGA-generated sequences -- zero antigen
exposure anywhere in its training -- which scored 0.8337 on task A against the
real model's 0.8414. If that same antigen-naive model also separates the
cohorts at 2.6x, then E88's result is a statement about sequence space, not
about recognition, and the specificity reading collapses. Run the whole E88
pipeline on `synthetic_data` and `shuffled_data` as well as the real model.

**2. Can it be turned into a feature?** Task A asks "is this annotation real",
and D answers it with typicality: how dense is this TCR's neighbourhood inside
its own cohort. That quantity is antigen-agnostic, which is why an OLGA-trained
model nearly matches it. But E88's statistic is different in kind -- it is
*contrastive*:

    contrast(t) = z_own(t) - max_{other} z_other(t)

"how much more typical is this TCR of its annotated epitope than of any other
epitope". A spurious annotation should look typical of TCRs in general and
*not* preferentially typical of its own cohort, so this should separate
validated from non-validated for a reason that is not typicality. The standing
`typicality_gate` (0.6378) exists precisely to catch a feature that only looks
new, and the pre-registered inclusion bar is +0.02 under nested LOEO.

The contrastive feature is label-free and computable on the private rows, since
it needs only cohort membership.

Run: .venv/bin/python scripts/80_contrastive_specificity.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/55_contrastive_specificity.md"
FAMILY_RADIUS = 0.12
SIGMA = 0.1
N_SUBSAMPLE = 110
REPEATS = 30
SEED = 0
TYPICALITY_GATE = 0.6378
ABLATIONS = ("synthetic_data", "shuffled_data")


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def switch_accuracy(Dm, ids, ep, epitopes, members, query, fam, panel, rng,
                    drop_nn=20):
    """E88's statistic, level (ii)+excision, on an arbitrary distance matrix."""
    q = query[np.isin(ep[query], panel)]
    true_ep = ep[q]
    hits, nulls = [], []
    for _ in range(REPEATS):
        pool = {e: rng.choice(members[e], size=min(N_SUBSAMPLE, len(members[e])),
                              replace=False) for e in panel}
        Zs = np.zeros((len(q), len(panel)))
        for j, e in enumerate(panel):
            ref = pool[e]
            dqr = Dm[np.ix_(q, ref)]
            keep = (q[:, None] != ref[None, :]) & (fam[q][:, None] != fam[ref][None, :])
            if drop_nn:
                np.put_along_axis(keep, np.argsort(dqr, axis=1)[:, :drop_nn],
                                  False, axis=1)
            k = np.exp(-dqr / SIGMA) * keep
            n = keep.sum(axis=1)
            Zs[:, j] = np.where(n > 0, k.sum(axis=1) / np.maximum(n, 1), np.nan)
        Zs = (Zs - np.nanmean(Zs, 0)) / (np.nanstd(Zs, 0) + 1e-12)
        picked = np.array(panel)[np.argmax(Zs, axis=1)]
        hits.append(float((picked == true_ep).mean()))
        nulls.append(float((picked == rng.permutation(true_ep)).mean()))
    return float(np.mean(hits)), float(np.mean(nulls))


def contrast_feature(Dm, ids, ep, epitopes, members, fam, drop_nn=0):
    """z_own - max_other, for EVERY TCR (label-free, so private-computable)."""
    n_all = len(ids)
    Z = np.zeros((n_all, len(epitopes)))
    for j, e in enumerate(epitopes):
        ref = members[e]
        dqr = Dm[:, ref]
        keep = (np.arange(n_all)[:, None] != ref[None, :]) \
            & (fam[:, None] != fam[ref][None, :])
        if drop_nn:
            np.put_along_axis(keep, np.argsort(dqr, axis=1)[:, :drop_nn],
                              False, axis=1)
        k = np.exp(-dqr / SIGMA) * keep
        cnt = keep.sum(axis=1)
        Z[:, j] = np.where(cnt > 0, k.sum(axis=1) / np.maximum(cnt, 1), np.nan)
    Z = (Z - np.nanmean(Z, 0)) / (np.nanstd(Z, 0) + 1e-12)
    own = np.array([list(epitopes).index(e) for e in ep])
    z_own = Z[np.arange(n_all), own]
    other = Z.copy()
    other[np.arange(n_all), own] = -np.inf
    return pd.Series(z_own - other.max(axis=1), index=ids), \
        pd.Series(z_own - np.nanmean(np.where(np.isfinite(other), other, np.nan), 1),
                  index=ids)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    ids = test["tcr_id"].to_numpy()
    ep = test["epitope_aa"].to_numpy()
    epitopes = sorted(set(ep))
    members = {e: np.where(ep == e)[0] for e in epitopes}
    lab = labels.set_index("tcr_id")
    y = pd.Series(ids).map(lab["pos"]).to_numpy()
    query = np.where(y == 1)[0]
    a2 = sorted({e for e, h in zip(ep, test["hla_short"]) if h == "A2"})

    def dist_from(vectors):
        W = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
        D = 1.0 - W @ W.T
        np.fill_diagonal(D, 0.0)
        return D

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    pos_of = {t: i for i, t in enumerate(cached["tcr_id"])}
    D_real = dist_from(cached["vectors"][[pos_of[t] for t in ids]])
    _, fam = connected_components(sp.csr_matrix(D_real <= FAMILY_RADIUS),
                                  directed=False)

    lines = ["# E89 — Is E88's specificity antigen selection, and can it be a feature?",
             "", "Generated by `scripts/80_contrastive_specificity.py`.", ""]

    # ---------- 1. the antigen-naive control ----------
    print("=== 1. does an antigen-naive model separate the cohorts too? ===")
    print("  E81's control: a SCEPTR trained purely on OLGA sequences has zero")
    print("  antigen exposure. If it matches, E88 is about sequence space.")
    import sceptr
    from sceptr import variant as V
    from libtcrlm.schema.tcr import TravGene, TrbvGene
    sceptr.disable_hardware_acceleration()
    slash = lambda g: g.replace("_", "/") if isinstance(g, str) else None  # noqa: E731
    frame = pd.DataFrame({
        "TRAV": [slash(x) for x in test["TRAV_IMGT"]],
        "CDR3A": test["cdr3_alpha_aa"].tolist(),
        "TRAJ": [slash(x) for x in test["TRAJ_IMGT"]],
        "TRBV": [slash(x) for x in test["TRBV_IMGT"]],
        "CDR3B": test["cdr3_beta_aa"].tolist(),
        "TRBJ": [slash(x) for x in test["TRBJ_IMGT"]]})
    for column, vocab in (("TRAV", {g.name for g in TravGene}),
                          ("TRBV", {g.name for g in TrbvGene})):
        base = frame[column].fillna("").str.split("*").str[0]
        frame.loc[~base.isin(vocab), column] = None

    mats = {"real SCEPTR (default)": D_real}
    for name in ABLATIONS:
        print(f"  embedding {name} ...")
        mats[name] = dist_from(getattr(V, name)().calc_vector_representations(frame))

    print(f"\n  {'representation':<26}{'switch acc':>12}{'null':>8}{'ratio':>8}"
          f"{'task A (D-style)':>18}")
    ctrl_rows = []
    for name, Dm in mats.items():
        rng = np.random.default_rng(SEED)
        acc, null = switch_accuracy(Dm, ids, ep, epitopes, members, query, fam,
                                    a2, rng)
        # task A on the same matrix, for the comparison E81 drew
        topk = {}
        for e in epitopes:
            mem = members[e]
            d = Dm[np.ix_(mem, mem)].copy()
            np.fill_diagonal(d, np.inf)
            for i, t in enumerate(ids[mem]):
                topk[t] = -np.sort(d[i])[:8].mean()
        a_auc = hyp.mean_auc(labels, pd.Series(topk))
        print(f"  {name:<26}{acc:12.3f}{null:8.3f}{acc/null:8.2f}{a_auc:18.4f}")
        ctrl_rows.append([name, acc, null, acc / null, a_auc])

    lines += ["## 1. Antigen selection, or sequence-space geometry?", "",
              "E81's control: a SCEPTR trained purely on OLGA-generated sequences has",
              "zero antigen exposure anywhere in its training. HLA-A2 panel, family",
              "excluded, 20 nearest neighbours excised from every cohort.", "",
              "| representation | switch accuracy | null | ratio | task A AUC |",
              "|---|---:|---:|---:|---:|"]
    for r in ctrl_rows:
        lines.append(f"| {r[0]} | **{r[1]:.3f}** | {r[2]:.3f} | **{r[3]:.2f}** "
                     f"| {r[4]:.4f} |")
    lines.append("")

    # ---------- 2. the contrastive feature on task A ----------
    print("\n=== 2. the contrastive feature on task A ===")
    print("  contrast(t) = z_own(t) - max_other z(t), label-free, private-computable")
    feats = {}
    for drop in (0, 20):
        cmax, cmean = contrast_feature(D_real, ids, ep, epitopes, members, fam,
                                       drop_nn=drop)
        feats[f"z_own - max_other (m={drop})"] = cmax
        feats[f"z_own - mean_other (m={drop})"] = cmean
    # the reference points
    topk = {}
    for e in epitopes:
        mem = members[e]
        d = D_real[np.ix_(mem, mem)].copy()
        np.fill_diagonal(d, np.inf)
        for i, t in enumerate(ids[mem]):
            topk[t] = -np.sort(d[i])[:8].mean()
    typicality = pd.Series(topk)
    print(f"  {'feature':<34}{'task A AUC':>12}{'ρ w/ typicality':>18}"
          f"{'gate':>8}")
    feat_rows = []
    for name, series in feats.items():
        auc = hyp.mean_auc(labels, series)
        rho = float(np.mean([
            spearmanr(series.reindex(ids[members[e]]),
                      typicality.reindex(ids[members[e]])).statistic
            for e in epitopes]))
        passes = "pass" if abs(rho) < 0.6 else "FAIL"
        print(f"  {name:<34}{auc:12.4f}{rho:+18.3f}{passes:>8}")
        feat_rows.append([name, auc, rho, passes])
    print(f"  {'typicality itself (top-8 density)':<34}"
          f"{hyp.mean_auc(labels, typicality):12.4f}{1.0:+18.3f}")

    lines += ["## 2. The contrastive feature on task A", "",
              "`contrast(t) = z_own(t) - max_other z(t)` -- how much more typical a",
              "TCR is of its annotated epitope than of any other. Label-free, so it",
              "is computable on the private rows.", "",
              "| feature | task A AUC | ρ with typicality | typicality gate |",
              "|---|---:|---:|---|"]
    for r in feat_rows:
        lines.append(f"| {r[0]} | **{r[1]:.4f}** | {r[2]:+.3f} | {r[3]} |")
    lines.append(f"| typicality itself (top-8 density) "
                 f"| {hyp.mean_auc(labels, typicality):.4f} | +1.000 | — |")
    lines.append("")

    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
