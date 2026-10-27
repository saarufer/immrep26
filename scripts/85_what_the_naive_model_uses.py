#!/usr/bin/env python
"""E89 left one question open: what is the antigen-naive model actually using?

E88 found that the SCEPTR distance space places a validated TCR nearer its own
epitope's cohort than the other seven at **2.6x** the permutation null,
surviving family, donor, publication, the 40 nearest neighbours in every
cohort, and the removal of V and J genes. E89 then found that a SCEPTR trained
purely on OLGA-generated sequences — zero antigen exposure — does it at **2.67x**,
slightly better than the real model. The separation is therefore a property of
sequence space rather than of learned recognition.

That is a negative result about the representation and a positive one about the
data, and it leaves the obvious question unanswered: **sequence space along
which axis?** An OLGA-trained model has learned the shape of V(D)J
recombination. If that is enough to tell GILGFVFTL's repertoire from
NLVPMVATV's, the difference must live in quantities a generative model can
see — CDR3 length, amino-acid composition, generation probability, gene usage.

So the SCEPTR distance is replaced with explicit, interpretable feature spaces
and E88's whole pipeline is re-run on each: same HLA-A2 panel, same family
exclusion, same 20-nearest-neighbour excision from every cohort, same
permutation null. Whatever recovers the 2.6x is what the naive model is using.

A second, blunter view is reported alongside: the per-epitope means themselves,
since if CDR3β is two residues longer on one epitope than another then the
whole effect may be legible without any model at all.

Run: .venv/bin/python scripts/85_what_the_naive_model_uses.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/62_what_the_naive_model_uses.md"
FAMILY_RADIUS = 0.12
SIGMA = 0.1
N_SUBSAMPLE = 110
REPEATS = 30
DROP_NN = 20
SEED = 0
AA = "ACDEFGHIKLMNPQRSTVWY"


def zdist(F: np.ndarray) -> np.ndarray:
    """Euclidean distance on column-standardised features, rescaled so the
    median off-diagonal distance matches SCEPTR cosine's, which keeps a single
    kernel bandwidth meaningful across very different feature spaces."""
    Z = (F - F.mean(0)) / (F.std(0) + 1e-12)
    D = np.sqrt(np.maximum(
        ((Z ** 2).sum(1)[:, None] + (Z ** 2).sum(1)[None, :] - 2 * Z @ Z.T), 0))
    iu = np.triu_indices(len(D), 1)
    med = np.median(D[iu])
    D = D * (0.45 / (med + 1e-12))      # 0.45 ~ median SCEPTR cosine distance
    np.fill_diagonal(D, 0.0)
    return D


def main() -> int:
    test = load_test().copy()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    ids = test["tcr_id"].to_numpy()
    ep = test["epitope_aa"].to_numpy()
    epitopes = sorted(set(ep))
    members = {e: np.where(ep == e)[0] for e in epitopes}
    lab = labels.set_index("tcr_id")
    y = pd.Series(ids).map(lab["pos"]).to_numpy()
    query = np.where(y == 1)[0]
    a2 = sorted({e for e, h in zip(ep, test["hla_short"]) if h == "A2"})

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = cached["vectors"]
    W = W / np.linalg.norm(W, axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(cached["tcr_id"])}
    X = W[[po[t] for t in ids]]
    D_sceptr = 1.0 - X @ X.T
    np.fill_diagonal(D_sceptr, 0.0)
    _, fam = connected_components(sp.csr_matrix(D_sceptr <= FAMILY_RADIUS),
                                  directed=False)

    # ---------- the explicit feature spaces ----------
    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    la = test["cdr3_alpha_aa"].str.len().to_numpy(float)
    lb = test["cdr3_beta_aa"].str.len().to_numpy(float)

    def comp(col):
        s = test[col].fillna("")
        return np.array([[r.count(a) / max(len(r), 1) for a in AA] for r in s])

    def onehot(cols):
        return pd.get_dummies(test[cols].astype(str)).to_numpy(float)

    spaces = {
        "CDR3 lengths only (2 dim)": np.column_stack([la, lb]),
        "OLGA Pgen only (2 dim)": pgen.reindex(ids)[
            ["log10_pgen_alpha", "log10_pgen_beta"]].to_numpy(float),
        "lengths + Pgen (4 dim)": np.column_stack(
            [la, lb, pgen.reindex(ids)[["log10_pgen_alpha",
                                        "log10_pgen_beta"]].to_numpy(float)]),
        "AA composition (40 dim)": np.hstack([comp("cdr3_alpha_aa"),
                                              comp("cdr3_beta_aa")]),
        "V/J gene usage (one-hot)": onehot(["TRAV_IMGT", "TRAJ_IMGT",
                                            "TRBV_IMGT", "TRBJ_IMGT"]),
        "all of the above": np.hstack([
            np.column_stack([la, lb]),
            pgen.reindex(ids)[["log10_pgen_alpha",
                               "log10_pgen_beta"]].to_numpy(float),
            comp("cdr3_alpha_aa"), comp("cdr3_beta_aa"),
            onehot(["TRAV_IMGT", "TRAJ_IMGT", "TRBV_IMGT", "TRBJ_IMGT"])]),
    }

    def switch(Dm, rng):
        q = query[np.isin(ep[query], a2)]
        true_ep = ep[q]
        hits, nulls = [], []
        for _ in range(REPEATS):
            pool = {e: rng.choice(members[e],
                                  size=min(N_SUBSAMPLE, len(members[e])),
                                  replace=False) for e in a2}
            Z = np.zeros((len(q), len(a2)))
            for j, e in enumerate(a2):
                ref = pool[e]
                dqr = Dm[np.ix_(q, ref)]
                keep = (q[:, None] != ref[None, :]) \
                    & (fam[q][:, None] != fam[ref][None, :])
                np.put_along_axis(keep, np.argsort(dqr, axis=1)[:, :DROP_NN],
                                  False, axis=1)
                k = np.exp(-dqr / SIGMA) * keep
                n = keep.sum(axis=1)
                Z[:, j] = np.where(n > 0, k.sum(axis=1) / np.maximum(n, 1), np.nan)
            Z = (Z - np.nanmean(Z, 0)) / (np.nanstd(Z, 0) + 1e-12)
            picked = np.array(a2)[np.argmax(Z, axis=1)]
            hits.append(float((picked == true_ep).mean()))
            nulls.append(float((picked == rng.permutation(true_ep)).mean()))
        return float(np.mean(hits)), float(np.mean(nulls))

    print("=== what recovers E88's 2.6x? HLA-A2 panel, family excluded, m=20 ===")
    print(f"  {'feature space':<30}{'dim':>5}{'switch acc':>12}{'null':>8}{'ratio':>8}")
    rows = []
    rng = np.random.default_rng(SEED)
    acc, null = switch(D_sceptr, rng)
    print(f"  {'SCEPTR (reference, E88)':<30}{X.shape[1]:5d}{acc:12.3f}{null:8.3f}"
          f"{acc/null:8.2f}")
    rows.append(["SCEPTR (reference, E88)", X.shape[1], acc, null, acc / null])
    ref_ratio = acc / null
    for name, F in spaces.items():
        F = np.nan_to_num(F, nan=np.nanmean(F) if np.isnan(F).any() else 0.0)
        rng = np.random.default_rng(SEED)
        a, n = switch(zdist(F), rng)
        print(f"  {name:<30}{F.shape[1]:5d}{a:12.3f}{n:8.3f}{a/n:8.2f}"
              f"   {(a/n - 1) / (ref_ratio - 1):6.0%} of SCEPTR's excess")
        rows.append([name, F.shape[1], a, n, a / n])

    lines = ["# E95 — What is the antigen-naive model using?", "",
             "E89 showed a SCEPTR trained purely on OLGA sequences separates the",
             "cohorts slightly *better* than the real model (2.67x vs 2.60x), so",
             "E88's effect is sequence-space geometry. This asks along which axis,",
             "by replacing the learned representation with explicit features and",
             "re-running E88's pipeline unchanged.",
             "Generated by `scripts/85_what_the_naive_model_uses.py`.", "",
             "| feature space | dim | switch accuracy | null | ratio | share of SCEPTR's excess |",
             "|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        share = "—" if r[0].startswith("SCEPTR") else \
            f"{(r[4] - 1) / (ref_ratio - 1):.0%}"
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:.3f} | {r[3]:.3f} "
                     f"| **{r[4]:.2f}** | {share} |")
    lines.append("")

    # ---------- the blunt view: per-epitope means ----------
    print("\n=== the blunt view: per-epitope means of each quantity ===")
    tab = pd.DataFrame({
        "epitope": ep, "CDR3a len": la, "CDR3b len": lb,
        "log10 Pgen a": pgen.reindex(ids)["log10_pgen_alpha"].to_numpy(),
        "log10 Pgen b": pgen.reindex(ids)["log10_pgen_beta"].to_numpy(),
    }).groupby("epitope").agg(["mean", "std"])
    sp_of = test.drop_duplicates("epitope_aa").set_index("epitope_aa")["split"]
    print(f"  {'epitope':<12}{'split':<9}{'CDR3a':>14}{'CDR3b':>14}"
          f"{'log10 Pgen a':>16}{'log10 Pgen b':>16}")
    blunt = []
    for e in epitopes:
        r = tab.loc[e]
        print(f"  {e:<12}{sp_of[e]:<9}"
              f"{r[('CDR3a len','mean')]:8.2f}±{r[('CDR3a len','std')]:<5.2f}"
              f"{r[('CDR3b len','mean')]:8.2f}±{r[('CDR3b len','std')]:<5.2f}"
              f"{r[('log10 Pgen a','mean')]:10.2f}±{r[('log10 Pgen a','std')]:<5.2f}"
              f"{r[('log10 Pgen b','mean')]:10.2f}±{r[('log10 Pgen b','std')]:<5.2f}")
        blunt.append([e, sp_of[e], r[('CDR3a len','mean')], r[('CDR3a len','std')],
                      r[('CDR3b len','mean')], r[('CDR3b len','std')],
                      r[('log10 Pgen a','mean')], r[('log10 Pgen b','mean')]])
    # how separable is a single scalar? one-way effect size across the A2 panel
    from scipy import stats as st
    print("\n  one-way ANOVA across the six HLA-A2 cohorts (is a single scalar enough?)")
    for nm, v in (("CDR3a length", la), ("CDR3b length", lb),
                  ("log10 Pgen alpha",
                   pgen.reindex(ids)["log10_pgen_alpha"].to_numpy()),
                  ("log10 Pgen beta",
                   pgen.reindex(ids)["log10_pgen_beta"].to_numpy())):
        groups = [v[members[e]] for e in a2]
        F, p = st.f_oneway(*groups)
        grand = np.concatenate(groups)
        ss_b = sum(len(g) * (g.mean() - grand.mean()) ** 2 for g in groups)
        eta = ss_b / ((grand - grand.mean()) ** 2).sum()
        print(f"    {nm:<18}F {F:8.2f}  p {p:.2e}   eta^2 {eta:.4f}")

    lines += ["## The blunt view: per-epitope means", "",
              "| epitope | split | CDR3α len | CDR3β len | log10 Pgen α | log10 Pgen β |",
              "|---|---|---:|---:|---:|---:|"]
    for b in blunt:
        lines.append(f"| {b[0]} | {b[1]} | {b[2]:.2f} ± {b[3]:.2f} "
                     f"| {b[4]:.2f} ± {b[5]:.2f} | {b[6]:.2f} | {b[7]:.2f} |")
    lines.append("")

    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
