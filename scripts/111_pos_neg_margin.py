#!/usr/bin/env python
"""E126 -- Antigen-conditioned positive/negative boundary (Minervina exo-TCR).

The collaborator's M_e(t) = d(t, N_e^-) - d(t, N_e^+): positive means closer to
measured binders than to measured non-binders. D reads only the positive side,
so the negative side is, on its face, a channel D does not have.

Panels available from Zenodo 10.5281/zenodo.14885005 (`pos_neg.zip`), which
covers two of the three private epitopes and three public ones for validation:

  LLWa  LLWNGPMAV   private   beta partners   24,595 pos / 336,356 neg
  LTDb  LTDEMIAQY   private   alpha partners   3,549 pos / 649,373 neg
  GILb  GILGFVFTL   public    alpha partners   4,019 pos / 448,477 neg
  TTDb  TTDPSFLGRY  public    alpha partners     533 pos / 377,252 neg
  YLQb  YLQPRTFLL   public    alpha partners     828 pos / 379,657 neg

LTDEMIAQY matters most: its competition reference pool holds **5 rows**, which
E84 established is far below the ~60 where the pool term does anything, and
here there are 3,549 experimentally validated alpha chains for it.

Two design points that decide whether the number means anything:

**Size matching.** Nearest-neighbour distance falls as a reference set grows,
and the negative sets are 80-1200x larger than the positive ones. An unmatched
margin would measure set size, not biology. Positives and negatives are
therefore subsampled to a common n per epitope, over several seeds.

**The typicality confound.** A 649k-sequence "true negative" set is close to
the four donors' whole repertoire minus the binders, so "my CDR3 appears in
the negative set" largely means "my CDR3 is a common sequence in these donors".
That is [[typicality]], the gravity well that has already absorbed four
features in this project. The standing typicality gate is the check, and the
exact-match signal is reported beside the margin so the two can be told apart.

Single chains against a fixed exo-TCR anchor, matched to natural alphabeta
competition receptors on one chain only -- so a match is not a receptor match.

Run: .venv/bin/python scripts/111_pos_neg_margin.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pwseqdist as pw
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

PANEL = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
             "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/pn/pos_neg")
REPORT = ROOT / "reports/89_pos_neg_margin.md"
UC = ["CDR3.amino.acid.sequence.x", "bestVGene.x", "donor.x", "frac"]
PANELS = [("LLWa", "LLWNGPMAV", "cdr3_beta_aa"),
          ("LTDb", "LTDEMIAQY", "cdr3_alpha_aa"),
          ("GILb", "GILGFVFTL", "cdr3_alpha_aa"),
          ("TTDb", "TTDPSFLGRY", "cdr3_alpha_aa"),
          ("YLQb", "YLQPRTFLL", "cdr3_alpha_aa")]
SEEDS = (0, 1, 2)
KNN = 5
MAXN = 4000          # cap per side, for runtime
# Pre-registered, as reports/87_candidate_admission.md requires. Floor 0.5098,
# best ever 0.5606. A real negative-side channel on 2/3 private cohorts should
# land clearly above the best-ever if it is new; predict 0.56-0.62.
PREDICTED_RESIDUAL = (0.56, 0.62)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def knn_dist(queries: list[str], refs: list[str], k: int) -> np.ndarray:
    """Mean distance to the k nearest references, tcrdist on CDR3."""
    d = pw.apply_pairwise_rect(metric=pw.metrics.nb_vector_tcrdist,
                               seqs1=queries, seqs2=refs,
                               use_numba=True, ncpus=1)
    k = min(k, d.shape[1])
    return np.sort(d, axis=1)[:, :k].mean(axis=1)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    lab = scored.set_index("tcr_id")["pos"]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)

    margin = pd.Series(index=test["tcr_id"], dtype=float)
    exact = pd.Series(0.0, index=test["tcr_id"])
    rows = []
    print(f"=== margin M = d(t,neg) - d(t,pos), size-matched, k={KNN}, "
          f"{len(SEEDS)} seeds ===")
    print(f"  {'epitope':<12}{'split':<9}{'n':>5}{'matched n':>10}"
          f"{'AUC margin':>12}{'AUC exact':>11}{'base':>7}")
    for tag, ep, ccol in PANELS:
        p = pd.read_csv(PANEL / f"{tag}_true_pos.tsv", sep="\t",
                        usecols=UC, low_memory=False)
        n = pd.read_csv(PANEL / f"{tag}_true_neg.tsv", sep="\t",
                        usecols=UC, low_memory=False)
        ps = p["CDR3.amino.acid.sequence.x"].dropna().unique().tolist()
        ns = n["CDR3.amino.acid.sequence.x"].dropna().unique().tolist()
        comp = test[test.epitope_aa == ep]
        q = comp[ccol].fillna("").tolist()
        keep = [i for i, s in enumerate(q) if len(s) >= 5]
        qq = [q[i] for i in keep]
        ids = comp["tcr_id"].to_numpy()[keep]

        m = min(len(ps), len(ns), MAXN)
        acc = np.zeros(len(qq))
        for sd in SEEDS:
            rng = np.random.default_rng(sd)
            pa = list(rng.choice(ps, m, replace=False))
            na = list(rng.choice(ns, m, replace=False))
            acc += knn_dist(qq, na, KNN) - knn_dist(qq, pa, KNN)
        margin.loc[ids] = acc / len(SEEDS)

        # exact-match signal, for comparison with the margin
        sp, sn = set(ps), set(ns)
        cq = comp[ccol]
        ip, inn = cq.isin(sp).to_numpy(), cq.isin(sn).to_numpy()
        exact.loc[comp["tcr_id"]] = np.where(ip & ~inn, 1.0,
                                             np.where(inn & ~ip, -1.0, 0.0))

        vis = [t for t in ids if t in lab.index]
        am = ae = bs = np.nan
        if len(vis) > 10:
            y = lab.reindex(vis).to_numpy()
            if len(set(y)) > 1:
                am = roc_auc_score(y, margin.reindex(vis).to_numpy())
                ae = roc_auc_score(y, exact.reindex(vis).to_numpy())
                bs = y.mean()
        rows.append([ep, comp["split"].iloc[0], len(comp), m, am, ae, bs])
        f = lambda v: "   n/a" if np.isnan(v) else f"{v:.4f}"  # noqa: E731
        print(f"  {ep:<12}{comp['split'].iloc[0]:<9}{len(comp):5d}{m:10d}"
              f"{f(am):>12}{f(ae):>11}"
              f"{'  n/a' if np.isnan(bs) else f'{bs:7.3f}'}")

    pub = sorted(scored["epitope_aa"].unique())
    mm = {"margin": margin, "exact match": exact}
    print("\n=== gates ===")
    print(f"  pre-registered residual-on-D: {PREDICTED_RESIDUAL[0]:.2f}"
          f"-{PREDICTED_RESIDUAL[1]:.2f}")
    gate_rows = []
    for nm, s in mm.items():
        filled = s.copy()
        for ep, g in test.groupby("epitope_aa"):
            i = g["tcr_id"]
            v = filled.reindex(i)
            filled.loc[i] = v.fillna(v.mean() if v.notna().any() else 0.0).to_numpy()
        tg = typicality_gate(test, scored, filled)
        pv, _ = private_variance_gate(test, filled)
        macro = hyp.mean_auc(labels, filled)
        res = pd.Series(index=filled.index, dtype=float)
        leak = []
        for ep, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in filled.index]
            x = dr.reindex(i).rank(pct=True).to_numpy()
            y = filled.reindex(i).rank(pct=True).to_numpy()
            r = y - np.polyval(np.polyfit(x, y, 1), x)
            res.loc[i] = r
            leak.append(abs(spearmanr(r, x).statistic))
        rr = hyp.mean_auc(labels, res)
        resid, mx = max(rr, 1 - rr), max(leak)
        rho = tg["rho_with_typicality"]
        print(f"  {nm:<13} macro {macro:.4f}  typ rho {rho:+.3f} "
              f"{'pass' if abs(rho) < 0.6378 else 'FAIL':<5} "
              f" priv {'pass' if pv else 'FAIL'}"
              f"  resid {resid:.4f} leak {mx:.3f} "
              f"{'VALID' if mx < 0.1 else 'INVALID'}")
        gate_rows.append([nm, macro, rho, pv, resid, mx])

        sr = hyp.rank_within(filled, test)
        WS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
        per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
        vals = {w: float(np.mean([per[w][e] for e in pub])) for w in WS}
        base = vals[0.0]
        held = []
        for e in pub:
            trn = [x for x in pub if x != e]
            pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
            held.append(per[pick][e])
        nested = float(np.mean(held))
        bw = max(WS, key=lambda w: vals[w])
        print(f"  {'':<13} blend: oracle w={bw:.2f} {vals[bw]:.4f} "
              f"({vals[bw]-base:+.4f})   nested {nested:.4f} ({nested-base:+.4f})")
        gate_rows[-1] += [bw, vals[bw] - base, nested - base]

    lines = ["# E126 — Antigen-conditioned positive/negative boundary", "",
             "Generated by `scripts/111_pos_neg_margin.py`.", "",
             "## Per-epitope, where labels exist", "",
             "| epitope | split | rows | matched n/side | AUC margin | AUC exact | base rate |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        f = lambda v: "n/a" if np.isnan(v) else f"{v:.4f}"  # noqa: E731
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {f(r[4])} "
                     f"| {f(r[5])} | {f(r[6])} |")
    lines += ["", "## Gates", "",
              "| feature | macro | typicality ρ | private var | residual on D | leak | oracle w | oracle Δ | nested Δ |",
              "|---|---:|---:|---|---:|---:|---:|---:|---:|"]
    for g in gate_rows:
        lines.append(f"| {g[0]} | {g[1]:.4f} | {g[2]:+.3f} "
                     f"| {'pass' if g[3] else 'FAIL'} | **{g[4]:.4f}** "
                     f"| {g[5]:.3f} | {g[6]:.2f} | {g[7]:+.4f} | **{g[8]:+.4f}** |")
    lines += ["", f"Pre-registered residual prediction: "
              f"{PREDICTED_RESIDUAL[0]:.2f}–{PREDICTED_RESIDUAL[1]:.2f}.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
