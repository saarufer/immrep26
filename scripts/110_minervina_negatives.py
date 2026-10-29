#!/usr/bin/env python
"""E125 -- Real experimental negatives from Minervina's exo-TCR panels.

The collaborator's E125 proposal, narrowed to the one part that can touch the
private split. Of the sources proposed:

  P-TEAM (NLV single-mutant panel)  NLVPMVATV is PUBLIC, and the authors'
                                    own caveat applies: a TCR selected for
                                    NLV that fails to respond to a mutant is
                                    not an NLV negative. No Kaggle value.
  Atlas Yellow Fever variants       the 44,431-row Atlas table already
                                    fetched for E122 has NO binding-label
                                    column at all -- it is positives only,
                                    and 67% of it is vdjdb_2025_12_29.
  Gaevert (mouse H-2Db)             not a competition MHC.
  Minervina exo-TCR                 GILGFVFTL (public, 670 rows) and
                                    LLWNGPMAV (PRIVATE, 229 rows).  <-- this

LLWNGPMAV is the yellow fever NS4B(214-222) epitope, which is why Minervina's
YF17D work covers it. It is also the one private epitope whose positive pool is
alive (191 rows; CIN has 11, LTD has 5), so adding positives there is saturated
-- E122 took it 191->203 for nothing. An experimentally measured *negative*
reference is a channel D does not read at all.

Admission standard (reports/87_candidate_admission.md):

  Q1 channel      real functional/binding negatives for the target pHLA. D
                  reads only positive pools; "far from known binders" and
                  "close to known non-binders" are different statements.
                  Predicted residual-on-D: stated below before it is computed.
  Q2 variance     measured per cohort here, before any weight is chosen.
  Q3 independence audited: the shared LLW sequences entered VDJdb from
                  PMID:28103239 (2017) and PMID:32081129 (2020), not from
                  Minervina 2025. Independent assay, overlapping donors.
  Q4 constraint   information, not noise.

The hard limitation, stated up front: these are SINGLE chains measured against
a fixed exo-TCR anchor (LLWa = beta partners on a fixed alpha; GILb = alpha
partners on a fixed beta). Competition rows are natural alphabeta pairs.
Matching on one chain is not matching a receptor, and the two epitopes use
opposite anchors, so GIL does not cleanly validate LLW.

Run: .venv/bin/python scripts/110_minervina_negatives.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

SCRATCH = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
               "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/"
               "cis/close_in_seq_pos_neg")
REPORT = ROOT / "reports/88_minervina_negatives.md"
PANELS = [("LLWa", "LLWNGPMAV", "cdr3b", "cdr3_beta_aa"),
          ("GILb", "GILGFVFTL", "cdr3a", "cdr3_alpha_aa")]
# Stated BEFORE running, as the standard requires. Floor 0.5098, best ever
# 0.5606. A 78%-concordant proxy covering a third of a cohort should, if it is
# real and not already in D, land high -- predict 0.58-0.62 on GIL.
PREDICTED_RESIDUAL = (0.58, 0.62)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def panel(tag: str) -> pd.DataFrame:
    p = pd.read_csv(SCRATCH / f"{tag}_close_positives.csv")
    n = pd.read_csv(SCRATCH / f"{tag}_close_negatives.csv")
    return pd.concat([p, n], ignore_index=True)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    lab = labels.dropna(subset=["pos"]).set_index("tcr_id")["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)

    print("=== 1. coverage and label concordance ===")
    feat = pd.Series(index=test["tcr_id"], dtype=float)
    cov_rows = []
    for tag, ep, chain, ccol in PANELS:
        d = panel(tag)
        frac = d.groupby(chain)["positive"].mean()
        nrec = d.groupby(chain).size()
        comp = test[test.epitope_aa == ep]
        hit = comp[comp[ccol].isin(frac.index)]
        feat.loc[hit["tcr_id"]] = hit[ccol].map(frac).to_numpy()
        vis = hit[hit.tcr_id.isin(lab.index)]
        agree = np.nan
        if len(vis):
            a = lab.reindex(vis.tcr_id).to_numpy()
            b = (vis[ccol].map(frac) > 0.5).astype(int).to_numpy()
            agree = float((a == b).mean())
        cov_rows.append([ep, comp["split"].iloc[0], len(comp), len(hit),
                         len(hit) / len(comp), len(vis), agree,
                         float(nrec.reindex(hit[ccol]).mean())])
        print(f"  {ep:<12}{comp['split'].iloc[0]:<9}{len(comp):5d} rows"
              f"  matched {len(hit):4d} ({len(hit)/len(comp)*100:5.1f}%)"
              f"  labelled {len(vis):4d}"
              f"  agreement {'n/a' if np.isnan(agree) else f'{agree*100:.1f}%'}")

    print("\n=== 2. Q2: within-cohort variance of the feature, before fitting ===")
    for tag, ep, chain, ccol in PANELS:
        s = feat.reindex(test[test.epitope_aa == ep]["tcr_id"]).dropna()
        print(f"  {ep:<12}n={len(s):4d}  mean {s.mean():.3f}  sd {s.std():.3f}"
              f"  distinct {s.nunique():4d}  at 0 {int((s==0).sum()):4d}"
              f"  at 1 {int((s==1).sum()):4d}")

    print("\n=== 3. the feature alone, on the cohort where labels exist ===")
    gil = test[test.epitope_aa == "GILGFVFTL"]
    v = feat.reindex(gil["tcr_id"]).dropna()
    v = v[v.index.isin(lab.index)]
    y = lab.reindex(v.index)
    from sklearn.metrics import roc_auc_score
    print(f"  GILGFVFTL, matched+labelled rows only (n={len(v)}, "
          f"{int(y.sum())} pos): AUC {roc_auc_score(y, v):.4f}")
    gl = lab.reindex(gil["tcr_id"]).dropna()
    print(f"  for reference, D on the same {len(v)} rows: "
          f"{roc_auc_score(y, dr.reindex(v.index)):.4f}")
    print(f"  D on the whole GILGFVFTL cohort ({len(gl)} labelled): "
          f"{roc_auc_score(gl, dr.reindex(gl.index)):.4f}")

    print("\n=== 4. gates, with the E112 leak check ===")
    print(f"  predicted residual-on-D, pre-registered: "
          f"{PREDICTED_RESIDUAL[0]:.2f}-{PREDICTED_RESIDUAL[1]:.2f}")
    # fill unmatched with the cohort mean so the feature is defined everywhere
    filled = feat.copy()
    for ep, g in test.groupby("epitope_aa"):
        ids = g["tcr_id"]
        m = filled.reindex(ids)
        filled.loc[ids] = m.fillna(m.mean() if m.notna().any() else 0.5).to_numpy()
    tg = typicality_gate(test, labels.dropna(subset=["pos"]), filled)
    pv, _ = private_variance_gate(test, filled)
    macro = hyp.mean_auc(labels, filled)
    leak = []
    res = pd.Series(index=filled.index, dtype=float)
    for ep, g in test.groupby("epitope_aa"):
        ids = [t for t in g["tcr_id"] if t in filled.index]
        x = dr.reindex(ids).rank(pct=True).to_numpy()
        yy = filled.reindex(ids).rank(pct=True).to_numpy()
        r = yy - np.polyval(np.polyfit(x, yy, 1), x)
        res.loc[ids] = r
        leak.append(abs(spearmanr(r, x).statistic))
    rm = hyp.mean_auc(labels, res)
    resid = max(rm, 1 - rm)
    print(f"  macro {macro:.4f} (inv {1-macro:.4f})   typicality rho "
          f"{tg['rho_with_typicality']:+.3f}"
          f"   private variance {'pass' if pv else 'FAIL'}")
    print(f"  residual on D {resid:.4f}   max |rho(resid,D)| {max(leak):.3f}"
          f"   {'VALID' if max(leak) < 0.1 else 'INVALID -- tie leakage'}")
    print(f"  floor 0.5098   best ever 0.5606   prediction "
          f"{'HIT' if PREDICTED_RESIDUAL[0] <= resid <= PREDICTED_RESIDUAL[1] else 'MISSED'}")

    print("\n=== 5. blended into D, nested LOEO ===")
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    base = hyp.mean_auc(labels, dr)
    WS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
    sr = hyp.rank_within(filled, test)
    per_w = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
    vals = {w: float(np.mean([per_w[w][e] for e in pub])) for w in WS}
    bw = max(WS, key=lambda w: vals[w])
    held = []
    for e in pub:
        tr = [x for x in pub if x != e]
        pick = max(WS, key=lambda w: np.mean([per_w[w][t] for t in tr]))
        held.append(per_w[pick][e])
    nested = float(np.mean(held))
    print(f"  D alone {base:.4f}")
    for w in WS:
        print(f"    w={w:.2f}  {vals[w]:.4f} ({vals[w]-base:+.4f})"
              f"   GIL {per_w[w]['GILGFVFTL']:.4f}")
    print(f"  oracle w={bw:.2f} {vals[bw]:.4f} ({vals[bw]-base:+.4f})"
          f"   nested {nested:.4f} ({nested-base:+.4f})")

    lines = ["# E125 — Real experimental negatives (Minervina exo-TCR)", "",
             "Generated by `scripts/110_minervina_negatives.py`.", "",
             "## 1. Coverage and concordance", "",
             "| epitope | split | rows | matched | % | labelled | agreement | records/seq |",
             "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in cov_rows:
        ag = "n/a" if np.isnan(r[6]) else f"{r[6]*100:.1f}%"
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]*100:.1f}% "
                     f"| {r[5]} | **{ag}** | {r[7]:.1f} |")
    lines += ["", "## 4. Gates", "",
              f"- macro **{macro:.4f}**, typicality ρ {tg['rho_with_typicality']:+.3f}, "
              f"private variance {'pass' if pv else 'FAIL'}",
              f"- residual on D **{resid:.4f}** (floor 0.5098, best ever 0.5606), "
              f"leak {max(leak):.3f}",
              f"- pre-registered prediction {PREDICTED_RESIDUAL[0]:.2f}–"
              f"{PREDICTED_RESIDUAL[1]:.2f}: "
              f"**{'HIT' if PREDICTED_RESIDUAL[0] <= resid <= PREDICTED_RESIDUAL[1] else 'MISSED'}**",
              "", "## 5. Blend", "", f"D alone **{base:.4f}**; oracle w={bw:.2f} "
              f"{vals[bw]:.4f} ({vals[bw]-base:+.4f}); nested **{nested:.4f} "
              f"({nested-base:+.4f})**.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
