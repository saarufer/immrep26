#!/usr/bin/env python
"""E138 -- Selection from generation to the mature repertoire.

A collaborator's proposal, the one of their four that passes all four admission
questions: OTS is a channel D does not read, the quantity is per-TCR so every
row including all 602 private ones gets a value, and it is independent of the
labelling process.

**The naive stratum does not exist, so the quantity is renamed.** The proposal
asks for Q = log P_observed_naive - log P_gen and warns, correctly, that not all
OTS records are naive. Checking OTS's own metadata: all **384** `includes_naive`
units are COVID patients, and the healthy + no-cancer + untreated subset
contains **zero** naive-labelled units. There is no healthy naive stratum in
OTS. What is constructible is a resting healthy *peripheral* repertoire -- 215
units, 425,409 sequences, after dropping the deliberately stimulated and
functionally sorted strata (il12/il18/cd28-stimulated, antigen-stimulated,
cd154+, tetramer-sorted) which would import antigen-driven expansion directly.

So the measured quantity is

    Q(t) = log P_healthy_peripheral(t) - log P_gen(t)

a selection factor over the healthy *mature* repertoire, not the naive one.
That is a weaker claim than proposed and it is stated as such. It also carries a
residual risk worth naming: a mature repertoire is partly shaped by past
exposure to common antigens, and four of the five public epitopes are from CMV,
EBV, influenza and SARS-CoV-2. That is not circular with the VDJdb annotation
process, so it is not a Q3 violation, but it is a route by which the feature
could carry real signal for a reason different from thymic selection.

**Estimator.** A full paired alphabeta sequence probability is far too sparse,
so the selection factor is estimated the way SONIA does, as a product of
feature marginals -- V gene, J gene, CDR3 length, and (position, amino acid)
indexed from both ends -- with p_gen from a large OLGA sample and
p_repertoire from the OTS subset. Prior art is explicit (SONIA, repertoire
selection models); "selection exists after generation" is not the claim. The
claim under test is narrower: does mature-repertoire selection explain
revalidation residual beyond generation probability and typicality?

**Pre-registered residual-on-D prediction: 0.51-0.55**, i.e. near the floor.
The reasoning is on record before the run: Q is built from compositional
marginals, and E135 has just shown that a positional-composition null recovers
**95.6%** of the real density's above-chance signal, so composition is already
almost entirely inside D. E78's generation-background correction also failed.
Floor is 0.5098 and the project record is 0.5609.

Run: .venv/bin/python scripts/119_repertoire_selection.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

OTS = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
           "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/ots_healthy.tsv")
REPORT = ROOT / "reports/97_repertoire_selection.md"
N_GEN = 400_000
PSEUDO = 1.0
SEED = 0
PREDICTED_RESIDUAL = (0.51, 0.55)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def gene(g):
    if not isinstance(g, str) or not g:
        return None
    return g.split(",")[0].split("*")[0].strip() or None


def marginals(V, J, C):
    """SONIA-style feature counts: V, J, length, (pos,aa) from both ends."""
    m = {"V": Counter(), "J": Counter(), "L": Counter(), "PA": Counter()}
    for v, j, c in zip(V, J, C):
        if not isinstance(c, str) or len(c) < 4:
            continue
        if v:
            m["V"][v] += 1
        if j:
            m["J"][j] += 1
        m["L"][len(c)] += 1
        for i, aa in enumerate(c):
            m["PA"][(i, aa)] += 1
            m["PA"][(i - len(c), aa)] += 1
    return m


def logratio(m_rep, m_gen):
    """log p_rep(f) - log p_gen(f) per feature, with add-one smoothing."""
    out = {}
    for k in m_rep:
        tr = sum(m_rep[k].values()) + PSEUDO * len(set(m_rep[k]) | set(m_gen[k]))
        tg = sum(m_gen[k].values()) + PSEUDO * len(set(m_rep[k]) | set(m_gen[k]))
        keys = set(m_rep[k]) | set(m_gen[k])
        out[k] = {f: np.log((m_rep[k][f] + PSEUDO) / tr)
                  - np.log((m_gen[k][f] + PSEUDO) / tg) for f in keys}
    return out


def score(V, J, C, LR):
    vals = []
    for v, j, c in zip(V, J, C):
        if not isinstance(c, str) or len(c) < 4:
            vals.append(np.nan)
            continue
        q = 0.0
        q += LR["V"].get(v, 0.0) if v else 0.0
        q += LR["J"].get(j, 0.0) if j else 0.0
        q += LR["L"].get(len(c), 0.0)
        for i, aa in enumerate(c):
            q += LR["PA"].get((i, aa), 0.0) + LR["PA"].get((i - len(c), aa), 0.0)
        vals.append(q)
    return np.array(vals)


def main() -> int:
    from tcrdist.pgen import OlgaModel

    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    print("=== the repertoire reference ===")
    R = pd.read_csv(OTS, sep="\t", low_memory=False).dropna(
        subset=["junction_aa_alpha", "junction_aa_beta"])
    for ch in ("alpha", "beta"):
        R[f"v_{ch}"] = R[f"v_call_{ch}"].map(gene)
        R[f"j_{ch}"] = R[f"j_call_{ch}"].map(gene)
    R = R.drop_duplicates(subset=["junction_aa_alpha", "junction_aa_beta"])
    print(f"  OTS healthy resting peripheral, deduplicated: {len(R):,} paired receptors")
    print(f"  (no healthy NAIVE stratum exists in OTS: all 384 'includes_naive'")
    print(f"   units are COVID patients, so this is a MATURE repertoire)")

    print("\n=== the generative background ===")
    rng = np.random.default_rng(SEED)
    gen = {}
    for ch, folder, rec in (("alpha", "human_T_alpha", "VJ"),
                            ("beta", "human_T_beta", "VDJ")):
        m = OlgaModel(chain_folder=folder, recomb_type=rec)
        seqs = []
        while len(seqs) < N_GEN:
            seqs += [s for s in m.gen_cdr3s(n=100_000) if s]
        gen[ch] = seqs[:N_GEN]
        print(f"  {ch}: {len(gen[ch]):,} OLGA CDR3s")
    print("  note: OLGA emits CDR3s without V/J, so the V and J marginals are")
    print("  taken from the model's own germline usage distribution instead.")
    gv = {}
    for ch, folder, rec in (("alpha", "human_T_alpha", "VJ"),
                            ("beta", "human_T_beta", "VDJ")):
        m = OlgaModel(chain_folder=folder, recomb_type=rec)
        gm = m.generative_model
        gd = m.genomic_data
        pv = np.asarray(getattr(gm, "PV", None) if hasattr(gm, "PV") else
                        gm.PVJ.sum(axis=1) if hasattr(gm, "PVJ") else None,
                        dtype=float)
        pj = np.asarray(gm.PVJ.sum(axis=0) if hasattr(gm, "PVJ") else
                        gm.PDJ.sum(axis=0), dtype=float)
        vn = [gene(x[0]) for x in gd.genV]
        jn = [gene(x[0]) for x in gd.genJ]
        gv[ch] = (dict(zip(vn, pv / pv.sum())), dict(zip(jn, pj / pj.sum())))
        print(f"  {ch}: {len(vn)} V genes, {len(jn)} J genes from the OLGA model")

    print("\n=== building the selection factor ===")
    Q = {}
    for ch in ("alpha", "beta"):
        m_rep = marginals(R[f"v_{ch}"], R[f"j_{ch}"], R[f"junction_aa_{ch}"])
        m_gen = marginals([None] * len(gen[ch]), [None] * len(gen[ch]), gen[ch])
        # germline marginals from the OLGA model, scaled to the same total
        tot = sum(m_rep["V"].values())
        for f, p in gv[ch][0].items():
            if f:
                m_gen["V"][f] = p * tot
        for f, p in gv[ch][1].items():
            if f:
                m_gen["J"][f] = p * tot
        LR = logratio(m_rep, m_gen)
        col = "cdr3_alpha_aa" if ch == "alpha" else "cdr3_beta_aa"
        vc = "TRAV_IMGT" if ch == "alpha" else "TRBV_IMGT"
        jc = "TRAJ_IMGT" if ch == "alpha" else "TRBJ_IMGT"
        Q[ch] = score(test[vc].map(gene), test[jc].map(gene), test[col], LR)
        print(f"  Q_{ch}: mean {np.nanmean(Q[ch]):+.2f}  sd {np.nanstd(Q[ch]):.2f}"
              f"  missing {int(np.isnan(Q[ch]).sum())}")
    ids = test["tcr_id"]
    F = {"Q_alpha": pd.Series(Q["alpha"], index=ids),
         "Q_beta": pd.Series(Q["beta"], index=ids),
         "Q_paired": pd.Series(Q["alpha"] + Q["beta"], index=ids)}
    # the control the proposal is measured against
    pg = ROOT / "data/interim/pgen.csv"
    if pg.exists():
        P = pd.read_csv(pg).set_index("tcr_id")
        for c in P.columns:
            if "log" in c.lower() or "pgen" in c.lower():
                F[f"control {c}"] = P[c].reindex(ids)

    print("\n=== gates, with the E112 leak check ===")
    print(f"  pre-registered residual for Q_paired: {PREDICTED_RESIDUAL[0]:.2f}"
          f"-{PREDICTED_RESIDUAL[1]:.2f}  (floor 0.5098, record 0.5609)")
    print(f"  {'quantity':<26}{'macro':>8}{'typ rho':>9}{'typ':>6}{'priv':>6}"
          f"{'resid':>8}{'leak':>7}{'valid':>7}")
    rows = []
    for nm, s in F.items():
        v = s.astype(float).copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.median()).to_numpy()
        v = v.fillna(v.median())
        tg = typicality_gate(test, scored, v)
        pv, _ = private_variance_gate(test, v)
        m = hyp.mean_auc(labels, v)
        res = pd.Series(index=v.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in v.index]
            x = dr.reindex(i).rank(pct=True).to_numpy()
            y = v.reindex(i).rank(pct=True).to_numpy()
            r = y - np.polyval(np.polyfit(x, y, 1), x)
            res.loc[i] = r
            leak.append(abs(spearmanr(r, x).statistic))
        rr = hyp.mean_auc(labels, res)
        resid, mx = max(rr, 1 - rr), max(leak)
        rho = tg["rho_with_typicality"]
        print(f"  {nm:<26}{max(m,1-m):8.4f}{rho:+9.3f}"
              f"{'pass' if abs(rho)<0.6378 else 'FAIL':>6}"
              f"{'pass' if pv else 'FAIL':>6}{resid:8.4f}{mx:7.3f}"
              f"{'yes' if mx<0.1 else 'NO':>7}")
        rows.append([nm, max(m, 1-m), rho, pv, resid, mx])

    print("\n=== blended into D, five-epitope LOEO ===")
    pub = sorted(scored["epitope_aa"].unique())
    base = hyp.mean_auc(labels, dr)
    WS = [0.0, 0.05, 0.10, 0.15, 0.20]
    print(f"  D alone {base:.4f}")
    blends = []
    for nm, s in F.items():
        v = s.astype(float).copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.median()).to_numpy()
        v = v.fillna(v.median())
        best = None
        for sign in (1.0, -1.0):
            sr = hyp.rank_within(sign * v, test)
            per = {w: hyp.per_epitope(labels, (1-w)*dr + w*sr) for w in WS}
            vals = {w: float(np.mean([per[w][e] for e in pub])) for w in WS}
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pick][e])
            nested = float(np.mean(held))
            bw = max(WS, key=lambda w: vals[w])
            if best is None or vals[bw] > best[2]:
                best = (sign, bw, vals[bw], nested)
        sign, bw, ov, nested = best
        print(f"  {nm:<26}sign {int(sign):+d}  oracle w={bw:.2f} {ov:.4f} "
              f"({ov-base:+.4f})   nested {nested:.4f} ({nested-base:+.4f})")
        blends.append([nm, sign, bw, ov-base, nested-base])

    lines = ["# E138 — Generation-to-mature-repertoire selection", "",
             "Generated by `scripts/119_repertoire_selection.py`.", "",
             "**The naive stratum does not exist in OTS.** All 384",
             "`includes_naive` units are COVID patients and the healthy +",
             "no-cancer + untreated subset contains zero naive-labelled units, so",
             "the quantity measured is",
             "`Q(t) = log P_healthy_peripheral(t) − log P_gen(t)` — a **mature**",
             "repertoire selection factor, not a naive one.", "",
             f"Reference: {len(R):,} deduplicated paired receptors from 215 OTS",
             "units (healthy, no cancer, untreated, stimulated and functionally",
             "sorted strata dropped). Background: "
             f"{N_GEN:,} OLGA CDR3s per chain plus the model's own V/J usage.", "",
             f"Pre-registered residual-on-D prediction: "
             f"{PREDICTED_RESIDUAL[0]:.2f}–{PREDICTED_RESIDUAL[1]:.2f}.", "",
             "| quantity | macro | typicality ρ | typ gate | private var | residual on D | leak | valid |",
             "|---|---:|---:|---|---|---:|---:|---|"]
    for nm, m, rho, pv, resid, mx in rows:
        lines.append(f"| `{nm}` | {m:.4f} | {rho:+.3f} "
                     f"| {'pass' if abs(rho)<0.6378 else 'FAIL'} "
                     f"| {'pass' if pv else 'FAIL'} | **{resid:.4f}** "
                     f"| {mx:.3f} | {'yes' if mx<0.1 else 'NO'} |")
    lines += ["", "| quantity | sign | oracle w | oracle Δ | nested Δ |",
              "|---|---:|---:|---:|---:|"]
    for nm, sg, bw, od, nd in blends:
        lines.append(f"| `{nm}` | {int(sg):+d} | {bw:.2f} | {od:+.4f} | **{nd:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
