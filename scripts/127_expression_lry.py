#!/usr/bin/env python
"""E149 -- TCR surface-expression framework residues (the LRY score).

An external review's highest-priority proposal: same annotation, but some TCRs
may be harder to express and assemble in the revalidation system, so they fail
for a reason that has nothing to do with recognition. Thomas et al.
(Nat Commun, PMC6773850) found three variable-framework substitutions -- the
"LRY" modification -- that raise TCR surface expression 2-6x, and the effect
persists in Jurkat cells **lacking endogenous TCR**, so it is not competition
with an endogenous receptor:

    alpha IMGT 96   Pro -> Leu
    beta  IMGT  9   Ser -> Arg
    beta  IMGT 10   Asn -> Tyr

The review also spotted a real connection: `polar_frac`'s signal localises to
the germline V region rather than CDR3, so it asked whether that is partly a
receptor-intrinsic property.

**Two structural facts established before scoring, both of which reshape the
proposal.**

1. `full_seq_alpha_aa` / `full_seq_beta_aa` are **103-121 residues** and end at
   the J region (`...GGGTRLMVKP`, `...GPGTRLTVL`). They are V domains, **not the
   expressed construct** -- there is no constant region at all. So the review's
   step 1, checking constant-region replacement and framework optimisation in
   the actual screening construct, **cannot be done from competition data**.
2. Holding the V gene fixed, framework variation is **0 of 2,636 rows (0.0%)**.
   One V gene maps to exactly one framework, 100% of the time. So any
   framework-based expression score is **exactly** a V-gene statistic, not
   approximately one.

That lands in territory measured twice: E111 established `polar_frac` is a
V-gene property but not one D already carries, and E113 found a
**label-supervised** V-gene prior transfers below chance on every fold. What
this proposal adds is a **literature-fixed** score -- the positions come from
Thomas et al., not from searching the competition labels -- which avoids
E113's Q3 problem.

**Position mapping, sanity-checked before use.** beta IMGT 9/10 taken at raw
N-terminal offsets 9/10; alpha IMGT 96 located 8 residues before the conserved
Cys at IMGT 104, i.e. the CDR3 start, which resolved for 2,636/2,636 rows. The
residue distributions match the paper's dichotomy: alpha 96 is P:1033 / L:664,
beta 9 is K:1031 / R:982 / S:490, beta 10 is Y:1182 / H:832.

**Evaluated by stratification, not residual regression.** This is a 4-level
discrete feature and the E112 tie artefact has now produced spurious
"records" five times, most recently four of them in E139. Residual-on-D is
reported with its leak check, and the stratified test inside D-quintiles is
the instrument that decides.

**Pre-registered: residual-on-D 0.51-0.54, nested -0.002 to +0.002.**

Run: .venv/bin/python scripts/127_expression_lry.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/104_expression_lry.md"
PREDICTED_RESID = (0.51, 0.54)
PREDICTED_NESTED = (-0.002, 0.002)
DOMINANT = {"a96": "L", "b9": "R", "b10": "Y"}
WEAK = {"a96": "P", "b9": "S", "b10": "N"}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    scored = labels.dropna(subset=["pos"]).copy()
    scored["pos"] = scored["pos"].astype(int)
    lab = scored.set_index("tcr_id")["pos"]
    pub = sorted(scored["epitope_aa"].unique())
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    # --- extract the three positions ---
    res = {}
    b = test["full_seq_beta_aa"]
    res["b9"] = [s[8] if isinstance(s, str) and len(s) > 8 else None for s in b]
    res["b10"] = [s[9] if isinstance(s, str) and len(s) > 9 else None for s in b]
    a, c3 = test["full_seq_alpha_aa"], test["cdr3_alpha_aa"]
    aa = []
    for s, c in zip(a, c3):
        i = s.rfind(c) if isinstance(s, str) and isinstance(c, str) else -1
        aa.append(s[i - 8] if i >= 8 else None)
    res["a96"] = aa
    R = pd.DataFrame(res, index=test["tcr_id"])

    print("=== the three LRY positions ===")
    for k in ("a96", "b9", "b10"):
        vc = R[k].value_counts()
        print(f"  {k:<5}dominant '{DOMINANT[k]}' in {int((R[k]==DOMINANT[k]).sum()):5d}"
              f"   weak '{WEAK[k]}' in {int((R[k]==WEAK[k]).sum()):5d}"
              f"   other {int(len(R)-vc.get(DOMINANT[k],0)-vc.get(WEAK[k],0)):5d}")

    F = {}
    F["lry_count"] = sum((R[k] == DOMINANT[k]).astype(float) for k in DOMINANT)
    F["lry_minus_weak"] = (sum((R[k] == DOMINANT[k]).astype(float) for k in DOMINANT)
                           - sum((R[k] == WEAK[k]).astype(float) for k in WEAK))
    for k in DOMINANT:
        F[f"dom_{k}"] = (R[k] == DOMINANT[k]).astype(float)

    print(f"\n=== Q2: does it vary inside each private cohort? ===")
    print(f"  {'epitope':<12}{'split':<9}{'n':>5}{'lry_count levels':>18}"
          f"{'mean':>7}{'sd':>7}")
    for e, g in test.groupby("epitope_aa"):
        v = F["lry_count"].reindex(g["tcr_id"])
        print(f"  {e:<12}{g['split'].iloc[0]:<9}{len(g):5d}"
              f"{v.nunique():18d}{v.mean():7.2f}{v.std():7.2f}")

    print(f"\n=== the review's own question: is it just polar_frac again? ===")
    Eg = pd.read_csv(ROOT / "data/interim/interface_energetics.csv"
                     ).set_index("tcr_id")
    pf = Eg["polar_frac"].reindex(test["tcr_id"])
    for nm in ("lry_count", "lry_minus_weak"):
        r = float(np.mean([spearmanr(F[nm].reindex(g["tcr_id"]),
                                     pf.reindex(g["tcr_id"]),
                                     nan_policy="omit").statistic
                           for _, g in test.groupby("epitope_aa")]))
        rd = float(np.mean([spearmanr(F[nm].reindex(g["tcr_id"]),
                                      dr.reindex(g["tcr_id"])).statistic
                            for _, g in test.groupby("epitope_aa")]))
        print(f"  {nm:<16}rho with polar_frac {r:+.3f}   rho with D {rd:+.3f}")

    print(f"\n=== gates; floor 0.5098, record 0.5609, "
          f"pre-registered {PREDICTED_RESID[0]:.2f}-{PREDICTED_RESID[1]:.2f} ===")
    print(f"  {'quantity':<16}{'macro':>8}{'typ rho':>9}{'priv':>6}"
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
        mm = hyp.mean_auc(labels, v)
        r_ = pd.Series(index=v.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in v.index]
            xx = dr.reindex(i).rank(pct=True).to_numpy()
            yy = v.reindex(i).rank(pct=True).to_numpy()
            rr = yy - np.polyval(np.polyfit(xx, yy, 1), xx)
            r_.loc[i] = rr
            leak.append(abs(spearmanr(rr, xx).statistic))
        q = hyp.mean_auc(labels, r_)
        resid, mx = max(q, 1 - q), max(leak)
        print(f"  {nm:<16}{max(mm,1-mm):8.4f}{tg['rho_with_typicality']:+9.3f}"
              f"{'pass' if pv else 'FAIL':>6}{resid:8.4f}{mx:7.3f}"
              f"{'yes' if mx<0.1 else 'NO':>7}")
        rows.append([nm, max(mm, 1 - mm), resid, mx, pv])

    print("\n=== the instrument that decides: stratified inside D-quintiles ===")
    strat = {}
    for nm, s in F.items():
        acc = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in lab.index]
            if len(i) < 40:
                continue
            qq = pd.qcut(dr.reindex(i).rank(method="first"), 4, labels=False)
            for k in range(4):
                sel = [t for t, z in zip(i, qq) if z == k]
                yy = lab.reindex(sel).to_numpy()
                xx = s.reindex(sel).to_numpy()
                if len(set(yy)) < 2 or len(set(xx[~pd.isna(xx)])) < 2:
                    continue
                acc.append(roc_auc_score(yy, pd.Series(xx).fillna(
                    pd.Series(xx).median())))
        if acc:
            above = sum(1 for z in acc if z > 0.5)
            p = binomtest(above, len(acc), 0.5, alternative="greater").pvalue
            strat[nm] = (np.mean(acc), above, len(acc), p)
            print(f"  {nm:<16}{len(acc):2d} strata  mean AUC {np.mean(acc):.4f}"
                  f"  above 0.5: {above}/{len(acc)}  binom p {p:.3f}")

    print("\n=== blended into D, nested LOEO ===")
    base = hyp.mean_auc(labels, dr)
    WS = [0.0, 0.05, 0.10, 0.15]
    print(f"  D alone {base:.4f}   pre-registered nested "
          f"{PREDICTED_NESTED[0]:+.3f} to {PREDICTED_NESTED[1]:+.3f}")
    best = None
    for nm, s in F.items():
        v = s.astype(float).copy()
        for e, g in test.groupby("epitope_aa"):
            x = v.reindex(g["tcr_id"])
            v.loc[g["tcr_id"]] = x.fillna(x.median()).to_numpy()
        v = v.fillna(v.median())
        for sg in (1.0, -1.0):
            sr = hyp.rank_within(sg * v, test)
            per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr) for w in WS}
            vals = {w: float(np.mean([per[w][e] for e in pub])) for w in WS}
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pick = max(WS, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pick][e])
            nest = float(np.mean(held)) - base
            bw = max(WS, key=lambda w: vals[w])
            if sg > 0:
                print(f"  {nm:<16}sign +1 w={bw:.2f}  oracle "
                      f"{vals[bw]-base:+.4f}  nested {nest:+.4f}")
            best = nest if best is None else max(best, nest)

    print(f"\n=== verdict ===")
    rr = [r for r in rows if r[0] == "lry_count"][0]
    print(f"  lry_count residual {rr[2]:.4f} -> "
          f"{'HIT' if PREDICTED_RESID[0] <= rr[2] <= PREDICTED_RESID[1] else 'MISSED'}")
    print(f"  best nested {best:+.4f} -> "
          f"{'HIT' if PREDICTED_NESTED[0] <= best <= PREDICTED_NESTED[1] else 'MISSED'}")

    lines = ["# E149 — TCR surface-expression framework residues (LRY)", "",
             "Generated by `scripts/127_expression_lry.py`.", "",
             "Thomas et al. (PMC6773850): alpha IMGT 96 P→L, beta IMGT 9 S→R,",
             "beta IMGT 10 N→Y raise surface expression 2–6×, persisting in",
             "Jurkat cells lacking endogenous TCR.", "",
             "**Two facts established before scoring.** `full_seq_*` is the V",
             "domain only (103–121 residues, ending at J) — **no constant region**,",
             "so the actual screening construct cannot be checked from competition",
             "data. And holding the V gene fixed, framework variation is **0 of",
             "2,636 rows**, so any framework score is *exactly* a V-gene statistic.",
             "", f"Pre-registered: residual {PREDICTED_RESID[0]:.2f}–"
             f"{PREDICTED_RESID[1]:.2f}, nested {PREDICTED_NESTED[0]:+.3f} to "
             f"{PREDICTED_NESTED[1]:+.3f}.", "",
             "| quantity | macro | residual on D | leak | valid | private var |",
             "|---|---:|---:|---:|---|---|"]
    for nm, mm, resid, mx, pv in rows:
        lines.append(f"| `{nm}` | {mm:.4f} | {resid:.4f} | {mx:.3f} "
                     f"| {'yes' if mx<0.1 else 'NO'} | {'pass' if pv else 'FAIL'} |")
    lines += ["", "| quantity | strata | mean AUC inside D-quintiles | above 0.5 | binom p |",
              "|---|---:|---:|---:|---:|"]
    for nm, (m, ab, n, p) in strat.items():
        lines.append(f"| `{nm}` | {n} | {m:.4f} | {ab}/{n} | {p:.3f} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
