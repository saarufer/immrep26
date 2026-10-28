#!/usr/bin/env python
"""Is polar_frac's germline half just another V-gene statistic?

E109 decomposed `polar_frac` and found its information sits in the germline
loops facing the MHC (residual-on-D 0.5510) rather than in CDR3 facing the
peptide (0.5120, at the 0.5098 floor). That is a clean result, but it comes
with an unresolved objection: **0.5510 is exactly where E103's V-gene docking
conformity sat (0.5458), and that feature turned out to be reducible to the
whole-V-domain alignment term D already carries** — residualising it on the
V-domain term alone dropped it from 0.7272 to 0.5748.

If `pf_germline` is the same thing, then the one structural feature this
project found is, in its useful half, a third restatement of V-gene identity
and nothing more. Four tests, in increasing severity:

  1. **How much is between-gene versus within-gene?** If the quantity is a
     V-gene label, almost all its variance is between genes. Measured by the
     one-way variance ratio over TRAV and TRBV, with cohort held fixed.
  2. **Does it survive replacing each TCR by its gene mean?** The gene-mean
     version is a pure V-gene statistic by construction; the residual is what
     is left of the individual receptor.
  3. **Does it survive residualising on the V-domain term**, the test that
     deflated E103's conformity.
  4. **Does the whole-interface `pf_all` survive all of the above?** E109 found
     the whole beats every part because the parts oppose each other; if that
     survives, the opposition is between individual receptors, not between
     V genes.

No selection step anywhere, so E106's permutation null does not apply.

Run: .venv/bin/python scripts/101_is_polar_a_vgene_statistic.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/78_polar_vgene_test.md"


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])

    def residual_auc(s, ref):
        out = pd.Series(index=s.index, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            ids = [t for t in g["tcr_id"] if t in s.index]
            if len(ids) < 10:
                continue
            x = ref.reindex(ids).rank(pct=True).to_numpy()
            y = s.reindex(ids).rank(pct=True).to_numpy()
            out.loc[ids] = y - np.polyval(np.polyfit(x, y, 1), x)
        m = hyp.mean_auc(labels, out)
        return max(m, 1 - m), out

    E = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    PD = pd.read_csv(ROOT / "data/interim/polar_frac_parts.csv").set_index("tcr_id")
    feats = {}
    for nm, src in (("pf_all", E["polar_frac"]),
                    ("pf_germline", PD["pf_germline"]),
                    ("pf_cdr3", PD["pf_cdr3"]),
                    ("pf_mhc_only", PD["pf_mhc_only"])):
        s = src.reindex(test["tcr_id"]).replace([np.inf, -np.inf], np.nan)
        feats[nm] = s.fillna(s.median())
    prof = pd.read_csv(ROOT / "data/interim/germline_docking.csv").set_index("tcr_id")
    gcols = {"a": [x for x in prof.columns if x.startswith(("a_A1b", "a_A2b"))],
             "b": [x for x in prof.columns if x.startswith(("b_A1b", "b_A2b"))]}
    conf = {}
    for tag, gc in (("a", "TRAV_IMGT"), ("b", "TRBV_IMGT")):
        M = prof[gcols[tag]].to_numpy(float)
        M = M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-9)
        ix = {t: i for i, t in enumerate(prof.index)}
        gs = test.set_index("tcr_id")[gc].reindex(prof.index)
        out = np.full(len(prof), np.nan)
        for gn, sub in gs.groupby(gs):
            ii = np.array([ix[t] for t in sub.index])
            if len(ii) < 5:
                continue
            for k in ii:
                o = ii[ii != k]
                md = M[o].mean(0)
                md /= np.linalg.norm(md) + 1e-9
                out[k] = float(M[k] @ md)
        conf[tag] = pd.Series(out, index=prof.index)
    cf = ((conf["a"] + conf["b"]) / 2)
    feats["E103 conformity"] = cf.reindex(test["tcr_id"]).fillna(cf.median())
    vd = pd.Series(rival.vdomain_topk(test))

    # ---------- 1. between-gene versus within-gene variance ----------
    print("=== 1. how much of each quantity is the V gene? (cohort held fixed) ===")
    print(f"  {'feature':<18}{'eta2 TRAV':>11}{'eta2 TRBV':>11}"
          f"{'eta2 both':>11}{'eta2 epitope':>14}")
    var_rows = []
    for nm, s in feats.items():
        res = []
        for gcol in ("TRAV_IMGT", "TRBV_IMGT", None, "epitope_aa"):
            tot, within = 0.0, 0.0
            for e, g in test.groupby("epitope_aa"):
                v = s.reindex(g["tcr_id"]).to_numpy()
                if gcol is None:
                    key = (g["TRAV_IMGT"].astype(str) + "|"
                           + g["TRBV_IMGT"].astype(str)).to_numpy()
                elif gcol == "epitope_aa":
                    key = None
                else:
                    key = g[gcol].astype(str).to_numpy()
                tot += ((v - v.mean()) ** 2).sum()
                if key is None:
                    continue
                for k in np.unique(key):
                    m = key == k
                    if m.sum() > 1:
                        within += ((v[m] - v[m].mean()) ** 2).sum()
                    # singleton groups contribute zero within-variance
            if gcol == "epitope_aa":
                v = s.reindex(test["tcr_id"]).to_numpy()
                grand = ((v - v.mean()) ** 2).sum()
                res.append(1 - tot / grand)
            else:
                res.append(1 - within / tot)
        print(f"  {nm:<18}{res[0]:11.3f}{res[1]:11.3f}{res[2]:11.3f}{res[3]:14.3f}")
        var_rows.append([nm] + res)

    # ---------- 2-4. the three deflations ----------
    print("\n=== 2-4. what survives each deflation? (residual-on-D) ===")
    print(f"  {'feature':<18}{'raw':>8}{'minus gene mean':>18}"
          f"{'minus V-domain':>16}{'minus both':>12}")
    rows = []
    for nm, s in feats.items():
        raw, _ = residual_auc(s, dr)
        # gene-mean version: replace each TCR by the mean of its (TRAV,TRBV) cell
        gm = pd.Series(index=s.index, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            key = (g["TRAV_IMGT"].astype(str) + "|" + g["TRBV_IMGT"].astype(str))
            v = s.reindex(g["tcr_id"])
            gm.loc[g["tcr_id"]] = v.groupby(key.values).transform("mean").to_numpy()
        res_gm = s - gm
        a_gm, _ = residual_auc(res_gm, dr)
        _, res_vd = residual_auc(s, vd)
        a_vd = residual_auc(res_vd, dr)[0]
        _, r2 = residual_auc(res_gm, vd)
        a_both = residual_auc(r2, dr)[0]
        print(f"  {nm:<18}{raw:8.4f}{a_gm:18.4f}{a_vd:16.4f}{a_both:12.4f}")
        rows.append([nm, raw, a_gm, a_vd, a_both])

    lines = ["# E111 — Is polar_frac's germline half a V-gene statistic?", "",
             "Generated by `scripts/101_is_polar_a_vgene_statistic.py`. No",
             "selection step, so E106's permutation null does not apply.", "",
             "## 1. How much of each quantity is the V gene?", "",
             "Variance explained with the cohort held fixed; the last column is",
             "the between-epitope share for comparison.", "",
             "| feature | η² TRAV | η² TRBV | η² both | η² epitope |",
             "|---|---:|---:|---:|---:|"]
    for r in var_rows:
        lines.append(f"| `{r[0]}` | {r[1]:.3f} | {r[2]:.3f} | {r[3]:.3f} "
                     f"| {r[4]:.3f} |")
    lines += ["", "## 2–4. What survives each deflation", "",
              "All columns are residual-on-D. `minus gene mean` subtracts the",
              "(TRAV, TRBV) cell mean within each cohort, which removes every",
              "pure V-gene statistic by construction. `minus V-domain` is the",
              "test that deflated E103's conformity from 0.7272 to 0.5748.", "",
              "| feature | raw | − gene mean | − V-domain | − both |",
              "|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| `{r[0]}` | **{r[1]:.4f}** | {r[2]:.4f} | {r[3]:.4f} "
                     f"| **{r[4]:.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
