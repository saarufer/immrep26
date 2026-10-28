#!/usr/bin/env python
"""Where in the interface does polar_frac's information live?

E104 found that the N/O share of a TCR's interface atoms carries ~0.059 of AUC
that D does not have, and E107 confirmed it is not reducible to CDR3
composition — three sequence-only analogues sit at chance. So the fold decides
something the sequence does not, and the question is **which part of the
interface decides it**.

Two decompositions, each answering a specific earlier claim:

  **By partner.** TCR interface atoms split into those facing the peptide and
  those facing the MHC. E36 concluded that the peptide-specific component of
  any structural score is at chance — an anchor-matched decoy margin of +0.0061
  against decoy-choice noise of 0.036. If polar_frac's signal sits in the
  peptide-facing half, that is in tension with E36 and worth knowing; if it
  sits in the MHC-facing half, it is about docking rather than recognition and
  E36 stands.

  **By loop.** CDR3 is somatically generated, CDR1 and CDR2 are germline. E103
  built a V-gene docking conformity feature and found it reducible to the
  V-domain term D already holds. If polar_frac's signal is in CDR3, it is about
  this particular receptor; if in the germline loops, it is another V-gene
  statistic and should behave like E103's.

CDR3 residues are located by matching the annotated CDR3 string against the
chain's own residue sequence, which fails for a handful of structures and is
reported rather than imputed. Every quantity is scored by residual-on-D, the
E103 diagnostic, since raw AUC has repeatedly been the wrong screen here.

No selection step anywhere, so E106's permutation null does not apply.

Run: .venv/bin/python scripts/99_polar_frac_decomposition.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/76_polar_frac_decomposition.md"
CACHE = ROOT / "data/interim/polar_frac_parts.csv"
PRED = ROOT / "data/interim/boltz/out5/boltz_results_inputs/predictions"
NEAR = 9.0          # same interface definition as E104
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def parse(path: Path):
    xyz, el, ch, ri, rn = [], [], [], [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            e = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if e == "H":
                continue
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            el.append(e)
            ch.append(line[21])
            ri.append(int(line[22:26]))
            rn.append(line[17:20].strip())
    return (np.asarray(xyz, float), np.asarray(el), np.asarray(ch),
            np.asarray(ri), np.asarray(rn))


def chain_seq(ch, ri, rn, which):
    m = ch == which
    order = {}
    for i, r in zip(np.where(m)[0], ri[m]):
        order.setdefault(r, rn[i])
    keys = sorted(order)
    return "".join(AA3.get(order[k], "X") for k in keys), keys


def features(path: Path, cdr3a: str, cdr3b: str):
    xyz, el, ch, ri, rn = parse(path)
    if len(xyz) == 0:
        return None
    pol = np.isin(el, ["N", "O"])
    tcr = np.isin(ch, ["A", "B"])
    pep = ch == "P"
    mhc = ch == "M"
    if not tcr.any() or not pep.any() or not mhc.any():
        return None
    tp, tm = cKDTree(xyz[pep]), cKDTree(xyz[mhc])
    ti = np.where(tcr)[0]
    d_pep = np.array([len(x) for x in tp.query_ball_point(xyz[ti], r=NEAR)])
    d_mhc = np.array([len(x) for x in tm.query_ball_point(xyz[ti], r=NEAR)])
    iface = (d_pep > 0) | (d_mhc > 0)
    out = {}
    pol_i = pol[ti]

    def f(mask):
        return float(pol_i[mask].mean()) if mask.sum() else np.nan

    out["pf_all"] = f(iface)
    out["pf_pep_facing"] = f(d_pep > 0)
    out["pf_mhc_only"] = f((d_mhc > 0) & (d_pep == 0))
    out["n_pep_facing"] = int((d_pep > 0).sum())

    # CDR3 versus germline loops, located by sequence match
    cdr3_mask = np.zeros(len(ti), dtype=bool)
    found = 0
    for which, cdr3 in (("A", cdr3a), ("B", cdr3b)):
        seq, keys = chain_seq(ch, ri, rn, which)
        pos = seq.find(str(cdr3))
        if pos < 0:
            continue
        found += 1
        want = set(keys[pos:pos + len(str(cdr3))])
        cdr3_mask |= (ch[ti] == which) & np.isin(ri[ti], list(want))
    out["cdr3_located"] = found
    out["pf_cdr3"] = f(iface & cdr3_mask)
    out["pf_germline"] = f(iface & ~cdr3_mask)
    return out


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    if CACHE.exists():
        print(f"loading {CACHE.name}")
        F = pd.read_csv(CACHE).set_index("tcr_id")
    else:
        rows = []
        sub = test[["tcr_id", "cdr3_alpha_aa", "cdr3_beta_aa"]].values
        for i, (t, ca, cb) in enumerate(sub):
            p = PRED / t / f"{t}_model_0.pdb"
            if not p.exists():
                continue
            r = features(p, ca, cb)
            if r:
                rows.append({"tcr_id": t, **r})
            if (i + 1) % 400 == 0:
                print(f"  {i+1}/{len(sub)}")
        F = pd.DataFrame(rows).set_index("tcr_id")
        F.to_csv(CACHE)
        print(f"wrote {CACHE.name} ({len(F)})")

    print(f"\n{len(F)} structures")
    print(f"  CDR3 located on both chains: "
          f"{(F['cdr3_located'] == 2).mean():.1%}, one chain "
          f"{(F['cdr3_located'] == 1).mean():.1%}, neither "
          f"{(F['cdr3_located'] == 0).mean():.1%}")
    print(f"  peptide-facing interface atoms: median "
          f"{F['n_pep_facing'].median():.0f} of ~326 interface atoms")

    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])

    def residual_auc(s):
        res = pd.Series(index=s.index, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            ids = [t for t in g["tcr_id"] if t in s.index]
            if len(ids) < 10:
                continue
            x = dr.reindex(ids).rank(pct=True).to_numpy()
            y = s.reindex(ids).rank(pct=True).to_numpy()
            res.loc[ids] = y - np.polyval(np.polyfit(x, y, 1), x)
        m = hyp.mean_auc(labels, res)
        return max(m, 1 - m)

    print("\n=== the decomposition, scored by residual-on-D ===")
    print("  E104's whole-interface polar_frac reaches 0.5586; a term already")
    print("  inside D reaches 0.5098.")
    print(f"  {'part':<28}{'macro':>8}{'inv':>8}{'rho D':>9}"
          f"{'resid on D':>12}{'coverage':>10}")
    parts = ["pf_all", "pf_pep_facing", "pf_mhc_only", "pf_cdr3", "pf_germline"]
    rows = []
    for nm in parts:
        s = F[nm].reindex(test["tcr_id"]).replace([np.inf, -np.inf], np.nan)
        cov = float(s.notna().mean())
        s = s.fillna(s.median())
        m = hyp.mean_auc(labels, s)
        rd = float(np.mean([spearmanr(s.reindex(g["tcr_id"]),
                                      dr.reindex(g["tcr_id"])).statistic
                            for _, g in test.groupby("epitope_aa")]))
        ra = residual_auc(s)
        print(f"  {nm:<28}{m:8.4f}{1-m:8.4f}{rd:+9.3f}{ra:12.4f}{cov:10.1%}")
        rows.append([nm, m, rd, ra, cov])

    print("\n=== are the parts independent of each other? ===")
    print(f"  {'':<16}" + "".join(f"{p.replace('pf_',''):>14}" for p in parts))
    cormat = []
    for a in parts:
        sa = F[a].reindex(test["tcr_id"])
        sa = sa.fillna(sa.median())
        row = []
        for b in parts:
            sb = F[b].reindex(test["tcr_id"])
            sb = sb.fillna(sb.median())
            row.append(float(np.mean([
                spearmanr(sa.reindex(g["tcr_id"]),
                          sb.reindex(g["tcr_id"])).statistic
                for _, g in test.groupby("epitope_aa")])))
        print(f"  {a.replace('pf_',''):<16}" + "".join(f"{v:+14.3f}" for v in row))
        cormat.append(row)

    # does the peptide-facing part survive removing the MHC-facing part?
    print("\n=== does either half survive removing the other? ===")
    surv = []
    for a, b in (("pf_pep_facing", "pf_mhc_only"), ("pf_mhc_only", "pf_pep_facing"),
                 ("pf_cdr3", "pf_germline"), ("pf_germline", "pf_cdr3")):
        sa = F[a].reindex(test["tcr_id"])
        sa = sa.fillna(sa.median())
        sb = F[b].reindex(test["tcr_id"])
        sb = sb.fillna(sb.median())
        res = pd.Series(index=sa.index, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            ids = list(g["tcr_id"])
            x = sb.reindex(ids).rank(pct=True).to_numpy()
            y = sa.reindex(ids).rank(pct=True).to_numpy()
            res.loc[ids] = y - np.polyval(np.polyfit(x, y, 1), x)
        ra = residual_auc(res)
        print(f"  {a:<16} after removing {b:<16} resid-on-D {ra:.4f}"
              f"   (alone {residual_auc(sa):.4f})")
        surv.append([a, b, ra, residual_auc(sa)])

    lines = ["# E109 — Where in the interface polar_frac's information lives", "",
             "Generated by `scripts/99_polar_frac_decomposition.py`. No selection",
             "step, so E106's permutation null does not apply.", "",
             f"{len(F)} structures; CDR3 located on both chains in "
             f"{(F['cdr3_located'] == 2).mean():.1%}; median "
             f"{F['n_pep_facing'].median():.0f} peptide-facing atoms.", "",
             "## Scored by residual-on-D", "",
             "E104's whole-interface value is 0.5586; a term already inside D is",
             "0.5098.", "",
             "| part | macro | inverted | ρ with D | residual on D | coverage |",
             "|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| `{r[0]}` | {r[1]:.4f} | {1-r[1]:.4f} | {r[2]:+.3f} "
                     f"| **{r[3]:.4f}** | {r[4]:.1%} |")
    lines += ["", "## Mutual correlation of the parts", "",
              "| | " + " | ".join(p.replace("pf_", "") for p in parts) + " |",
              "|---|" + "---:|" * len(parts)]
    for a, row in zip(parts, cormat):
        lines.append(f"| {a.replace('pf_','')} | "
                     + " | ".join(f"{v:+.3f}" for v in row) + " |")
    lines += ["", "## Does either half survive removing the other?", "",
              "| part | after removing | residual on D | alone |",
              "|---|---|---:|---:|"]
    for r in surv:
        lines.append(f"| `{r[0]}` | `{r[1]}` | **{r[2]:.4f}** | {r[3]:.4f} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
