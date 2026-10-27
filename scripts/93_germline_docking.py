#!/usr/bin/env python
"""Germline-encoded docking: a structural criterion with a literature basis.

Everything structural tried here so far has been a confidence scalar from the
folding model (E14, E19, E28, E36) or a contact count treated as coordinates
(E83). None of it asks the question structural immunology actually asks of a
TCR-pMHC complex: **does this receptor dock the way its V genes are built to
dock?**

The literature basis is the germline-encoded interaction codon
(Marrack & Kappler; Garcia and colleagues; Feng et al. 2007). TCR V-gene
segments carry evolutionarily conserved CDR1/CDR2 contacts with the MHC
helices, and class I complexes show a canonical diagonal polarity: the
**alpha chain over the MHC alpha-2 helix, the beta chain over alpha-1**. The
relationship is V-gene specific — different TRBV segments favour different
contact positions on the alpha-1 helix — and departures from it are rare and
noteworthy rather than routine.

That gives two features of a kind this project has not built:

  **polarity**   does chain A sit on alpha-2 and chain B on alpha-1, or not
  **conformity** does this TCR's MHC-contact profile match the *modal profile
                 of other TCRs carrying the same V gene*

Conformity is the sharper one and is the structural analogue of everything the
sequence side has been doing: instead of "is this TCR typical of its cohort",
it asks "is this complex docked the way this V gene docks". It is label-free,
it varies within every cohort, and its reference set is grouped by V gene
rather than by epitope, so it is not another epitope-cohort statistic.

Two cautions carried from earlier results. E14 established that Boltz
confidence mostly reports how canonical a TCR is, and E36 that the
peptide-specific component of any structural score is at chance — so a
docking-geometry feature should be expected to behave like a TCR property, and
that is fine for task A, which asks whether an annotation is real, not which
peptide it binds. And E75 found the contact footprint nearly identical across
epitopes, so the variance this depends on must be *within* a cohort; that is
measured first and the experiment stops if it is not there.

Run: .venv/bin/python scripts/93_germline_docking.py
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

REPORT = ROOT / "reports/70_germline_docking.md"
CACHE = ROOT / "data/interim/germline_docking.csv"
PRED = ROOT / "data/interim/boltz/out5/boltz_results_inputs/predictions"
CONTACT = 5.0          # heavy-atom contact cutoff, the usual 5 A
A1 = range(57, 85)     # MHC class I alpha-1 helix, standard numbering
A2 = range(138, 177)   # alpha-2 helix
NBIN = 12              # contact profile resolution along each helix


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def parse(path: Path):
    """Heavy-atom coordinates and residue numbers per chain, straight from the
    PDB text -- Biopython's parser is 30x slower and nothing here needs it."""
    xyz = {"A": [], "B": [], "M": [], "P": []}
    res = {"A": [], "B": [], "M": [], "P": []}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[76:78].strip() == "H":
                continue
            ch = line[21]
            if ch not in xyz:
                continue
            xyz[ch].append((float(line[30:38]), float(line[38:46]),
                            float(line[46:54])))
            res[ch].append(int(line[22:26]))
    return ({k: np.asarray(v, dtype=float) for k, v in xyz.items()},
            {k: np.asarray(v, dtype=int) for k, v in res.items()})


def contact_profile(xyz, res, chain: str) -> dict:
    """Per-MHC-residue contact counts for one TCR chain, plus helix sums."""
    if len(xyz[chain]) == 0 or len(xyz["M"]) == 0:
        return None
    tree = cKDTree(xyz["M"])
    pairs = tree.query_ball_point(xyz[chain], r=CONTACT)
    hit = np.zeros(200, dtype=float)
    for lst in pairs:
        for j in lst:
            r = res["M"][j]
            if 0 < r < 200:
                hit[r] += 1
    return hit


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())

    if CACHE.exists():
        print(f"loading cached contact profiles from {CACHE.name}")
        prof = pd.read_csv(CACHE).set_index("tcr_id")
    else:
        rows = []
        ids = test["tcr_id"].tolist()
        for i, t in enumerate(ids):
            p = PRED / t / f"{t}_model_0.pdb"
            if not p.exists():
                continue
            xyz, res = parse(p)
            ha = contact_profile(xyz, res, "A")
            hb = contact_profile(xyz, res, "B")
            if ha is None or hb is None:
                continue
            rec = {"tcr_id": t}
            for tag, h in (("a", ha), ("b", hb)):
                rec[f"{tag}_a1"] = float(h[list(A1)].sum())
                rec[f"{tag}_a2"] = float(h[list(A2)].sum())
                rec[f"{tag}_tot"] = float(h.sum())
                for k, chunk in enumerate(np.array_split(np.array(list(A1)), NBIN)):
                    rec[f"{tag}_A1b{k}"] = float(h[chunk].sum())
                for k, chunk in enumerate(np.array_split(np.array(list(A2)), NBIN)):
                    rec[f"{tag}_A2b{k}"] = float(h[chunk].sum())
            rows.append(rec)
            if (i + 1) % 400 == 0:
                print(f"  parsed {i+1}/{len(ids)}")
        prof = pd.DataFrame(rows).set_index("tcr_id")
        prof.to_csv(CACHE)
        print(f"wrote {CACHE.name}  ({len(prof)} structures)")

    test = test[test["tcr_id"].isin(prof.index)].copy()
    print(f"\n{len(prof)} structures with contact profiles")

    # ---------- 0. is there any within-cohort variance to work with? ----------
    pol_num = prof["a_a2"] + prof["b_a1"]
    pol_den = prof["a_a1"] + prof["a_a2"] + prof["b_a1"] + prof["b_a2"] + 1.0
    polarity = pol_num / pol_den
    print("\n=== 0. canonical polarity: alpha on alpha-2, beta on alpha-1 ===")
    print(f"  {'epitope':<12}{'split':<9}{'mean':>8}{'sd':>8}{'min':>8}{'max':>8}"
          f"{'% canonical':>13}")
    sp_of = test.drop_duplicates("epitope_aa").set_index("epitope_aa")["split"]
    var_rows = []
    for e, g in test.groupby("epitope_aa"):
        v = polarity.reindex(g["tcr_id"]).to_numpy()
        print(f"  {e:<12}{sp_of[e]:<9}{v.mean():8.3f}{v.std():8.3f}"
              f"{v.min():8.3f}{v.max():8.3f}{(v > 0.5).mean()*100:13.1f}")
        var_rows.append([e, sp_of[e], v.mean(), v.std(), v.min(), v.max(),
                         float((v > 0.5).mean())])

    # ---------- 1. V-gene conformity ----------
    print("\n=== 1. V-gene docking conformity ===")
    print("  For each V gene, the modal contact profile over all TCRs carrying")
    print("  it; each TCR scored by cosine similarity to its own V gene's mode.")
    cols = {"a": [c for c in prof.columns if c.startswith(("a_A1b", "a_A2b"))],
            "b": [c for c in prof.columns if c.startswith(("b_A1b", "b_A2b"))]}
    conf = {}
    gene_n = {}
    for tag, gcol in (("a", "TRAV_IMGT"), ("b", "TRBV_IMGT")):
        M = prof[cols[tag]].to_numpy(float)
        M = M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-9)
        idx = {t: i for i, t in enumerate(prof.index)}
        genes = test.set_index("tcr_id")[gcol].reindex(prof.index)
        out = np.full(len(prof), np.nan)
        for gname, sub in genes.groupby(genes):
            ii = np.array([idx[t] for t in sub.index])
            if len(ii) < 5:
                continue
            gene_n[(tag, gname)] = len(ii)
            for k in ii:                       # leave-one-out mode
                others = ii[ii != k]
                mode = M[others].mean(axis=0)
                mode /= np.linalg.norm(mode) + 1e-9
                out[k] = float(M[k] @ mode)
        conf[tag] = pd.Series(out, index=prof.index)
        ok = np.isfinite(out)
        print(f"  {gcol:<12} genes with >=5 members: "
              f"{sum(1 for (t, _) in gene_n if t == tag):3d}   "
              f"covered {ok.mean():.1%}   conformity mean {np.nanmean(out):.3f} "
              f"sd {np.nanstd(out):.3f}")

    feats = {
        "canonical polarity": polarity,
        "TRAV docking conformity": conf["a"],
        "TRBV docking conformity": conf["b"],
        "conformity (both chains)": (conf["a"] + conf["b"]) / 2,
        "alpha contacts on MHC": prof["a_tot"],
        "beta contacts on MHC": prof["b_tot"],
    }

    # ---------- 2. score them ----------
    from immrep.validation import private_variance_gate, typicality_gate
    print("\n=== 2. on task A ===")
    print(f"  {'feature':<28}{'macro':>8}{'inv':>8}{'typ rho':>10}{'gates':>8}"
          f"{'coverage':>10}")
    rows = []
    for nm, s in feats.items():
        s = pd.Series(s).replace([np.inf, -np.inf], np.nan)
        cov = float(s.notna().mean())
        s = s.fillna(s.median())
        feats[nm] = s
        g = typicality_gate(test, labels.dropna(subset=["pos"]), s)
        pv, _ = private_variance_gate(test, s)
        m = hyp.mean_auc(labels, s)
        rho = g["rho_with_typicality"]
        ok = "pass" if (abs(rho) < 0.6378 and pv) else "FAIL"
        print(f"  {nm:<28}{m:8.4f}{1-m:8.4f}{rho:+10.3f}{ok:>8}{cov:10.1%}")
        rows.append([nm, m, rho, ok, cov])

    # ---------- 3. decorrelation and blending ----------
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))
    print(f"\n=== 3. blended into D ({base:.4f}) ===")
    blend = []
    for nm, s in feats.items():
        rd = float(np.mean([spearmanr(s.reindex(g["tcr_id"]),
                                      dr.reindex(g["tcr_id"])).statistic
                            for _, g in test.groupby("epitope_aa")]))
        signed = s if hyp.mean_auc(labels, s) >= 0.5 else -s
        sr = hyp.rank_within(signed, test)
        per_w = {0.0: hyp.per_epitope(labels, dr)}
        for w in (0.05, 0.10, 0.15, 0.20):
            per_w[w] = hyp.per_epitope(labels, (1 - w) * dr + w * sr)
        vals = {w: float(np.mean([per_w[w][e] for e in pub])) for w in per_w}
        bw = max([w for w in vals if w > 0], key=lambda w: vals[w])
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            pick = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
            held.append(per_w[pick][e])
        nested = float(np.mean(held))
        print(f"  {nm:<28} rho(D) {rd:+.3f}   oracle w={bw:.2f} {vals[bw]:.4f} "
              f"({vals[bw]-base:+.4f})   nested {nested:.4f} ({nested-base:+.4f})")
        blend.append([nm, rd, bw, vals[bw], vals[bw] - base, nested, nested - base])

    lines = ["# E103 — Germline-encoded docking", "",
             "Generated by `scripts/93_germline_docking.py` from the Boltz PDBs.",
             "", "## 0. Canonical polarity, and is there within-cohort variance?",
             "", "`(alpha on alpha-2 + beta on alpha-1) / all TCR-MHC helix contacts`",
             "", "| epitope | split | mean | sd | min | max | % canonical |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for r in var_rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:.3f} | {r[3]:.3f} | {r[4]:.3f} "
                     f"| {r[5]:.3f} | {r[6]*100:.1f} |")
    lines += ["", "## 2. On task A", "",
              "| feature | macro | inverted | ρ typicality | gates | coverage |",
              "|---|---:|---:|---:|---|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {1-r[1]:.4f} | {r[2]:+.3f} "
                     f"| {r[3]} | {r[4]:.1%} |")
    lines += ["", "## 3. Blended into D", "", f"D alone **{base:.4f}**.", "",
              "| feature | ρ with D | oracle w | oracle | vs D | nested | vs D |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for b in blend:
        lines.append(f"| {b[0]} | {b[1]:+.3f} | {b[2]:.2f} | {b[3]:.4f} "
                     f"| {b[4]:+.4f} | {b[5]:.4f} | **{b[6]:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
