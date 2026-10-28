#!/usr/bin/env python
"""Direction (b): the interface read from the pMHC side.

Every structural measurement in this project reads the TCR — its confidence,
its contacts, its interface polarity. The pMHC side has never been looked at,
and it is not symmetric with the TCR side in an interesting way: **within one
epitope the peptide and the MHC are identical molecules**, so any variation in
a pMHC-side quantity comes purely from *which part of an unchanging surface a
particular receptor covers*. It is a readout of docking geometry from the
opposite bank.

Computed for each complex, over pMHC atoms within 9 A of the TCR:

  polarity of the covered pMHC surface, split into peptide and MHC
  how much of the peptide is covered at all
  which peptide positions are covered (centroid and spread along the groove)
  the polarity *difference* across the interface, TCR side minus pMHC side

The last is the one with a reason behind it: a well-made interface pairs polar
with polar and apolar with apolar, so the mismatch is a complementarity
measure, and it is the one quantity here that is not a property of either
partner alone.

Scored by residual-on-D throughout, per E103, since raw AUC has repeatedly
been the wrong screen. No selection step, so E106's null does not apply.

Run: .venv/bin/python scripts/102_pmhc_side_polarity.py
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

REPORT = ROOT / "reports/79_pmhc_side_polarity.md"
CACHE = ROOT / "data/interim/pmhc_side.csv"
PRED = ROOT / "data/interim/boltz/out5/boltz_results_inputs/predictions"
NEAR = 9.0


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def parse(path: Path):
    xyz, el, ch, ri = [], [], [], []
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
    return (np.asarray(xyz, float), np.asarray(el), np.asarray(ch),
            np.asarray(ri))


def features(path: Path):
    xyz, el, ch, ri = parse(path)
    if len(xyz) == 0:
        return None
    pol = np.isin(el, ["N", "O"])
    tcr = np.isin(ch, ["A", "B"])
    pep, mhc = ch == "P", ch == "M"
    if not (tcr.any() and pep.any() and mhc.any()):
        return None
    tt = cKDTree(xyz[tcr])
    out = {}
    cov = {}
    for tag, m in (("pep", pep), ("mhc", mhc)):
        idx = np.where(m)[0]
        near = np.array([len(x) > 0 for x in
                         tt.query_ball_point(xyz[idx], r=NEAR)])
        cov[tag] = (idx, near)
        out[f"pmhc_pol_{tag}"] = (float(pol[idx][near].mean())
                                  if near.any() else np.nan)
        out[f"pmhc_n_{tag}"] = int(near.sum())
    # how far along the groove the covered peptide sits, and how spread
    pidx, pnear = cov["pep"]
    if pnear.any():
        rr = ri[pidx][pnear].astype(float)
        out["pep_cover_centroid"] = float(rr.mean())
        out["pep_cover_spread"] = float(rr.std())
        out["pep_frac_covered"] = float(len(np.unique(rr))
                                        / max(len(np.unique(ri[pidx])), 1))
    else:
        out["pep_cover_centroid"] = out["pep_cover_spread"] = np.nan
        out["pep_frac_covered"] = 0.0
    # complementarity: TCR-side polarity minus pMHC-side polarity
    tidx = np.where(tcr)[0]
    tp = cKDTree(xyz[pep | mhc])
    tnear = np.array([len(x) > 0 for x in
                      tp.query_ball_point(xyz[tidx], r=NEAR)])
    tcr_pol = float(pol[tidx][tnear].mean()) if tnear.any() else np.nan
    both = np.concatenate([cov["pep"][0][cov["pep"][1]],
                           cov["mhc"][0][cov["mhc"][1]]])
    pmhc_pol = float(pol[both].mean()) if len(both) else np.nan
    out["polarity_mismatch"] = tcr_pol - pmhc_pol
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
        ids = test["tcr_id"].tolist()
        for i, t in enumerate(ids):
            p = PRED / t / f"{t}_model_0.pdb"
            if not p.exists():
                continue
            r = features(p)
            if r:
                rows.append({"tcr_id": t, **r})
            if (i + 1) % 400 == 0:
                print(f"  {i+1}/{len(ids)}")
        F = pd.DataFrame(rows).set_index("tcr_id")
        F.to_csv(CACHE)
        print(f"wrote {CACHE.name} ({len(F)})")

    print(f"\n{len(F)} structures")
    print("\n=== within-epitope variance: is there anything to read? ===")
    print(f"  {'epitope':<12}{'split':<9}" + "".join(
        f"{c.replace('pmhc_','').replace('_',''):>14}"
        for c in ("pmhc_pol_pep", "pmhc_pol_mhc", "polarity_mismatch",
                  "pep_frac_covered")))
    for e, g in test.groupby("epitope_aa"):
        row = [F[c].reindex(g["tcr_id"]).std() for c in
               ("pmhc_pol_pep", "pmhc_pol_mhc", "polarity_mismatch",
                "pep_frac_covered")]
        print(f"  {e:<12}{g['split'].mode().iloc[0]:<9}"
              + "".join(f"{v:14.4f}" for v in row))

    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))

    def residual_auc(s):
        """Residual-on-D with a validity check.

        E103 introduced this diagnostic without one, and E112 found the gap:
        for a feature with few distinct values, percentile-rank gives every
        tied row the same rank, so subtracting a linear fit on D leaves a
        residual that is monotone in D *within each tie group* and smuggles D
        back in. `pep_frac_covered` has 7 values with 91.5% at 1.0 and reached
        |rho(residual, D)| = 0.806, producing a spurious 0.6626. Every other
        feature ever gated this way audits clean at <= 0.041, so no earlier
        conclusion moves -- but the check is now returned, not assumed.
        """
        out = pd.Series(index=s.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            ids = [t for t in g["tcr_id"] if t in s.index]
            if len(ids) < 10:
                continue
            x = dr.reindex(ids).rank(pct=True).to_numpy()
            y = s.reindex(ids).rank(pct=True).to_numpy()
            r = y - np.polyval(np.polyfit(x, y, 1), x)
            out.loc[ids] = r
            leak.append(abs(spearmanr(r, x).statistic))
        m = hyp.mean_auc(labels, out)
        return max(m, 1 - m), float(np.mean(leak))

    print("\n=== residual-on-D gate (floor 0.5098; best so far 0.5586) ===")
    print(f"  {'feature':<22}{'macro':>8}{'inv':>8}{'rho D':>9}"
          f"{'resid on D':>12}{'leak':>9}")
    rows = []
    for nm in F.columns:
        s = F[nm].reindex(test["tcr_id"]).replace([np.inf, -np.inf], np.nan)
        if s.nunique() < 5:
            continue
        s = s.fillna(s.median())
        m = hyp.mean_auc(labels, s)
        rd = float(np.mean([spearmanr(s.reindex(g["tcr_id"]),
                                      dr.reindex(g["tcr_id"])).statistic
                            for _, g in test.groupby("epitope_aa")]))
        ra, leak = residual_auc(s)
        flag = "INVALID" if leak > 0.30 else ("suspect" if leak > 0.15 else "")
        print(f"  {nm:<22}{m:8.4f}{1-m:8.4f}{rd:+9.3f}{ra:12.4f}"
              f"{leak:9.3f} {flag}")
        rows.append([nm, m, rd, ra, leak, flag])

    carry = [r for r in rows if r[3] > 0.5586 and not r[5]]
    print(f"\n  clearing the best-so-far 0.5586: {len(carry)}"
          + (f" -> {[r[0] for r in carry]}" if carry else " -- none"))
    blend = []
    for r in (carry if carry else sorted(rows, key=lambda x: -x[3])[:2]):
        s = F[r[0]].reindex(test["tcr_id"]).replace([np.inf, -np.inf], np.nan)
        s = s.fillna(s.median())
        sv = s if hyp.mean_auc(labels, s) >= 0.5 else -s
        sr = rank(sv)
        per_w = {0.0: hyp.per_epitope(labels, dr)}
        for w in (0.05, 0.10, 0.15, 0.20):
            per_w[w] = hyp.per_epitope(labels, (1 - w) * dr + w * sr)
        vals = {w: float(np.mean([per_w[w][e] for e in pub])) for w in per_w}
        bw = max([w for w in vals if w > 0], key=lambda w: vals[w])
        held = []
        for e in pub:
            tr = [x for x in pub if x != e]
            pk = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
            held.append(per_w[pk][e])
        n = float(np.mean(held))
        tag = "cleared" if carry else "best two anyway"
        print(f"  [{tag}] {r[0]:<22} oracle w={bw:.2f} {vals[bw]:.4f} "
              f"({vals[bw]-base:+.4f})  nested {n:.4f} ({n-base:+.4f})")
        blend.append([r[0], bw, vals[bw], vals[bw] - base, n, n - base])

    lines = ["# E112 — The interface read from the pMHC side", "",
             "Generated by `scripts/102_pmhc_side_polarity.py`. Within one",
             "epitope the peptide and MHC are identical molecules, so any",
             "variation here is purely *which part of an unchanging surface a",
             "particular receptor covers*.", "",
             "## Residual-on-D gate", "",
             "Floor 0.5098; the best any candidate has reached is 0.5586.", "",
             "| feature | macro | inverted | ρ with D | residual on D |",
             "|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| `{r[0]}` | {r[1]:.4f} | {1-r[1]:.4f} | {r[2]:+.3f} "
                     f"| **{r[3]:.4f}** |")
    if blend:
        lines += ["", "## Blended into D", "", f"D alone {base:.4f}.", "",
                  "| feature | oracle w | oracle | vs D | nested | vs D |",
                  "|---|---:|---:|---:|---:|---:|"]
        for b in blend:
            lines.append(f"| `{b[0]}` | {b[1]:.2f} | {b[2]:.4f} | {b[3]:+.4f} "
                         f"| {b[4]:.4f} | **{b[5]:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
