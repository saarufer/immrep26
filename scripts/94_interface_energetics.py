#!/usr/bin/env python
"""Interface energetics: physics rather than geometry, gated before blending.

E83 and E103 read the TCR-pMHC interface as *geometry* — which peptide
positions are touched, which MHC helix each chain sits on. Neither asks whether
the interface is physically well made. Buried surface area, hydrogen bonds,
salt bridges, charge complementarity and packing density are the quantities a
structural biologist would use to judge whether a modelled complex is a real
one, and they are not contact counts: two interfaces can touch the same
residues and bury very different area, or make the same contacts with opposite
electrostatics.

**The procedural change E103 earned.** Six candidates in a row cleared the
`typicality_gate` and delivered nothing, and the reason was visible in one
number: residualising each on D within each epitope left 0.535–0.557, against
0.5098 for a term already inside D. So this run computes the residual-on-D AUC
**first**, before any blending, and only a feature that clears the 0.557 set by
E83's fingerprint — the best so far — is worth carrying further.

SASA is computed by Shrake-Rupley restricted to atoms within 9 A of the
partner, since no atom further than twice the probe-inflated radius can change
burial; that makes the differential exact while touching a few hundred atoms
instead of several thousand.

Run: .venv/bin/python scripts/94_interface_energetics.py
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
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/71_interface_energetics.md"
CACHE = ROOT / "data/interim/interface_energetics.csv"
PRED = ROOT / "data/interim/boltz/out5/boltz_results_inputs/predictions"
PROBE = 1.4
NPTS = 92
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
HB_MAX = 3.5        # heavy-atom donor-acceptor distance; the PDBs carry no H
SB_MAX = 4.0
POS_AT = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ"),
          ("HIS", "ND1"), ("HIS", "NE2")}
NEG_AT = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sphere_points(n: int) -> np.ndarray:
    """Golden-spiral points on the unit sphere."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


SPH = sphere_points(NPTS)


def sasa(xyz: np.ndarray, rad: np.ndarray) -> float:
    """Shrake-Rupley total solvent-accessible area."""
    if len(xyz) == 0:
        return 0.0
    rp = rad + PROBE
    tree = cKDTree(xyz)
    nbr = tree.query_ball_point(xyz, r=(rp.max() + rp.max()))
    total = 0.0
    for i, nb in enumerate(nbr):
        nb = [j for j in nb if j != i]
        pts = xyz[i] + rp[i] * SPH
        if nb:
            d = np.linalg.norm(pts[:, None, :] - xyz[nb][None, :, :], axis=2)
            free = (d >= rp[nb][None, :]).all(axis=1)
        else:
            free = np.ones(NPTS, dtype=bool)
        total += 4 * np.pi * rp[i] ** 2 * free.mean()
    return float(total)


def parse(path: Path):
    rec = {"xyz": [], "el": [], "ch": [], "rn": [], "an": [], "ri": []}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if el == "H":
                continue
            rec["xyz"].append((float(line[30:38]), float(line[38:46]),
                               float(line[46:54])))
            rec["el"].append(el)
            rec["ch"].append(line[21])
            rec["rn"].append(line[17:20].strip())
            rec["an"].append(line[12:16].strip())
            rec["ri"].append(int(line[22:26]))
    return {k: np.asarray(v) for k, v in rec.items()}


def features(a) -> dict:
    xyz = a["xyz"].astype(float)
    rad = np.array([VDW.get(e, 1.7) for e in a["el"]])
    tcr = np.isin(a["ch"], ["A", "B"])
    pmhc = np.isin(a["ch"], ["M", "P"])
    if tcr.sum() == 0 or pmhc.sum() == 0:
        return None
    # only atoms that can change burial need SASA
    t1 = cKDTree(xyz[pmhc])
    t2 = cKDTree(xyz[tcr])
    near_t = np.zeros(len(xyz), dtype=bool)
    near_p = np.zeros(len(xyz), dtype=bool)
    near_t[np.where(tcr)[0]] = np.array(
        [len(x) > 0 for x in t1.query_ball_point(xyz[tcr], r=9.0)])
    near_p[np.where(pmhc)[0]] = np.array(
        [len(x) > 0 for x in t2.query_ball_point(xyz[pmhc], r=9.0)])
    sel = near_t | near_p
    s_cplx = sasa(xyz[sel], rad[sel])
    s_t = sasa(xyz[sel & tcr], rad[sel & tcr])
    s_p = sasa(xyz[sel & pmhc], rad[sel & pmhc])
    bsa = s_t + s_p - s_cplx

    # polar contacts across the interface
    pol = np.isin(a["el"], ["N", "O"])
    ti = np.where(tcr & pol)[0]
    pi = np.where(pmhc & pol)[0]
    hb = 0
    if len(ti) and len(pi):
        d = np.linalg.norm(xyz[ti][:, None, :] - xyz[pi][None, :, :], axis=2)
        hb = int((d <= HB_MAX).sum())
    key = list(zip(a["rn"], a["an"]))
    pos = np.array([k in POS_AT for k in key])
    neg = np.array([k in NEG_AT for k in key])
    sb = 0
    for pa, nb_ in ((tcr & pos, pmhc & neg), (tcr & neg, pmhc & pos)):
        ia, ib = np.where(pa)[0], np.where(nb_)[0]
        if len(ia) and len(ib):
            d = np.linalg.norm(xyz[ia][:, None, :] - xyz[ib][None, :, :], axis=2)
            sb += int((d <= SB_MAX).sum())

    # interface composition and packing
    t_if = np.where(tcr & near_t)[0]
    p_if = np.where(pmhc & near_p)[0]
    q_t = int(pos[t_if].sum()) - int(neg[t_if].sum())
    q_p = int(pos[p_if].sum()) - int(neg[p_if].sum())
    tree_p = cKDTree(xyz[pmhc])
    pack = float(np.mean([len(x) for x in
                          tree_p.query_ball_point(xyz[t_if], r=6.0)])) \
        if len(t_if) else 0.0
    return {"bsa": bsa, "hbonds": hb, "saltbridges": sb,
            "charge_product": -q_t * q_p, "polar_frac": float(pol[t_if].mean())
            if len(t_if) else 0.0, "packing": pack,
            "n_if_tcr": len(t_if), "n_if_pmhc": len(p_if),
            "bsa_per_atom": bsa / max(len(t_if), 1)}


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())

    if CACHE.exists():
        print(f"loading {CACHE.name}")
        F = pd.read_csv(CACHE).set_index("tcr_id")
    else:
        rows, ids = [], test["tcr_id"].tolist()
        for i, t in enumerate(ids):
            p = PRED / t / f"{t}_model_0.pdb"
            if not p.exists():
                continue
            f = features(parse(p))
            if f is None:
                continue
            rows.append({"tcr_id": t, **f})
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(ids)}")
        F = pd.DataFrame(rows).set_index("tcr_id")
        F.to_csv(CACHE)
        print(f"wrote {CACHE.name} ({len(F)})")

    print(f"\n{len(F)} structures")
    print(F.describe().loc[["mean", "std", "min", "max"]].T.to_string())

    # ---------- the gate, run FIRST ----------
    d_sub = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                        ).set_index("tcr_id")["prediction"]
    dr = hyp.rank_within(d_sub, test)
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))

    def residual_auc(s: pd.Series) -> float:
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

    print("\n=== residual-on-D gate, run BEFORE any blending ===")
    print("  benchmark: E83's fingerprint 0.5570 is the best any candidate has")
    print("  reached; a term already inside D reaches 0.5098.")
    print(f"  {'feature':<18}{'macro':>8}{'inv':>8}{'rho D':>9}"
          f"{'resid on D':>12}{'typ rho':>10}{'gates':>8}")
    rows = []
    for nm in F.columns:
        s = F[nm].replace([np.inf, -np.inf], np.nan)
        s = s.fillna(s.median())
        m = hyp.mean_auc(labels, s)
        rd = float(np.mean([spearmanr(s.reindex(g["tcr_id"]),
                                      dr.reindex(g["tcr_id"])).statistic
                            for _, g in test.groupby("epitope_aa")]))
        ra = residual_auc(s)
        g_ = typicality_gate(test, labels.dropna(subset=["pos"]), s)
        pv, _ = private_variance_gate(test, s)
        ok = "pass" if (abs(g_["rho_with_typicality"]) < 0.6378 and pv) else "FAIL"
        print(f"  {nm:<18}{m:8.4f}{1-m:8.4f}{rd:+9.3f}{ra:12.4f}"
              f"{g_['rho_with_typicality']:+10.3f}{ok:>8}")
        rows.append([nm, m, rd, ra, g_["rho_with_typicality"], ok])

    best = max(rows, key=lambda r: r[3])
    print(f"\n  best residual: {best[0]} at {best[3]:.4f}")
    carry = [r for r in rows if r[3] > 0.5570]
    print(f"  clearing E83's 0.5570: {len(carry)}"
          + (f" -> {[r[0] for r in carry]}" if carry else " -- none"))

    # ---------- only blend what cleared ----------
    blend = []
    if carry:
        print(f"\n=== blending the {len(carry)} that cleared ===")
        for r in carry:
            s = F[r[0]].replace([np.inf, -np.inf], np.nan)
            s = s.fillna(s.median())
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
                pk = max(per_w, key=lambda w: np.mean([per_w[w][t] for t in tr]))
                held.append(per_w[pk][e])
            nested = float(np.mean(held))
            print(f"  {r[0]:<18} oracle w={bw:.2f} {vals[bw]:.4f} "
                  f"({vals[bw]-base:+.4f})  nested {nested:.4f} ({nested-base:+.4f})")
            blend.append([r[0], bw, vals[bw], vals[bw] - base, nested, nested - base])
    else:
        print("\n  nothing cleared the gate, so nothing is blended -- which is")
        print("  the point of running the gate first.")

    lines = ["# E104 — Interface energetics, gated before blending", "",
             "Generated by `scripts/94_interface_energetics.py` from the Boltz PDBs.",
             "Buried surface area by Shrake-Rupley restricted to atoms within 9 A",
             "of the partner, hydrogen bonds as N/O pairs within 3.5 A (the PDBs",
             "carry no hydrogens), salt bridges as charged groups within 4 A.", "",
             "## The residual-on-D gate, run first", "",
             "E103's benchmark: the best candidate so far (E83's fingerprint)",
             "reaches **0.5570** after D is removed; a term already inside D",
             "reaches 0.5098.", "",
             "| feature | macro | inverted | ρ with D | **residual on D** | ρ typicality | gates |",
             "|---|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {1-r[1]:.4f} | {r[2]:+.3f} "
                     f"| **{r[3]:.4f}** | {r[4]:+.3f} | {r[5]} |")
    lines += ["", f"Best residual: `{best[0]}` at **{best[3]:.4f}**. "
              f"Clearing 0.5570: **{len(carry)}**.", ""]
    if blend:
        lines += ["## Blended into D", "", f"D alone {base:.4f}.", "",
                  "| feature | oracle w | oracle | vs D | nested | vs D |",
                  "|---|---:|---:|---:|---:|---:|"]
        for b in blend:
            lines.append(f"| {b[0]} | {b[1]:.2f} | {b[2]:.4f} | {b[3]:+.4f} "
                         f"| {b[4]:.4f} | **{b[5]:+.4f}** |")
        lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
