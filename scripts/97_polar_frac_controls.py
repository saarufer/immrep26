#!/usr/bin/env python
"""The E104 controls, as a script instead of a shell heredoc.

Three analyses supporting E103 and E104 were run as inline commands and left no
retained code — the same process failure E60 criticised in a collaborator's
ledger and that E77–E82 repeated here. They are reconstructed as one script so
the numbers in the log can be reproduced:

  1. **Which gravity well does germline docking conformity fall into?** E103
     reported it correlates with everything at +0.37 to +0.46 and loses 60% of
     its margin when residualised on the V-domain term alone.

  2. **Is `polar_frac` CDR3 composition in structural clothing?** E104's central
     scientific claim is that it is not — three sequence-only analogues sit at
     chance. This is the control that makes "the structure earns something the
     sequence cannot" a measurement rather than an assertion.

  3. **Residual-on-D for every component**, E103's diagnostic table.

Unlike the stacking results, none of these involves a selection step, so the
E106 permutation null does not apply to them: they are measured on all 1,782
scorable rows with fixed definitions.

Run: .venv/bin/python scripts/97_polar_frac_controls.py
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

REPORT = ROOT / "reports/74_polar_frac_controls.md"
EPS = 1e-9
POLAR_RES = set("STNQCYHKRDEW")
HYDRO_RES = set("AVLIMFWCY")
# heavy-atom and N+O counts per residue, for the composition-only analogue
NO = {'A': 2, 'R': 5, 'N': 4, 'D': 4, 'C': 2, 'Q': 4, 'E': 4, 'G': 2, 'H': 4,
      'I': 2, 'L': 2, 'K': 3, 'M': 2, 'F': 2, 'P': 2, 'S': 3, 'T': 3, 'W': 3,
      'Y': 3, 'V': 2}
HV = {'A': 5, 'R': 11, 'N': 8, 'D': 8, 'C': 6, 'Q': 9, 'E': 9, 'G': 4, 'H': 10,
      'I': 8, 'L': 8, 'K': 9, 'M': 8, 'F': 11, 'P': 7, 'S': 6, 'T': 7, 'W': 14,
      'Y': 12, 'V': 7}


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

    def wrho(a, b):
        return float(np.mean([spearmanr(a.reindex(g["tcr_id"]),
                                        b.reindex(g["tcr_id"])).statistic
                              for _, g in test.groupby("epitope_aa")]))

    def residual_on(s, ref):
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

    # ---------- the components ----------
    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}
    dens, tn, lf, vv = {}, {}, {}, {}
    fz = np.load(ROOT / "data/interim/contact_fingerprints.npz", allow_pickle=True)
    rfp = fz["fingerprint"].astype(float)
    wh = {t: i for i, t in enumerate(fz["tcr_id"])}
    nzf = rfp / (np.linalg.norm(rfp, axis=1, keepdims=True) + 1e-9)
    for e, g in test.groupby("epitope_aa"):
        idx = np.array([po[t] for t in g["tcr_id"]])
        D = 1.0 - W[idx] @ W[idx].T
        np.fill_diagonal(D, np.inf)
        sd = np.sort(D, axis=1)
        dens.update(dict(zip(g["tcr_id"], -sd[:, :8].mean(axis=1))))
        r1, r2 = np.maximum(sd[:, 0], EPS), np.maximum(sd[:, 1], EPS)
        tn.update(dict(zip(g["tcr_id"],
                           1.0 / np.log(np.maximum(r2 / r1, 1 + EPS)))))
        la, lb = g["cdr3_alpha_aa"].str.len(), g["cdr3_beta_aa"].str.len()
        for t, v in zip(g["tcr_id"], (la.map(la.value_counts(normalize=True))
                                      * lb.map(lb.value_counts(normalize=True))
                                      ).to_numpy()):
            lf[t] = v
        idl = [t for t in g["tcr_id"] if t in wh]
        B = nzf[[wh[t] for t in idl]]
        dd = 1.0 - B @ B.T
        np.fill_diagonal(dd, np.inf)
        for i, t in enumerate(idl):
            vv[t] = -float(np.sort(dd[i])[:8].mean())
    dens, tn, lf = pd.Series(dens), pd.Series(tn), pd.Series(lf)
    fpd = pd.Series(vv)
    E = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    pf = E["polar_frac"].reindex(test["tcr_id"])
    pf = pf.fillna(pf.median())
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
    cf = cf.fillna(cf.median())
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")

    # ---------- 1. which well does conformity fall into? ----------
    print("=== 1. germline conformity: correlated with what? (E103) ===")
    vd = rival.vdomain_topk(test)
    rows1 = []
    for nm, ref in (("V-domain top-k (m50)", vd), ("SCEPTR density", dens),
                    ("pae_tcr_pep_min", boltz["pae_tcr_pep_min"]),
                    ("submitted D", dr)):
        r = wrho(cf, pd.Series(ref))
        print(f"  {nm:<24}{r:+8.3f}   [its own macro "
              f"{hyp.mean_auc(labels, pd.Series(ref)):.4f}]")
        rows1.append([nm, r, hyp.mean_auc(labels, pd.Series(ref))])
    rv, _ = residual_on(cf, pd.Series(vd))
    print(f"  conformity macro {hyp.mean_auc(labels, cf):.4f}; "
          f"residualised on m50 alone -> {rv:.4f}")

    # ---------- 2. is polar_frac CDR3 composition? ----------
    print("\n=== 2. is polar_frac CDR3 composition? (E104's central control) ===")
    def frac(seq, S):
        s = str(seq)
        return sum(ch in S for ch in s) / max(len(s), 1)

    def nof(s):
        s = str(s)
        return sum(NO.get(ch, 2) for ch in s) / max(
            sum(HV.get(ch, 8) for ch in s), 1)

    idx = test.set_index("tcr_id")
    seqs = {
        "CDR3 polar residue fraction": idx.apply(
            lambda r: .5 * frac(r.cdr3_alpha_aa, POLAR_RES)
            + .5 * frac(r.cdr3_beta_aa, POLAR_RES), axis=1),
        "CDR3 hydrophobic fraction": idx.apply(
            lambda r: .5 * frac(r.cdr3_alpha_aa, HYDRO_RES)
            + .5 * frac(r.cdr3_beta_aa, HYDRO_RES), axis=1),
        "CDR3 implied N/O atom fraction": idx.apply(
            lambda r: .5 * nof(r.cdr3_alpha_aa) + .5 * nof(r.cdr3_beta_aa),
            axis=1),
    }
    print(f"  {'analogue':<32}{'rho w/ pf':>11}{'macro':>9}{'resid on D':>12}")
    rows2 = []
    for nm, s in seqs.items():
        r = wrho(pf, s)
        m = hyp.mean_auc(labels, s)
        ra, _ = residual_on(s, dr)
        print(f"  {nm:<32}{r:+11.3f}{m:9.4f}{ra:12.4f}")
        rows2.append([nm, r, m, ra])
    mpf = hyp.mean_auc(labels, pf)
    rpf, _ = residual_on(pf, dr)
    print(f"  {'structural polar_frac':<32}{'--':>11}{mpf:9.4f}{rpf:12.4f}")
    best = max(seqs, key=lambda k: abs(wrho(pf, seqs[k])))
    m2, res = residual_on(pf, seqs[best])
    r2, _ = residual_on(res, dr)
    print(f"\n  polar_frac residualised on '{best}':")
    print(f"    macro {m2:.4f}   resid-on-D {r2:.4f}  "
          f"(raw {max(mpf,1-mpf):.4f} / {rpf:.4f})")

    # ---------- 3. residual-on-D for every component ----------
    print("\n=== 3. residual-on-D, every component (E103's diagnostic) ===")
    from immrep.validation import typicality as typf
    comps = {"structure fingerprint (E83)": fpd, "lenfreq (E96)": lf,
             "TwoNN LID (E100)": -tn, "germline conformity (E103)": cf,
             "polar_frac (E104)": -pf,
             "pae_tcr_pep_min (inside D)": boltz["pae_tcr_pep_min"],
             "typicality reference": typf(test)}
    print(f"  {'feature':<32}{'macro':>8}{'rho D':>9}{'resid on D':>12}")
    rows3 = []
    for nm, s in comps.items():
        s = pd.Series(s).reindex(test["tcr_id"])
        s = s.replace([np.inf, -np.inf], np.nan).fillna(s.median())
        m = hyp.mean_auc(labels, s)
        if m < 0.5:
            s, m = -s, 1 - m
        r = wrho(s, dr)
        ra, _ = residual_on(s, dr)
        print(f"  {nm:<32}{m:8.4f}{r:+9.3f}{ra:12.4f}")
        rows3.append([nm, m, r, ra])

    lines = ["# E107 — The E103/E104 controls, reconstructed as code", "",
             "Generated by `scripts/97_polar_frac_controls.py`. These three",
             "analyses were originally run as inline shell commands with no",
             "retained script. None involves a selection step, so E106's",
             "permutation null does not apply: all are measured on the 1,782",
             "scorable rows with fixed definitions.", "",
             "## 1. What germline conformity is correlated with", "",
             "| reference | ρ with conformity | reference's own macro |",
             "|---|---:|---:|"]
    for r in rows1:
        lines.append(f"| {r[0]} | {r[1]:+.3f} | {r[2]:.4f} |")
    lines += ["", f"Conformity macro {hyp.mean_auc(labels, cf):.4f}; residualised "
              f"on the V-domain term alone it drops to **{rv:.4f}**.", "",
              "## 2. Is polar_frac CDR3 composition?", "",
              "| sequence analogue | ρ with polar_frac | its macro | residual on D |",
              "|---|---:|---:|---:|"]
    for r in rows2:
        lines.append(f"| {r[0]} | {r[1]:+.3f} | {r[2]:.4f} | {r[3]:.4f} |")
    lines += [f"| **structural polar_frac** | — | **{mpf:.4f}** | **{rpf:.4f}** |",
              "", f"Residualised on `{best}`, polar_frac becomes *stronger*: "
              f"macro {m2:.4f}, residual-on-D {r2:.4f}.", "",
              "## 3. Residual-on-D for every component", "",
              "| feature | macro | ρ with D | residual on D |",
              "|---|---:|---:|---:|"]
    for r in rows3:
        lines.append(f"| {r[0]} | {r[1]:.4f} | {r[2]:+.3f} | **{r[3]:.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
