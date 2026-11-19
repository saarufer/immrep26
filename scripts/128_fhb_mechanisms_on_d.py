#!/usr/bin/env python
"""E150 -- FHB's two mechanisms applied to D.

A teammate's `lumos/research-intelligence` branch reaches public macro 0.8748
against D's 0.8600, with private-row rho(D) = 0.892 and a valid residual-on-D
of **0.6216**, the highest this project has measured. Their code exposes two
mechanisms, both of which this project never tried.

**Correction to my own earlier reading.** The +0.0148 gap to D is mostly their
*anchor*, not the mechanism: their `EXPERIMENT_REGISTRY.csv` records anchor
`DH` at **0.8702** and their best variant at 0.8749, so the mechanism itself is
worth **+0.0047** on their anchor, measured in-sample on five epitopes.

**Mechanism 1 -- smooth physicochemical moments.** For each chain, three
properties (Kyte-Doolittle hydropathy /4.5, charge +1 for KRH and -1 for DE,
aromaticity for FWY) are projected onto four cosine modes over the normalised
position, giving 4x3 = 12 numbers per chain. CDR3s are trimmed `[2:-2]` to drop
the conserved C and the terminal F/W. Features are `[a, b, a*b, |a-b|, a*p,
b*p]`. Note `a*p` and `b*p` are **redundant within an epitope** -- the peptide
is constant there, so they are `a` and `b` rescaled by a fixed vector -- so the
effective content is `a, b, a*b, |a-b|` = 48 dimensions.

This connects directly to our own E127/E135: a positional-composition null
recovers **95.6%** of the density's above-chance signal, and mode k=0 is exactly
the compositional mean while k=1..3 are its gradients along the chain. So these
features parametrise the layer where E135 showed the signal lives, through a
smooth low-dimensional lens that D's distance metric does not expose directly.

**Mechanism 2 -- leave-one-source-out median consensus.** Their correction head
is a bounded pairwise model: standardise, build positive-negative pairs *within
each source epitope*, keep only pairs with `|anchor_a - anchor_b| <= 0.15` (the
pairs the anchor cannot already order), and fit `score = anchor + 0.05*tanh(Xw)`
under a pairwise logistic loss with L2 0.01. Then each training epitope is
dropped in turn, one head fitted per deletion, and the **median** of the
resulting deltas is taken -- plus a stricter variant keeping only coordinates
where all heads agree in sign.

This differs from everything we tried. E143 took a median over D's four
*components*; this takes a median over heads fitted with one *training epitope*
deleted, which targets the cross-epitope zero-shot structure the private split
actually has. E147 conditioned on component disagreement; this conditions on
anchor proximity.

**Pre-registered.** Mechanism 1 as a label-free density space: residual-on-D
**0.52-0.57**, nested ~0. Mechanism 2 with D as anchor: nested **-0.002 to
+0.005**, the upper end being what the mechanism bought on their own anchor.

Run: .venv/bin/python scripts/128_fhb_mechanisms_on_d.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/105_fhb_mechanisms_on_d.md"
HYD = dict(zip("ACDEFGHIKLMNPQRSTVWY",
               [1.8, 2.5, -3.5, -3.5, 2.8, -.4, -3.2, 4.5, -3.9, 3.8,
                1.9, -3.5, -1.6, -3.5, -4.5, -.8, -.7, 4.2, -.9, -1.3]))
PAIR_WINDOW = 0.15
CAP = 0.05
L2 = 0.01
PREDICTED_RESID = (0.52, 0.57)
PREDICTED_NESTED = (-0.002, 0.005)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def moments(seq: str) -> np.ndarray:
    v = np.array([[HYD[c] / 4.5,
                   float(c in "KRH") - float(c in "DE"),
                   float(c in "FWY")] for c in seq])
    pos = (np.arange(len(seq)) + .5) / len(seq)
    return np.cos(np.pi * np.arange(4)[:, None] * pos) @ v / len(seq)


def spectral(alpha, beta, peptide):
    f = lambda s: moments(s[2:-2] if len(s) > 4 else s).ravel()  # noqa: E731
    a = np.array([f(s) for s in alpha])
    b = np.array([f(s) for s in beta])
    p = np.array([moments(s).ravel() for s in peptide])
    return np.c_[a, b, a * b, np.abs(a - b), a * p, b * p], a, b


def make_pairs(labels, eps, tr, anchor):
    L, R, W = [], [], []
    for src in np.unique(eps[tr]):
        ix = tr[eps[tr] == src]
        pos, neg = ix[labels[ix] == 1], ix[labels[ix] == 0]
        if len(pos) == 0 or len(neg) == 0:
            continue
        n = len(pos) * len(neg)
        ch = np.linspace(0, n - 1, min(4096, n), dtype=int)
        a, b = pos[ch // len(neg)], neg[ch % len(neg)]
        keep = np.abs(anchor[a] - anchor[b]) <= PAIR_WINDOW
        a, b = a[keep], b[keep]
        if len(a) == 0:
            continue
        L.extend(a); R.extend(b); W.extend(np.full(len(a), 1 / len(a)))
    W = np.array(W, dtype=float)
    return np.array(L), np.array(R), W / W.sum()


def objective(w, x, anchor, L, R, W):
    h = np.tanh(x @ w)
    score = anchor + CAP * h
    diff = 5 * (score[L] - score[R])
    ds = np.zeros(len(anchor))
    coef = -.25 * W * expit(-diff)
    np.add.at(ds, L, coef)
    np.add.at(ds, R, -coef)
    return (float(W @ np.logaddexp(0, -diff) + L2 * w @ w),
            x.T @ (ds * (1 - h * h)) + 2 * L2 * w)


def fit_head(x, anchor, labels, eps, tr, te):
    xx = StandardScaler().fit(x[tr]).transform(x)
    L, R, W = make_pairs(labels, eps, tr, anchor)
    if len(L) == 0:
        return anchor[te].copy()
    opt = minimize(objective, np.zeros(xx.shape[1]),
                   args=(xx, anchor, L, R, W), jac=True,
                   method="L-BFGS-B", options={"maxiter": 150})
    return anchor[te] + CAP * np.tanh(xx[te] @ opt.x)


def main() -> int:
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    test = load_test().reset_index(drop=True)
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    sc = labels.dropna(subset=["pos"]).copy()
    sc["pos"] = sc["pos"].astype(int)
    lab = sc.set_index("tcr_id")["pos"]
    pub = sorted(sc["epitope_aa"].unique())
    dr = hyp.rank_within(
        pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                    ).set_index("tcr_id")["prediction"], test)

    X, A, B = spectral(test["cdr3_alpha_aa"].fillna("CXX").tolist(),
                       test["cdr3_beta_aa"].fillna("CXX").tolist(),
                       test["epitope_aa"].tolist())
    print(f"=== spectral features: {X.shape[1]} dims "
          f"(effective within-epitope: {A.shape[1]*4}) ===")

    # ---------- Mechanism 1: label-free density in spectral space ----------
    print("\n=== mechanism 1: spectral space as a density space (label-free) ===")
    print(f"  pre-registered residual {PREDICTED_RESID[0]:.2f}-{PREDICTED_RESID[1]:.2f}")
    eff = np.c_[A, B, A * B, np.abs(A - B)]
    ids = test["tcr_id"]
    feats = {}
    for nm, M in (("spectral_density", eff), ("ab_coupling_density", A * B)):
        v = pd.Series(index=ids, dtype=float)
        for e, g in test.groupby("epitope_aa"):
            i = g.index.to_numpy()
            Z = M[i]
            Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
            Dm = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
            np.fill_diagonal(Dm, np.inf)
            k = min(8, len(i) - 1)
            v.loc[g["tcr_id"]] = -np.sort(Dm, axis=1)[:, :k].mean(1)
        feats[nm] = v
    print(f"  {'feature':<22}{'macro':>8}{'typ rho':>9}{'priv':>6}"
          f"{'resid':>8}{'leak':>7}{'valid':>7}{'nested':>9}")
    m1 = []
    for nm, v in feats.items():
        tg = typicality_gate(test, sc, v)
        pv, _ = private_variance_gate(test, v)
        mm = hyp.mean_auc(labels, v)
        res = pd.Series(index=v.index, dtype=float)
        leak = []
        for e, g in test.groupby("epitope_aa"):
            i = [t for t in g["tcr_id"] if t in v.index]
            xx = dr.reindex(i).rank(pct=True).to_numpy()
            yy = v.reindex(i).rank(pct=True).to_numpy()
            r = yy - np.polyval(np.polyfit(xx, yy, 1), xx)
            res.loc[i] = r
            leak.append(abs(spearmanr(r, xx).statistic))
        q = hyp.mean_auc(labels, res)
        base = hyp.mean_auc(labels, dr)
        best = None
        for sg in (1.0, -1.0):
            sr = hyp.rank_within(sg * v, test)
            per = {w: hyp.per_epitope(labels, (1 - w) * dr + w * sr)
                   for w in (0.0, .05, .10, .15)}
            held = []
            for e in pub:
                trn = [x for x in pub if x != e]
                pk = max(per, key=lambda w: np.mean([per[w][t] for t in trn]))
                held.append(per[pk][e])
            n = float(np.mean(held)) - base
            best = n if best is None else max(best, n)
        print(f"  {nm:<22}{max(mm,1-mm):8.4f}{tg['rho_with_typicality']:+9.3f}"
              f"{'pass' if pv else 'FAIL':>6}{max(q,1-q):8.4f}{max(leak):7.3f}"
              f"{'yes' if max(leak)<0.1 else 'NO':>7}{best:+9.4f}")
        m1.append([nm, max(mm, 1-mm), max(q, 1-q), max(leak), best])

    # ---------- Mechanism 2: D-anchored correction, median consensus ----------
    print("\n=== mechanism 2: D as anchor, leave-one-source-out median consensus ===")
    print(f"  pre-registered nested {PREDICTED_NESTED[0]:+.3f} to "
          f"{PREDICTED_NESTED[1]:+.3f}  (their mechanism bought +0.0047 on their anchor)")
    eps = test["epitope_aa"].to_numpy()
    anchor = dr.reindex(test["tcr_id"]).to_numpy()
    yv = np.full(len(test), -1)
    m = lab.reindex(test["tcr_id"])
    yv[m.notna().to_numpy()] = m.dropna().astype(int).to_numpy()

    base_per = hyp.per_epitope(labels, dr)
    base = float(np.mean([base_per[e] for e in pub]))
    print(f"  D alone {base:.4f}")
    print(f"  {'variant':<22}" + "".join(f"{e[:3]:>9}" for e in pub) + f"{'nested':>10}")
    rows = []
    for variant in ("single_head", "median", "sign_consensus"):
        gains = {}
        for target in pub:
            te = np.flatnonzero(eps == target)
            srcs = [e for e in pub if e != target]
            lab_m = yv.copy()
            lab_m[np.isin(eps, [target])] = -1          # target labels hidden
            tr_all = np.flatnonzero(lab_m >= 0)
            if variant == "single_head":
                pred = fit_head(X, anchor, lab_m, eps, tr_all, te)
                out = pred
            else:
                heads = []
                for drop in srcs:
                    tt = tr_all[eps[tr_all] != drop]
                    heads.append(fit_head(X, anchor, lab_m, eps, tt, te))
                deltas = np.array([(h - anchor[te]) / CAP for h in heads])
                if variant == "median":
                    out = anchor[te] + CAP * np.median(deltas, axis=0)
                else:
                    s = np.sign(deltas)
                    agree = (np.all(s >= 0, axis=0) | np.all(s <= 0, axis=0))
                    out = anchor[te] + CAP * np.median(deltas * agree, axis=0)
            sub = pd.Series(out, index=test["tcr_id"].to_numpy()[te])
            ids_l = [t for t in sub.index if t in lab.index]
            yy = lab.reindex(ids_l).to_numpy()
            gains[target] = (roc_auc_score(yy, sub.reindex(ids_l))
                             - base_per[target])
        nested = float(np.mean([gains[e] for e in pub]))
        print(f"  {variant:<22}" + "".join(f"{gains[e]:+9.4f}" for e in pub)
              + f"{nested:+10.4f}")
        rows.append([variant, gains, nested])

    print(f"\n=== verdict ===")
    r1 = m1[0]
    print(f"  mechanism 1 residual {r1[2]:.4f} -> "
          f"{'HIT' if PREDICTED_RESID[0] <= r1[2] <= PREDICTED_RESID[1] else 'MISSED'}"
          f"   (leak {r1[3]:.3f} {'valid' if r1[3]<0.1 else 'INVALID'})")
    bn = max(r[2] for r in rows)
    print(f"  mechanism 2 best nested {bn:+.4f} -> "
          f"{'HIT' if PREDICTED_NESTED[0] <= bn <= PREDICTED_NESTED[1] else 'MISSED'}")

    lines = ["# E150 — FHB's two mechanisms applied to D", "",
             "Generated by `scripts/128_fhb_mechanisms_on_d.py`.", "",
             "Their anchor `DH` is 0.8702 and their best variant 0.8749, so the",
             "mechanism is worth **+0.0047** on their own anchor, not the +0.0148",
             "gap to D.", "",
             f"Pre-registered: mechanism 1 residual {PREDICTED_RESID[0]:.2f}–"
             f"{PREDICTED_RESID[1]:.2f}; mechanism 2 nested "
             f"{PREDICTED_NESTED[0]:+.3f} to {PREDICTED_NESTED[1]:+.3f}.", "",
             "## Mechanism 1 — spectral moments as a density space", "",
             "| feature | macro | residual on D | leak | valid | nested |",
             "|---|---:|---:|---:|---|---:|"]
    for nm, mm, rd, lk, nb in m1:
        lines.append(f"| `{nm}` | {mm:.4f} | **{rd:.4f}** | {lk:.3f} "
                     f"| {'yes' if lk<0.1 else 'NO'} | {nb:+.4f} |")
    lines += ["", "## Mechanism 2 — D-anchored correction with source-deletion consensus",
              "", f"D alone **{base:.4f}**.", "",
              "| variant | " + " | ".join(pub) + " | nested |",
              "|---|" + "---:|" * (len(pub) + 1)]
    for v, g, n in rows:
        lines.append(f"| `{v}` | " + " | ".join(f"{g[e]:+.4f}" for e in pub)
                     + f" | **{n:+.4f}** |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
