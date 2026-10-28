#!/usr/bin/env python
"""Direction (a): a V-gene prior on revalidation, learned and transferred.

E111 established that `polar_frac`'s information sits at the resolution of the
(TRAV, TRBV) pair — the within-gene residual carries nothing — and that the
V-domain *sequence* alignment D already holds does not capture it. The obvious
next question is whether the structural detour is needed at all: if the signal
is a V-gene property, estimate it **directly from labels** and see whether that
is cleaner.

This is the first supervised feature the project has built, so it needs
discipline the unsupervised ones did not:

  **Epitope confounding.** TRBV19 validates at 0.76 pooled, but TRBV19 is
  GILGFVFTL's public-motif gene and GILGFVFTL's base rate is 0.66. Every gene
  effect is therefore estimated after removing each epitope's own mean.

  **Transfer.** A gene effect is only usable if it holds on an epitope it was
  not fitted on. Fitted on four, applied to the fifth, every time.

  **Shrinkage.** A gene seen three times has a meaningless rate. Effects are
  shrunk toward zero by n/(n+k) with k chosen inside the training folds.

  **Coverage.** Single genes are seen in public for ~100% of private TCRs, but
  **(TRAV, TRBV) pairs for only 41–60%**, so the pair version is reported with
  its coverage and a fallback.

The comparison that decides it is against `polar_frac`, which reaches the same
gene-level information without using a single label.

Run: .venv/bin/python scripts/103_vgene_prior.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/80_vgene_prior.md"
SHRINK = [2, 5, 10, 20, 50]


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
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    lab = labels.dropna(subset=["pos"]).copy()
    lab["pos"] = lab["pos"].astype(int)
    ti = test.set_index("tcr_id")
    y = lab.set_index("tcr_id")["pos"]
    dr = rank(pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
                          ).set_index("tcr_id")["prediction"])
    base = float(np.mean([hyp.per_epitope(labels, dr)[e] for e in pub]))

    def keyser(cols):
        if len(cols) == 1:
            return ti[cols[0]].astype(str)
        return ti[cols[0]].astype(str) + "|" + ti[cols[1]].astype(str)

    def fit_effect(train_eps, cols, k):
        """Gene effect = mean of (label - its epitope's mean) over train epitopes."""
        sub = lab[lab["epitope_aa"].isin(train_eps)].copy()
        sub["centred"] = sub["pos"] - sub.groupby("epitope_aa")["pos"].transform("mean")
        key = keyser(cols).reindex(sub["tcr_id"]).to_numpy()
        df = pd.DataFrame({"k": key, "c": sub["centred"].to_numpy()})
        g = df.groupby("k")["c"].agg(["sum", "size"])
        return (g["sum"] / (g["size"] + k)).to_dict(), set(g.index)

    def apply_effect(eff, cols, ids):
        key = keyser(cols).reindex(ids)
        return pd.Series([eff.get(k, 0.0) for k in key], index=ids), \
            pd.Series([k in eff for k in key], index=ids)

    print(f"D {base:.4f}\n")
    print("=== transfer: fit on four epitopes, score the fifth ===")
    print(f"  {'genes used':<22}{'shrink k':>9}" + "".join(f"{e[:3]:>8}" for e in pub)
          + f"{'macro':>9}{'coverage':>10}")
    results = {}
    for cols, nm in ((["TRBV_IMGT"], "TRBV only"),
                     (["TRAV_IMGT"], "TRAV only"),
                     (["TRAV_IMGT", "TRBV_IMGT"], "TRAV|TRBV pair")):
        for k in SHRINK:
            per, cov = {}, []
            for e in pub:
                tr = [x for x in pub if x != e]
                eff, seen = fit_effect(tr, cols, k)
                ids = test[test["epitope_aa"] == e]["tcr_id"]
                v, ok = apply_effect(eff, cols, ids)
                cov.append(float(ok.mean()))
                yy = y.reindex(ids)
                m = yy.notna().to_numpy()
                if m.sum() < 10 or yy[m].nunique() < 2:
                    continue
                per[e] = roc_auc_score(yy[m], v.to_numpy()[m])
            mac = float(np.mean(list(per.values())))
            results[(nm, k)] = (per, mac, float(np.mean(cov)))
            print(f"  {nm:<22}{k:>9}" + "".join(f"{per.get(e, np.nan):8.3f}" for e in pub)
                  + f"{mac:9.4f}{np.mean(cov):10.1%}")
    best = max(results, key=lambda kk: results[kk][1])
    print(f"\n  best transferred: {best[0]}, shrink k={best[1]}, "
          f"macro {results[best][1]:.4f}")

    # ---------- the honest version: everything nested ----------
    print("\n=== nested: the gene set and k also chosen inside each fold ===")
    held, picks = [], []
    for e in pub:
        tr = [x for x in pub if x != e]
        best_in, bv = None, -1.0
        for cols, nm in ((["TRBV_IMGT"], "TRBV"), (["TRAV_IMGT"], "TRAV"),
                         (["TRAV_IMGT", "TRBV_IMGT"], "pair")):
            for k in SHRINK:
                inner = []
                for h in tr:
                    tr2 = [x for x in tr if x != h]
                    eff, _ = fit_effect(tr2, cols, k)
                    ids = test[test["epitope_aa"] == h]["tcr_id"]
                    v, _ = apply_effect(eff, cols, ids)
                    yy = y.reindex(ids)
                    m = yy.notna().to_numpy()
                    if m.sum() < 10 or yy[m].nunique() < 2:
                        continue
                        inner.append(roc_auc_score(yy[m], v.to_numpy()[m]))
                if inner and np.mean(inner) > bv:
                    bv, best_in = float(np.mean(inner)), (cols, nm, k)
        cols, nm, k = best_in
        eff, _ = fit_effect(tr, cols, k)
        ids = test[test["epitope_aa"] == e]["tcr_id"]
        v, _ = apply_effect(eff, cols, ids)
        yy = y.reindex(ids)
        m = yy.notna().to_numpy()
        held.append(roc_auc_score(yy[m], v.to_numpy()[m]))
        picks.append(f"{nm}/k={k}")
        print(f"  {e:<12} chose {nm:<5} k={k:<3} -> {held[-1]:.4f}")
    nested_alone = float(np.mean(held))
    print(f"  nested macro of the prior alone: {nested_alone:.4f}")

    # ---------- blended into D, with the prior rebuilt per fold ----------
    print("\n=== blended into D, the prior rebuilt inside every fold ===")
    E = pd.read_csv(ROOT / "data/interim/interface_energetics.csv").set_index("tcr_id")
    pf = E["polar_frac"].reindex(test["tcr_id"])
    pfr = rank(-pf.fillna(pf.median()))
    cols, k = ["TRAV_IMGT", "TRBV_IMGT"], best[1]
    WS = [0.0, 0.05, 0.10, 0.15, 0.20]
    rows = []
    for tag, use_cols in (("TRBV only", ["TRBV_IMGT"]),
                          ("TRAV|TRBV pair", ["TRAV_IMGT", "TRBV_IMGT"])):
        held2 = []
        for e in pub:
            tr = [x for x in pub if x != e]
            eff, _ = fit_effect(tr, use_cols, k)
            v, _ = apply_effect(eff, use_cols, test["tcr_id"])
            vr = rank(v)
            inner = {}
            for w in WS:
                s = (1 - w) * dr + w * vr
                per = hyp.per_epitope(labels, s)
                inner[w] = float(np.mean([per[t] for t in tr]))
            bw = max(inner, key=inner.get)
            per = hyp.per_epitope(labels, (1 - bw) * dr + bw * vr)
            held2.append(per[e])
        n = float(np.mean(held2))
        print(f"  {tag:<22} nested {n:.4f} ({n-base:+.4f} vs D)")
        rows.append([tag, n, n - base])
    # reference: polar_frac, label-free, same blend protocol
    held3 = []
    for e in pub:
        tr = [x for x in pub if x != e]
        inner = {}
        for w in WS:
            per = hyp.per_epitope(labels, (1 - w) * dr + w * pfr)
            inner[w] = float(np.mean([per[t] for t in tr]))
        bw = max(inner, key=inner.get)
        per = hyp.per_epitope(labels, (1 - bw) * dr + bw * pfr)
        held3.append(per[e])
    npf = float(np.mean(held3))
    print(f"  {'polar_frac (label-free)':<22} nested {npf:.4f} ({npf-base:+.4f} vs D)")

    rho = float(np.mean([spearmanr(
        apply_effect(fit_effect(pub, cols, k)[0], cols, g["tcr_id"])[0],
        pf.reindex(g["tcr_id"])).statistic
        for _, g in test.groupby("epitope_aa")]))
    print(f"\n  within-epitope rho(V-gene prior, polar_frac) = {rho:+.3f}")

    lines = ["# E113 — A V-gene prior on revalidation, learned and transferred", "",
             "Generated by `scripts/103_vgene_prior.py`. The project's first",
             "supervised feature; every gene effect is estimated after removing",
             "each epitope's own mean, fitted on four epitopes and applied to the",
             "fifth, shrunk by n/(n+k).", "",
             "## Transfer, with the gene set and k fixed", "",
             "| genes | k | " + " | ".join(e[:3] for e in pub) + " | macro | coverage |",
             "|---|---:|" + "---:|" * (len(pub) + 2)]
    for (nm, k_), (per, mac, cov) in results.items():
        lines.append(f"| {nm} | {k_} | "
                     + " | ".join(f"{per.get(e, float('nan')):.3f}" for e in pub)
                     + f" | **{mac:.4f}** | {cov:.0%} |")
    lines += ["", f"Fully nested, gene set and k chosen inside each fold: "
              f"**{nested_alone:.4f}** (picks: {', '.join(picks)}).", "",
              "## Blended into D, the prior rebuilt inside every fold", "",
              f"D alone {base:.4f}.", "",
              "| variant | nested | vs D |", "|---|---:|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.4f} | **{r[2]:+.4f}** |")
    lines += [f"| polar_frac (label-free reference) | {npf:.4f} | **{npf-base:+.4f}** |",
              "", f"Within-epitope ρ between the V-gene prior and `polar_frac`: "
              f"**{rho:+.3f}**.", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
