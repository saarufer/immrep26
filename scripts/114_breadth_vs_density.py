#!/usr/bin/env python
"""E128 -- Does TCR density predict cross-reactivity breadth?

Two opposing readings of the density mechanism have coexisted in this project
without ever being separated:

  high density = typical = a receptor the repertoire makes often, therefore
                 plausibly a **promiscuous** one, broad across peptides
  high density = sits in a validated cluster, therefore a **focused** one,
                 narrow and specific

Only functional panels measuring one receptor against many peptides can tell
them apart, and two exist:

  P-TEAM (Schubert 2024, Zenodo 10.5281/zenodo.11197942)
    human, HLA-A*02:01, **20 TCRs x 172 peptides** (NLVPMVATV + 171 single
    mutants). NLVPMVATV is a competition epitope with 366 cohort rows and
    25,554 paired VDJdb complexes, so the *paired* density -- the measure D
    actually uses -- is computable here.
  P-TEAM mouse APL panels, H-2Kb/SIINFEKL
    educated (16 TCRs) and naive (21 TCRs + OT-I), same 171-variant design.
    SIINFEKL has **12 paired complexes** in VDJdb, so the paired density is
    NOT computable; the mouse panels get a beta-only analogue and are reported
    separately. Mouse H-2Kb and human HLA-A*02:01 are never pooled.

**Power, stated before the result.** 3,440 measurements are not 3,440
independent observations: breadth is one number per receptor, so n=20 for the
human panel. Worse, the 20 are not independent -- TCRs 4-4/82-14 share CDR3B
`CASSLAPGATNEKLFF` with near-identical alpha, 52-10/73-14 share
`CAERLQTGANNLFF`, 77-14/84-3 share `CARNTGNQFYF`. Clonal duplicates are
collapsed and the effective n is reported. A correlation on ~17 points cannot
resolve a modest effect, and every estimate carries a bootstrap interval.

**Leakage discipline.** Every predictor -- density, Pgen, lengths, publicness --
is computed from sequence and from the competition/VDJdb cohorts only. The 172
variant activation values are used *solely* to form the breadth outcome.

Run: .venv/bin/python scripts/114_breadth_vs_density.py
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import load_test  # noqa: E402

PTEAM = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
             "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/pteam/"
             "TcrPrediction_MutatedAPLs-master/data")
VDJDB = ROOT / "data/external/vdjdb/vdjdb-2026-06-03/vdjdb.txt"
REPORT = ROOT / "reports/92_breadth_vs_density.md"
SIGMA = 0.1
THRESHOLDS = (0.3, 0.5, 0.7)
BOOT = 4000
SEED = 0


def clean_gene(g, default_allele="*01"):
    """'TRBV7-6*00(1526.7)' -> 'TRBV7-6*01'; takes the first of a list."""
    if pd.isna(g):
        return None
    s = str(g).split(",")[0].split(" or ")[0].strip()
    s = re.sub(r"\(.*?\)", "", s).strip()
    s = s.split()[0] if s.split() else s
    base = s.split("*")[0]
    return base + default_allele if base else None


def to_frame(rows):
    from libtcrlm.schema.tcr import TravGene, TrbvGene
    f = pd.DataFrame(rows, columns=["TRAV", "CDR3A", "TRAJ",
                                    "TRBV", "CDR3B", "TRBJ"])
    for col, vocab in (("TRAV", {g.name for g in TravGene}),
                       ("TRBV", {g.name for g in TrbvGene})):
        bad = ~f[col].astype(str).str.split("*").str[0].isin(vocab)
        f.loc[bad, col] = None
    return f


def boot_spearman(x, y, rng, n=BOOT):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 5:
        return np.nan, np.nan, np.nan, len(x)
    r = spearmanr(x, y).statistic
    bs = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if len(set(y[i])) > 1 and len(set(x[i])) > 1:
            bs.append(spearmanr(x[i], y[i]).statistic)
    lo, hi = np.percentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan)
    return r, lo, hi, len(x)


def main() -> int:
    import sceptr.variant as variants
    from tcrdist.pgen import OlgaModel

    rng = np.random.default_rng(SEED)
    x = pd.ExcelFile(PTEAM / "Affinity_prediction_cmv.xlsx")
    seqs = x.parse("TCR sequences")
    mean = x.parse("Mean").set_index("Peptide_ID")
    peps = x.parse("peptides")
    seqs["col"] = seqs["TCR id"].str.replace("TCR ", "", regex=False
                                             ).str.replace("-", "_")
    assert set(seqs["col"]) == set(mean.columns), "TCR/column mismatch"
    orig = int(peps.loc[peps["Amino acid"] == "original", "ID"].iloc[0])
    print(f"=== P-TEAM human: {len(seqs)} TCRs x {len(mean)} peptides, "
          f"original = ID {orig} ({peps.loc[peps.ID==orig,'Peptide'].iloc[0]}) ===")

    # ---- outcome: breadth, from the variant panel only ----
    base = mean.loc[orig]
    rel = mean.drop(index=orig).div(base, axis=1)
    print(f"  response to the original peptide: "
          f"{base.min():.1f}-{base.max():.1f} (not a % of max, so each TCR is "
          f"normalised to its own)")
    B = {t: (rel[c] >= t).mean() for t in THRESHOLDS for c in [None] if False}
    breadth = {t: (rel >= t).mean(axis=0) for t in THRESHOLDS}
    for t in THRESHOLDS:
        b = breadth[t]
        print(f"  breadth at >={t:.0%} of own baseline: mean {b.mean():.3f}  "
              f"range {b.min():.3f}-{b.max():.3f}")

    # ---- collapse clonal duplicates ----
    seqs["key"] = seqs["cdr3_a_aa"].astype(str) + "|" + seqs["cdr3_b_aa"].astype(str)
    dup_a = seqs.groupby("cdr3_a_aa").size()
    dup_b = seqs.groupby("cdr3_b_aa").size()
    shared = sorted(set(dup_a[dup_a > 1].index) | set(dup_b[dup_b > 1].index))
    print(f"\n  CDR3s shared by more than one of the 20: {len(shared)} -> {shared}")
    # cluster by shared alpha OR shared beta (single linkage)
    lab = list(range(len(seqs)))
    for col in ("cdr3_a_aa", "cdr3_b_aa"):
        for _, g in seqs.groupby(col):
            if len(g) > 1:
                m = min(lab[i] for i in g.index)
                for i in g.index:
                    lab[i] = m
    changed = True
    while changed:
        changed = False
        for col in ("cdr3_a_aa", "cdr3_b_aa"):
            for _, g in seqs.groupby(col):
                if len(g) > 1:
                    m = min(lab[i] for i in g.index)
                    for i in g.index:
                        if lab[i] != m:
                            lab[i] = m
                            changed = True
    seqs["clone"] = lab
    print(f"  effective n after collapsing clonal duplicates: "
          f"{seqs['clone'].nunique()} of {len(seqs)}")

    # ---- predictors, no variant labels used ----
    model = variants.default()
    frame = to_frame([(clean_gene(r.TRAV), r.cdr3_a_aa, clean_gene(r.TRAJ),
                       clean_gene(r.TRBV), r.cdr3_b_aa, clean_gene(r.TRBJ))
                      for r in seqs.itertuples()])
    Vq = model.calc_vector_representations(frame)
    Vq = Vq / np.linalg.norm(Vq, axis=1, keepdims=True)

    # reference cohort 1: the competition's NLVPMVATV rows
    test = load_test()
    nlv = test[test.epitope_aa == "NLVPMVATV"]
    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    pos = {t: i for i, t in enumerate(c["tcr_id"])}
    Rc = W[[pos[t] for t in nlv["tcr_id"] if t in pos]]

    # reference cohort 2: VDJdb paired NLVPMVATV complexes
    V = pd.read_csv(VDJDB, sep="\t", low_memory=False)
    v = V[(V["antigen.epitope"] == "NLVPMVATV") & (V["complex.id"] != 0)]
    a = v[v.gene == "TRA"].set_index("complex.id")
    b = v[v.gene == "TRB"].set_index("complex.id")
    both = a.index.intersection(b.index)
    vp = pd.DataFrame({"CDR3A": a.loc[both, "cdr3"].to_numpy(),
                       "TRAV": [clean_gene(g) for g in a.loc[both, "v.segm"]],
                       "TRAJ": [clean_gene(g) for g in a.loc[both, "j.segm"]],
                       "CDR3B": b.loc[both, "cdr3"].to_numpy(),
                       "TRBV": [clean_gene(g) for g in b.loc[both, "v.segm"]],
                       "TRBJ": [clean_gene(g) for g in b.loc[both, "j.segm"]]})
    vp = vp.drop_duplicates(subset=["CDR3A", "CDR3B"]).dropna(subset=["CDR3A", "CDR3B"])
    vp = vp.sample(min(4000, len(vp)), random_state=SEED)
    Vv = model.calc_vector_representations(
        to_frame(list(zip(vp.TRAV, vp.CDR3A, vp.TRAJ, vp.TRBV, vp.CDR3B, vp.TRBJ))))
    Vv = Vv / np.linalg.norm(Vv, axis=1, keepdims=True)
    print(f"\n  reference cohorts: competition NLV {len(Rc)}, "
          f"VDJdb paired NLV {len(Vv)} (deduplicated)")

    def density(Q, R):
        d = 1.0 - Q @ R.T
        # drop exact self-matches so a panel TCR present in VDJdb is not its
        # own neighbour
        d = np.where(d < 1e-9, np.nan, d)
        return np.nanmean(np.exp(-d / SIGMA), axis=1)

    P = pd.DataFrame(index=seqs.index)
    P["density_competition"] = density(Vq, Rc)
    P["density_vdjdb"] = density(Vq, Vv)
    k = min(20, len(Vv))
    dv = np.sort(np.where((1.0 - Vq @ Vv.T) < 1e-9, np.nan, 1.0 - Vq @ Vv.T), axis=1)
    P["knn20_vdjdb"] = -np.nanmean(dv[:, :k], axis=1)
    P["len_cdr3a"] = seqs["cdr3_a_aa"].str.len().to_numpy()
    P["len_cdr3b"] = seqs["cdr3_b_aa"].str.len().to_numpy()
    P["len_total"] = P["len_cdr3a"] + P["len_cdr3b"]
    ma = OlgaModel(chain_folder="human_T_alpha", recomb_type="VJ")
    mb = OlgaModel(chain_folder="human_T_beta", recomb_type="VDJ")
    with np.errstate(divide="ignore"):
        P["log_pgen_a"] = np.log10(np.array(
            ma.compute_aa_cdr3_pgens(seqs["cdr3_a_aa"].tolist()), float) + 1e-30)
        P["log_pgen_b"] = np.log10(np.array(
            mb.compute_aa_cdr3_pgens(seqs["cdr3_b_aa"].tolist()), float) + 1e-30)
    P["log_pgen_paired"] = P["log_pgen_a"] + P["log_pgen_b"]
    vb = set(V[V.gene == "TRB"]["cdr3"])
    va = set(V[V.gene == "TRA"]["cdr3"])
    cnt = V.groupby(["gene", "cdr3"]).size()
    P["publicness"] = [
        int(cnt.get(("TRA", r.cdr3_a_aa), 0)) + int(cnt.get(("TRB", r.cdr3_b_aa), 0))
        for r in seqs.itertuples()]

    print("\n=== breadth vs each predictor (human, HLA-A*02:01) ===")
    print("  positive rho = higher predictor goes with BROADER cross-reactivity")
    rows = []
    for t in THRESHOLDS:
        y = breadth[t].reindex(seqs["col"]).to_numpy()
        print(f"\n  --- breadth threshold >={t:.0%} "
              f"(mean breadth {np.mean(y):.3f}) ---")
        print(f"  {'predictor':<22}{'rho':>8}{'95% CI':>18}{'n':>4}"
              f"{'rho (clones)':>14}")
        for col in P.columns:
            r, lo, hi, n = boot_spearman(P[col], y, rng)
            # collapsed to one point per clone
            g = pd.DataFrame({"c": seqs["clone"], "x": P[col].to_numpy(), "y": y})
            gg = g.groupby("c").mean()
            rc, _, _, nc = boot_spearman(gg.x, gg.y, rng)
            print(f"  {col:<22}{r:+8.3f}  [{lo:+.3f},{hi:+.3f}]{n:5d}"
                  f"{rc:+14.3f}")
            rows.append([t, col, r, lo, hi, n, rc, nc])

    print("\n=== the two hypotheses, judged ===")
    main_rows = [r for r in rows if r[0] == 0.5 and "density" in r[1]]
    for _, col, r, lo, hi, n, rc, nc in main_rows:
        verdict = ("broader" if lo > 0 else
                   "narrower" if hi < 0 else "UNRESOLVED (CI spans 0)")
        print(f"  {col:<22}rho {r:+.3f} [{lo:+.3f},{hi:+.3f}]  -> {verdict}")

    lines = ["# E128 — Does density predict cross-reactivity breadth?", "",
             "Generated by `scripts/114_breadth_vs_density.py`.", "",
             f"P-TEAM human panel: {len(seqs)} TCRs × {len(mean)} peptides "
             f"(NLVPMVATV + 171 single mutants), HLA-A\\*02:01.",
             f"Breadth is one number per receptor, so **n={len(seqs)}**, and after "
             f"collapsing clonal duplicates **effective n="
             f"{seqs['clone'].nunique()}**. 3,440 measurements are not 3,440 "
             "independent observations.", "",
             "| threshold | predictor | ρ | 95% CI | ρ (clone-collapsed) |",
             "|---|---|---:|---|---:|"]
    for t, col, r, lo, hi, n, rc, nc in rows:
        lines.append(f"| ≥{t:.0%} | `{col}` | {r:+.3f} "
                     f"| [{lo:+.3f}, {hi:+.3f}] | {rc:+.3f} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
