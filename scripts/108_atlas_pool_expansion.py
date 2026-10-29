#!/usr/bin/env python
"""Does the TCR Motif Atlas add reference receptors where the private split needs them?

A collaborator proposed re-opening reference-pool expansion with a new source.
E51 measured pool expansion as one of the few routes that ever worked
(+0.0112 public), and E84 measured why it does not transfer: the pool term's
contribution shows a dose-response that **vanishes below roughly 60
same-epitope pool rows**, and the private epitopes hold **11, 191 and 5**.

So the only question worth asking of a new source is not "does it have new
records" but **"does it have new records for CINGVCWTV and LTDEMIAQY"**, the
two sitting in the dead zone. The TCR Motif Atlas (44,431 records, 178 studies,
13 databases, with `study` and `database` retained) is the first source that
can be checked that way.

After restricting to human paired-chain records and de-duplicating against
both the competition rows and D's existing 86,221-row pool — exactly and again
after normalising away the leading C and trailing F/W, which agree:

| epitope | split | Atlas | already held | **new** | source |
|---|---|---:|---:|---:|---|
| **LTDEMIAQY** | **private** | 262 | 134 | **128** | PMID:40628266 |
| LLWNGPMAV | private | 384 | 372 | 12 | PMID:32081129 |
| **CINGVCWTV** | **private** | 239 | 239 | **0** | — |
| GILGFVFTL | public | 2475 | 2374 | 101 | 6 studies |
| YLQPRTFLL | public | 409 | 342 | 67 | 4 studies |
| NLVPMVATV | public | 261 | 223 | 38 | 4 studies |
| GLCTLVAML | public | 384 | 366 | 18 | 2 studies |
| TTDPSFLGRY | public | 409 | 409 | 0 | — |

**LTDEMIAQY goes from 5 same-epitope pool rows to 133.** That is the one cell
where E84's curve says the pool term switches from dead to live, and it is a
private epitope.

**Provenance audited, per E92's lesson.** All 128 come from PMID:40628266 —
Drost et al., *Cell Genomics* 2025, "Benchmarking of T cell receptor-epitope
predictors with ePytope-TCR" — which **predates** the IMMREP26 revalidation
(PMID:42822443, 2026). These are reference receptors from a benchmarking
aggregation, not validated private-split rows, so no competition label enters
the pool.

**What can and cannot be measured.** LTDEMIAQY is private: the gain there is
unmeasurable and can only be extrapolated from E84's dose-response. What *is*
measurable is whether adding the Atlas records to the **public** pools changes
anything, which tests the mechanism rather than the bet. Both are reported,
and the extrapolation is labelled as such.

Run: .venv/bin/python scripts/108_atlas_pool_expansion.py
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

REPORT = ROOT / "reports/85_atlas_pool_expansion.md"
ATLAS = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
             "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/full.csv")
TOPK = 8


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
    pool_mod = load_module("pool", "scripts/33_pool_negative_control.py")
    pub = sorted(labels.dropna(subset=["pos"])["epitope_aa"].unique())
    eps = sorted(test["epitope_aa"].unique())
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731

    import sceptr
    sceptr.disable_hardware_acceleration()

    A = pd.read_csv(ATLAS, low_memory=False)
    A = A[(A["species"] == "HomoSapiens") & A["cdr3_TRA"].notna()
          & A["cdr3_TRB"].notna()]
    pool = pool_mod.load_pool()
    print(f"Atlas {len(A)} human paired rows; existing pool {len(pool)} rows")

    # ---------- build the query embeddings and both pools ----------
    query, extra_old, extra_new = {}, {}, {}
    for e, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        frame = pd.DataFrame({
            "cdr3a": cohort["cdr3_alpha_aa"], "cdr3b": cohort["cdr3_beta_aa"],
            "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
            "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"], "ep": e})
        query[e] = (cohort, pool_mod.embed(frame))
        keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
        same = pool[pool["ep"] == e]
        old = same[[(a, b) not in keys for a, b in
                    zip(same["cdr3a"], same["cdr3b"])]].drop_duplicates(
                        ["cdr3a", "cdr3b"])
        extra_old[e] = old
        held = keys | set(zip(old["cdr3a"], old["cdr3b"]))
        a = A[A["peptide"] == e]
        add = a[[(x, y) not in held for x, y in
                 zip(a["cdr3_TRA"], a["cdr3_TRB"])]].drop_duplicates(
                     ["cdr3_TRA", "cdr3_TRB"])
        if len(add):
            add = pd.DataFrame({"cdr3a": add["cdr3_TRA"], "cdr3b": add["cdr3_TRB"],
                                "trav": add["TRAV"], "traj": add["TRAJ"],
                                "trbv": add["TRBV"], "trbj": add["TRBJ"],
                                "ep": e})
        extra_new[e] = pd.concat([old, add]) if len(add) else old
        print(f"  {e:<12} cohort {len(cohort):4d}  pool {len(old):5d}"
              f" -> {len(extra_new[e]):5d}  (+{len(extra_new[e])-len(old)})")

    def pooled_topk(extra_by_ep) -> pd.Series:
        out = pd.Series(index=test["tcr_id"], dtype=float)
        for e, (cohort, vectors) in query.items():
            ex = extra_by_ep.get(e)
            ref = vectors if ex is None or ex.empty else np.vstack(
                [vectors, pool_mod.embed(ex)])
            d = 1.0 - vectors @ ref.T
            d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
            out.loc[cohort["tcr_id"]] = -np.sort(d, axis=1)[:, :TOPK].mean(axis=1)
        return out

    print("\nembedding the expanded pools ...")
    m30_old, m30_new = rank(pooled_topk(extra_old)), rank(pooled_topk(extra_new))

    # ---------- rebuild D with each ----------
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure = rank(boltz["pae_tcr_pep_min"])
    frozen = decomp.submitted_score(test)
    density3 = (frozen - 0.3 * structure) / 0.7
    vdomain = rank(rival.vdomain_topk(test))
    d_old = (density3 + m30_old + vdomain + structure) / 4
    d_new = (density3 + m30_new + vdomain + structure) / 4

    po, pn = hyp.per_epitope(labels, d_old), hyp.per_epitope(labels, d_new)
    print("\n=== effect on the PUBLIC epitopes (the measurable part) ===")
    print(f"  {'epitope':<12}{'new rows':>10}{'D now':>9}{'D + Atlas':>11}{'diff':>9}")
    rows = []
    for e in pub:
        add = len(extra_new[e]) - len(extra_old[e])
        print(f"  {e:<12}{add:10d}{po[e]:9.4f}{pn[e]:11.4f}{pn[e]-po[e]:+9.4f}")
        rows.append([e, add, po[e], pn[e], pn[e] - po[e]])
    mo = float(np.mean([po[e] for e in pub]))
    mn = float(np.mean([pn[e] for e in pub]))
    print(f"  {'MACRO':<12}{'':10}{mo:9.4f}{mn:11.4f}{mn-mo:+9.4f}")

    print("\n=== the private side: unmeasurable, so stated as a bet ===")
    print(f"  {'epitope':<12}{'pool before':>13}{'pool after':>12}"
          f"{'  E84 dose-response says':<28}")
    bet = []
    for e in ("CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY"):
        b, a_ = len(extra_old[e]), len(extra_new[e])
        verdict = ("still dead (<60)" if a_ < 60 else
                   "crosses into the live range")
        print(f"  {e:<12}{b:13d}{a_:12d}  {verdict}")
        bet.append([e, b, a_, verdict])

    sub = pd.DataFrame({"tcr_id": test["tcr_id"],
                        "prediction": d_new.reindex(test["tcr_id"]).to_numpy()})
    out = ROOT / "submissions/S7_D_atlas_pool.csv"
    sub.to_csv(out, index=False)
    import hashlib
    print(f"\nwrote {out.name}  sha256 "
          f"{hashlib.sha256(out.read_bytes()).hexdigest()[:16]}")

    lines = ["# E122 — The TCR Motif Atlas, audited where it matters", "",
             "Generated by `scripts/108_atlas_pool_expansion.py`.", "",
             "E84 measured that the pool term vanishes below ~60 same-epitope",
             "rows and the private epitopes hold 11, 191 and 5, so the only",
             "question for a new source is whether it supplies rows **there**.",
             "", "## New records after de-duplication", "",
             "| epitope | split | pool before | pool after | new |",
             "|---|---|---:|---:|---:|"]
    for e in eps:
        sp = test[test["epitope_aa"] == e]["split"].mode().iloc[0]
        lines.append(f"| {e} | {sp} | {len(extra_old[e])} | {len(extra_new[e])} "
                     f"| **{len(extra_new[e])-len(extra_old[e])}** |")
    lines += ["", "All 128 LTDEMIAQY additions come from PMID:40628266 (Drost et",
              "al., *Cell Genomics* 2025), which predates the IMMREP26",
              "revalidation (PMID:42822443, 2026) — reference receptors, not",
              "validated private rows.", "",
              "## Effect on the public epitopes, which is what can be measured", "",
              "| epitope | new rows | D | D + Atlas | diff |",
              "|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:.4f} | {r[3]:.4f} | {r[4]:+.4f} |")
    lines.append(f"| **MACRO** | | {mo:.4f} | {mn:.4f} | **{mn-mo:+.4f}** |")
    lines += ["", "## The private side, stated as a bet rather than a result", "",
              "| epitope | pool before | pool after | E84's dose-response |",
              "|---|---:|---:|---|"]
    for b in bet:
        lines.append(f"| {b[0]} | {b[1]} | {b[2]} | {b[3]} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
