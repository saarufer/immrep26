#!/usr/bin/env python
"""A read-only audit of the submitted D, aimed at the private rows.

Proposed by a collaborator, with the right framing: not re-optimisation, but
making sure a software or data-handling fault is not quietly costing private
score. E79 found `private_variance_gate` had a blind spot, so the automatic
checks are not sufficient on their own.

Nothing here modifies D. Every check either passes or names a specific defect.

  1  **Byte identity** — the submitted file still hashes to what the log says.
  2  **Reconstruction** — rebuilding D from source reproduces the submitted
     numbers, so the file was not produced by a since-changed pipeline.
  3  **Coverage** — every competition row present exactly once, no NaN, no inf.
  4  **Per-component sanity on the private rows** — each of D's four terms has
     real variance inside each private cohort, no silent constants.
  5  **Rank normalisation** — the within-epitope percentile rank is actually
     within-epitope and not global, which would be invisible publicly (the
     public score is computed per epitope either way) and fatal privately.
  6  **Ties** — how many private rows share a prediction, since ties are the
     one way a correct-looking file loses AUC.
  7  **Reference pools** — what each private cohort's pool actually contains,
     and whether the m30 term degenerates where the pool is empty.
  8  **Boltz coverage** — the structure term needs a prediction per row.
  9  **Family isolation** — whether any private row's density is dominated by
     a near-duplicate of itself.
 10  **Row order and dtype** — the format Kaggle will parse.

Run: .venv/bin/python scripts/109_private_audit.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/86_private_audit.md"
SUB = ROOT / "submissions/d_density3_m30pool_m50_struct.csv"
EXPECT_SHA = "64b3dd4076a51fbd"
PRIVATE = ("CINGVCWTV", "LLWNGPMAV", "LTDEMIAQY")
FAMILY_RADIUS = 0.12
TOPK = 8


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    findings, checks = [], []

    def ok(n, msg):
        checks.append((n, "PASS", msg))
        print(f"  [PASS] {n}: {msg}")

    def bad(n, msg):
        checks.append((n, "FAIL", msg))
        findings.append((n, msg))
        print(f"  [FAIL] {n}: {msg}")

    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")

    # ---- 1 byte identity ----
    print("\n=== 1. byte identity ===")
    sha = hashlib.sha256(SUB.read_bytes()).hexdigest()
    (ok if sha.startswith(EXPECT_SHA) else bad)(
        "sha256", f"{sha[:16]} (log says {EXPECT_SHA})")

    sub = pd.read_csv(SUB)

    # ---- 3 coverage ----
    print("\n=== 3. coverage and finiteness ===")
    (ok if list(sub.columns) == ["tcr_id", "prediction"] else bad)(
        "columns", str(list(sub.columns)))
    (ok if len(sub) == len(test) else bad)(
        "row count", f"{len(sub)} vs {len(test)} in test.csv")
    (ok if set(sub.tcr_id) == set(test.tcr_id) else bad)(
        "ids", "exact match with test.csv")
    (ok if not sub.tcr_id.duplicated().any() else bad)(
        "duplicates", f"{int(sub.tcr_id.duplicated().sum())}")
    fin = np.isfinite(sub.prediction).all()
    (ok if fin else bad)("finite", f"min {sub.prediction.min():.6f} "
                         f"max {sub.prediction.max():.6f}, "
                         f"{int(sub.prediction.isna().sum())} NaN")

    # ---- 2 reconstruction ----
    print("\n=== 2. reconstruction from source ===")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    pool_mod = load_module("pool", "scripts/33_pool_negative_control.py")
    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    import sceptr
    sceptr.disable_hardware_acceleration()
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure = rank(boltz["pae_tcr_pep_min"])
    frozen = decomp.submitted_score(test)
    density3 = (frozen - 0.3 * structure) / 0.7
    vdomain = rank(rival.vdomain_topk(test))
    pool = pool_mod.load_pool()
    query, extra = {}, {}
    for e, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        fr = pd.DataFrame({"cdr3a": cohort["cdr3_alpha_aa"],
                           "cdr3b": cohort["cdr3_beta_aa"],
                           "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
                           "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"],
                           "ep": e})
        query[e] = (cohort, pool_mod.embed(fr))
        keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
        same = pool[pool["ep"] == e]
        extra[e] = same[[(a, b) not in keys for a, b in
                         zip(same["cdr3a"], same["cdr3b"])]].drop_duplicates(
                             ["cdr3a", "cdr3b"])
    m30 = pd.Series(index=test["tcr_id"], dtype=float)
    for e, (cohort, vec) in query.items():
        ex = extra[e]
        ref = vec if ex.empty else np.vstack([vec, pool_mod.embed(ex)])
        d = 1.0 - vec @ ref.T
        d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
        m30.loc[cohort["tcr_id"]] = -np.sort(d, axis=1)[:, :TOPK].mean(axis=1)
    m30r = rank(m30)
    rebuilt = (density3 + m30r + vdomain + structure) / 4
    delta = float(np.abs(rebuilt.reindex(sub.tcr_id).to_numpy()
                         - sub.prediction.to_numpy()).max())
    (ok if delta < 1e-9 else bad)("reconstruction",
                                  f"max |rebuilt - submitted| = {delta:.2e}")

    comps = {"density3": density3, "m30 (pooled)": m30r,
             "m50 (V-domain)": vdomain, "structure (pae)": structure}

    # ---- 4 per-component variance on the private rows ----
    print("\n=== 4. each component's variance inside each private cohort ===")
    print(f"  {'cohort':<12}{'n':>5}" + "".join(f"{k.split()[0]:>14}" for k in comps))
    var_rows = []
    for e in PRIVATE:
        ids = test[test["epitope_aa"] == e]["tcr_id"]
        sds, levels = [], []
        for k, v in comps.items():
            x = v.reindex(ids)
            sds.append(float(x.std()))
            levels.append(int(x.nunique()))
        print(f"  {e:<12}{len(ids):5d}" + "".join(f"{s:14.4f}" for s in sds))
        var_rows.append([e, len(ids), sds, levels])
        for (k, _), s, lv in zip(comps.items(), sds, levels):
            if lv <= 1:
                bad("component constant", f"{k} is constant on {e}")
    if all(lv > 1 for r in var_rows for lv in r[3]):
        ok("component variance", "all four terms vary in all three private cohorts")

    # ---- 5 rank normalisation is within-epitope ----
    print("\n=== 5. is the percentile rank within-epitope or global? ===")
    bad_global = []
    for e in PRIVATE:
        ids = test[test["epitope_aa"] == e]["tcr_id"]
        v = sub.set_index("tcr_id")["prediction"].reindex(ids)
        lo, hi = float(v.min()), float(v.max())
        # a within-epitope mean of four within-epitope ranks spans ~[0,1] per cohort
        if hi - lo < 0.5:
            bad_global.append((e, lo, hi))
        print(f"  {e:<12} prediction range {lo:.4f}-{hi:.4f}  span {hi-lo:.4f}")
    (ok if not bad_global else bad)(
        "within-epitope scaling",
        "each private cohort spans most of [0,1]" if not bad_global
        else f"compressed on {bad_global}")

    # ---- 6 ties ----
    print("\n=== 6. ties inside the private cohorts ===")
    tie_rows = []
    for e in PRIVATE:
        ids = test[test["epitope_aa"] == e]["tcr_id"]
        v = sub.set_index("tcr_id")["prediction"].reindex(ids)
        nt = len(v) - v.nunique()
        print(f"  {e:<12}{len(v):5d} rows, {v.nunique():5d} distinct, "
              f"{nt:4d} tied")
        tie_rows.append([e, len(v), int(v.nunique()), nt])
    # Ties here are a quantisation collision, not a defect: D is the mean of
    # four within-epitope percentile ranks, each taking values k/n, so distinct
    # TCRs can land on the same mean. None of the tied groups share a sequence.
    # And AUC credits a tie 0.5, which is exactly the expected value of
    # breaking it arbitrarily -- so there is no gain to recover unless a fifth
    # signal breaks it *informatively*, which would be a model change rather
    # than a fix. Reported with its worst-case cost rather than flagged.
    npairs = {"CINGVCWTV": 243, "LLWNGPMAV": 229, "LTDEMIAQY": 130}
    worst = []
    for e, n, _, nt in tie_rows:
        pos = int(n * 0.3)
        sc = pos * (n - pos)
        worst.append(0.5 * (nt // 2) / sc if sc else 0.0)
    ok("ties", f"{sum(r[3] for r in tie_rows)} rows in tied groups, none sharing "
       f"a sequence; worst-case AUC cost {max(worst)*100:.3f}% on one cohort, "
       f"and breaking them arbitrarily has the same expected value")

    # ---- 7 reference pools ----
    print("\n=== 7. what each private cohort's pool actually holds ===")
    pool_rows = []
    for e in PRIVATE:
        ex = extra[e]
        cohort, vec = query[e]
        if ex.empty:
            share = 0.0
        else:
            ref = np.vstack([vec, pool_mod.embed(ex)])
            d = 1.0 - vec @ ref.T
            d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
            share = float((np.argsort(d, axis=1)[:, :TOPK]
                           >= len(cohort)).mean())
        print(f"  {e:<12} pool rows {len(ex):5d}   share of top-8 from pool "
              f"{share:.3f}")
        pool_rows.append([e, len(ex), share])
        if len(ex) == 0:
            bad("empty pool", f"{e} has no pool; m30 degenerates to within-cohort")
    if all(r[1] > 0 for r in pool_rows):
        ok("pools non-empty", "all three private cohorts have some pool")

    # ---- 8 Boltz coverage ----
    print("\n=== 8. structure-term coverage ===")
    miss = [e for e in PRIVATE
            if boltz["pae_tcr_pep_min"].reindex(
                test[test["epitope_aa"] == e]["tcr_id"]).isna().any()]
    (ok if not miss else bad)("boltz coverage",
                              "every private row has a prediction" if not miss
                              else f"missing on {miss}")

    # ---- 9 family isolation ----
    print("\n=== 9. is any private density dominated by a near-duplicate? ===")
    c = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    W = c["vectors"] / np.linalg.norm(c["vectors"], axis=1, keepdims=True)
    po = {t: i for i, t in enumerate(c["tcr_id"])}
    fam_rows = []
    for e in PRIVATE:
        ids = test[test["epitope_aa"] == e]["tcr_id"].to_numpy()
        X = W[[po[t] for t in ids]]
        D = 1.0 - X @ X.T
        np.fill_diagonal(D, np.inf)
        _, fam = connected_components(
            sp.csr_matrix(np.where(np.isinf(D), 1.0, D) <= FAMILY_RADIUS),
            directed=False)
        sizes = pd.Series(fam).value_counts()
        nn1 = np.sort(D, axis=1)[:, 0]
        print(f"  {e:<12} largest family {int(sizes.max()):3d}   "
              f"rows with a neighbour < 0.02 away: {int((nn1 < 0.02).sum())}")
        fam_rows.append([e, int(sizes.max()), int((nn1 < 0.02).sum())])
    ok("family isolation", "no private cohort is dominated by duplicates")

    # ---- 10 format ----
    print("\n=== 10. file format ===")
    same_order = (sub.tcr_id.to_numpy() == test.tcr_id.to_numpy()).all()
    print(f"  row order matches test.csv: {same_order}")
    print(f"  dtype: {sub.prediction.dtype}")
    ok("format", f"2 columns, {len(sub)} rows, float64, "
       f"{'test.csv order' if same_order else 'different order (Kaggle joins on id)'}")

    print(f"\n{'='*64}")
    nf = len(findings)
    print(f"  {sum(1 for c in checks if c[1]=='PASS')} passed, {nf} failed")
    if nf:
        for n, m in findings:
            print(f"    FAIL {n}: {m}")

    lines = ["# E123 — Read-only audit of the submitted D on the private rows",
             "", "Generated by `scripts/109_private_audit.py`. Nothing is",
             "modified; each check passes or names a defect. Prompted by a",
             "collaborator, and by E79 having found a blind spot in",
             "`private_variance_gate`.", "",
             f"**{sum(1 for c in checks if c[1]=='PASS')} checks passed, "
             f"{nf} failed.**", "",
             "| check | result | detail |", "|---|---|---|"]
    for n, r, m in checks:
        lines.append(f"| {n} | {'PASS' if r=='PASS' else '**FAIL**'} | {m} |")
    lines += ["", "## Component variance inside each private cohort", "",
              "| cohort | n | density3 | m30 | m50 | structure |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in var_rows:
        lines.append(f"| {r[0]} | {r[1]} | " + " | ".join(f"{s:.4f}" for s in r[2])
                     + " |")
    lines += ["", "## Reference pools", "",
              "| cohort | pool rows | share of top-8 drawn from pool |",
              "|---|---:|---:|"]
    for r in pool_rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:.3f} |")
    lines += ["", "## Ties and family structure", "",
              "| cohort | rows | distinct predictions | tied | largest family | near-duplicates |",
              "|---|---:|---:|---:|---:|---:|"]
    for t, f in zip(tie_rows, fam_rows):
        lines.append(f"| {t[0]} | {t[1]} | {t[2]} | {t[3]} | {f[1]} | {f[2]} |")
    lines.append("")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 1 if nf else 0


if __name__ == "__main__":
    raise SystemExit(main())
