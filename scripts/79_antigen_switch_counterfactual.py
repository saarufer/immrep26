#!/usr/bin/env python
"""The antigen-switch counterfactual: can the density mechanism pick the right
peptide, or only the right TCR?

Proposed by a collaborator, and it isolates the one claim this project has
never tested. Everything measured so far is **task A, revalidation**: given a
(TCR, peptide) annotation, is it real? D does 0.860 at that. **Task B,
specificity** is a different question: given one TCR and several peptides,
which one is the real target? A model can be excellent at A and useless at B,
because A can be won by asking only "is this TCR typical of the TCRs that were
validated against this peptide" -- a question about the TCR and the cohort, not
about the peptide.

So for one TCR t, compute the density against every epitope's reference cohort,

    z_GIL(t), z_GLC(t), z_NLV(t), ..., z_LTD(t)

and ask whether the annotated epitope comes out on top.

**The negatives caveat, which the proposal is right to insist on.** For the
seven non-annotated epitopes we do not know that t is a non-binder -- it was
never tested. That contaminates the measurement in one direction only: genuine
cross-reactivity can only make a real signal *look worse*. So a positive result
is unambiguous, and a chance-level result is nearly so, since cross-reactivity
across peptides from five different viruses would have to be implausibly common
to drag a real signal all the way to chance. The rate is estimated from VDJdb
rather than assumed. A within-TCR experimental panel would be better and does
not exist here: the competition data holds exactly **one** alpha-beta pair under
two epitopes (negative under both), and the PRP panel (16 peptides x 21 TCRs,
with real tested negatives) has no reference cohort for any of its peptides, so
no density can be computed against them.

**The design that makes the result interpretable is the exclusion ladder.**
A TCR sits inside its own epitope's cohort, so a naive version wins by finding
the query's own clonal relatives and says nothing about peptides. Three levels:

  (i)   exclude the query only
  (ii)  exclude the query and its single-linkage family (SCEPTR cosine <= 0.12)
  (iii) exclude the query, its family, and every TCR sharing a donor with it

If accuracy collapses down the ladder, then what looks like antigen specificity
is clonal and donor proximity -- which is E55's and E76's finding arriving from
a third direction. If it survives (iii), the mechanism really does carry
peptide-specific information and this project has been underselling it.

Every cohort is subsampled to a common size so that a large reference set
cannot win by density alone, and a label-permutation null gives the chance
level empirically rather than assuming 1/8.

The structural half of this question is already answered: E36 folded 300 TCRs
against anchor-matched decoy peptides (2,706 folds including the private
extension) and every target-minus-decoy margin collapsed to chance. This is the
density half.

Run: .venv/bin/python scripts/79_antigen_switch_counterfactual.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/54_antigen_switch_counterfactual.md"
VDJDB_FULL = ROOT / "data/external/vdjdb/vdjdb-2026-06-03/vdjdb_full.txt"
FAMILY_RADIUS = 0.12
SIGMA = 0.1
N_SUBSAMPLE = 110          # below the smallest cohort (LTDEMIAQY, 130) after exclusions
REPEATS = 30
SEED = 0


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    V = cached["vectors"]
    V = V / np.linalg.norm(V, axis=1, keepdims=True)
    pos_of = {t: i for i, t in enumerate(cached["tcr_id"])}

    ids = test["tcr_id"].to_numpy()
    ep = test["epitope_aa"].to_numpy()
    idx = np.array([pos_of[t] for t in ids])
    X = V[idx]
    Dm = 1.0 - X @ X.T                      # full 2636 x 2636 cosine distance
    np.fill_diagonal(Dm, 0.0)
    epitopes = sorted(set(ep))
    members = {e: np.where(ep == e)[0] for e in epitopes}
    print(f"{len(ids)} TCRs, {len(epitopes)} epitopes: "
          + ", ".join(f"{e[:4]} {len(members[e])}" for e in epitopes))

    # ---------- the exclusion ladder ----------
    # family: single-linkage components over the WHOLE set, so a query's
    # relatives are removed from every cohort, not only its own.
    _, fam = connected_components(sp.csr_matrix(Dm <= FAMILY_RADIUS),
                                  directed=False)
    print(f"single-linkage families over all 2636 at <= {FAMILY_RADIUS}: "
          f"{len(set(fam))} components, largest {np.bincount(fam).max()}")

    print("resolving donors from the VDJdb full table ...")
    cross = load_module("cross", "scripts/55_cross_donor_persistence.py")
    donor = cross.donor_sets(test)
    donor_list = [donor[t] for t in ids]
    have_donor = np.array([len(d) > 0 for d in donor_list])
    print(f"  donor resolved for {have_donor.mean():.1%} of TCRs")

    # Donor = reference.id + meta.subject.id, so excluding donors does NOT
    # exclude two subjects from the same publication -- same protocol, same
    # primers, same sorting gate. That is exactly the confound that could
    # produce a cohort-identification signal with no peptide information in it,
    # so the study is resolved separately and used as a fourth rung.
    vfull = pd.read_csv(VDJDB_FULL, sep="\t", low_memory=False)
    vfull = vfull[(vfull["species"] == "HomoSapiens")
                  & vfull["cdr3.alpha"].notna() & vfull["cdr3.beta"].notna()]
    ref_grouped = (vfull.groupby(["cdr3.alpha", "cdr3.beta", "antigen.epitope"])
                   ["reference.id"].agg(lambda s: frozenset(s.astype(str))))
    study_list = [ref_grouped.get((a, b, e), frozenset()) for a, b, e in
                  zip(test["cdr3_alpha_aa"], test["cdr3_beta_aa"],
                      test["epitope_aa"])]
    print(f"  study  resolved for "
          f"{np.mean([len(x) > 0 for x in study_list]):.1%} of TCRs")
    hla = test["hla_short"].to_numpy()
    print(f"  HLA: " + ", ".join(f"{h} {int((hla == h).sum())}"
                                 for h in sorted(set(hla))))

    # ---------- the queries ----------
    lab = labels.set_index("tcr_id")
    y = pd.Series(ids).map(lab["pos"]).to_numpy()
    query = np.where(y == 1)[0]              # validated binders only
    print(f"\nqueries: {len(query)} validated (pos=1) TCRs across "
          f"{len(set(ep[query]))} epitopes")

    # ---------- the cross-reactivity rate, from VDJdb ----------
    print("\n=== how often is one TCR annotated to more than one of the eight? ===")
    full = pd.read_csv(ROOT / "data/external/vdjdb/vdjdb-2026-06-03/vdjdb_full.txt",
                       sep="\t", low_memory=False)
    full = full[(full["species"] == "HomoSapiens")
                & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]
    sub = full[full["antigen.epitope"].isin(epitopes)]
    per_tcr = sub.groupby(["cdr3.alpha", "cdr3.beta"])["antigen.epitope"].nunique()
    xr = float((per_tcr > 1).mean())
    print(f"  {len(per_tcr)} distinct alpha-beta pairs annotated to any of the eight")
    print(f"  annotated to more than one: {int((per_tcr > 1).sum())} ({xr:.2%})")
    print("  -> the ceiling on how much 'never tested' can be contaminating")

    rng = np.random.default_rng(SEED)

    def run(level: str, permute: bool = False, panel=None,
            drop_nn: int = 0, Dmat=None) -> dict:
        """Mean over REPEATS subsamples of: argmax accuracy and true-epitope rank.

        `panel` restricts the candidate epitopes (and hence the queries) so that
        MHC restriction can be held fixed instead of being allowed to answer the
        question for us."""
        panel = list(epitopes) if panel is None else list(panel)
        Dl = Dm if Dmat is None else Dmat
        qmask = np.isin(ep[query], panel) | np.array(
            [e in panel for e in ep[query]])
        q = query[np.isin(ep[query], panel)] if len(panel) < len(epitopes) else query
        true_ep = ep[q].copy()
        if permute:
            true_ep = rng.permutation(true_ep)
        hits, ranks, per_ep_hits = [], [], {e: [] for e in panel}
        surviving = {j: [] for j in range(len(panel))}
        pred_marginal = {e: [] for e in panel}
        for _ in range(REPEATS):
            # one common subsample of each cohort for this repeat
            pool = {e: rng.choice(members[e], size=min(N_SUBSAMPLE, len(members[e])),
                                  replace=False) for e in panel}
            Zs = np.zeros((len(q), len(panel)))
            for j, e in enumerate(panel):
                ref = pool[e]
                # queries against this cohort, with the ladder applied
                dq_ref = Dl[np.ix_(q, ref)]
                kq = np.exp(-dq_ref / SIGMA)
                keep = np.ones_like(kq, dtype=bool)
                keep &= q[:, None] != ref[None, :]                      # (i)
                if level in ("family", "donor", "study"):
                    keep &= fam[q][:, None] != fam[ref][None, :]        # (ii)
                if level in ("donor", "study"):                         # (iii)
                    for a, qi in enumerate(q):
                        dq = donor_list[qi]
                        if not dq:
                            continue
                        for b, ri in enumerate(ref):
                            if dq & donor_list[ri]:
                                keep[a, b] = False
                if drop_nn:
                    # Symmetric excision: remove the query's drop_nn nearest
                    # TCRs in EVERY cohort, not only its own. This does not
                    # depend on a family radius, a donor field or a citation,
                    # and it removes the "my own relatives are sitting in my
                    # cohort" route by construction.
                    cut = np.argsort(dq_ref, axis=1)[:, :drop_nn]
                    np.put_along_axis(keep, cut, False, axis=1)
                if level == "study":                                    # (iv)
                    for a, qi in enumerate(q):
                        sq = study_list[qi]
                        if not sq:
                            continue
                        for b, ri in enumerate(ref):
                            if sq & study_list[ri]:
                                keep[a, b] = False
                n = keep.sum(axis=1)
                dens = np.where(n > 0, (kq * keep).sum(axis=1) / np.maximum(n, 1),
                                np.nan)
                Zs[:, j] = dens
                surviving[j].append(float(n.mean()))
            # Standardise each column over the QUERY set, not over the cohort's
            # own members. The first version standardised against each cohort's
            # internal leave-one-out densities, which is not symmetric -- cohort
            # members have their relatives present while queries have theirs
            # excluded by the ladder -- and it produced a degenerate predictor
            # that assigned 62% of everything to NLVPMVATV. Per-column
            # standardisation over the same 740-910 queries removes every
            # cohort-level scale difference by construction.
            Zs = (Zs - np.nanmean(Zs, axis=0)) / (np.nanstd(Zs, axis=0) + 1e-12)
            order = np.argsort(-Zs, axis=1)
            picked = np.array(panel)[order[:, 0]]
            hits.append(float((picked == true_ep).mean()))
            rk = np.array([list(np.array(panel)[order[a]]).index(true_ep[a]) + 1
                           for a in range(len(q))])
            ranks.append(float(rk.mean()))
            for e in panel:
                m = true_ep == e
                if m.any():
                    per_ep_hits[e].append(float((picked[m] == e).mean()))
                pred_marginal[e].append(float((picked == e).mean()))
        return {"acc": float(np.mean(hits)), "acc_sd": float(np.std(hits)),
                "rank": float(np.mean(ranks)),
                "per_ep": {e: float(np.mean(v)) if v else float("nan")
                           for e, v in per_ep_hits.items()},
                "n_query": int(len(q)), "n_panel": len(panel),
                "chance": 1.0 / len(panel),
                "pred": {e: float(np.mean(v)) for e, v in pred_marginal.items()},
                "surv": {panel[j]: float(np.mean(v))
                         for j, v in surviving.items()}}

    print(f"\n=== the exclusion ladder ({REPEATS} subsamples of "
          f"{N_SUBSAMPLE} per cohort) ===")
    print(f"  chance argmax = 1/{len(epitopes)} = {1/len(epitopes):.3f}, "
          f"chance mean rank = {(len(epitopes)+1)/2:.1f}")
    print(f"  {'exclusion level':<34}{'argmax acc':>12}{'± sd':>8}"
          f"{'mean rank of true':>20}")
    results = {}
    for level, label in (("self", "(i)   query only"),
                         ("family", "(ii)  + single-linkage family"),
                         ("donor", "(iii) + same donor"),
                         ("study", "(iv)  + same publication")):
        r = run(level)
        results[level] = r
        print(f"  {label:<34}{r['acc']:12.3f}{r['acc_sd']:8.3f}{r['rank']:20.2f}")
    null = run("study", permute=True)
    results["null"] = null
    print(f"  {'permuted-label null (level iv)':<34}{null['acc']:12.3f}"
          f"{null['acc_sd']:8.3f}{null['rank']:20.2f}")

    print(f"\n  how degenerate is the predictor? "
          f"(share of queries assigned to each epitope, level iv)")
    print("  " + "  ".join(f"{e[:4]} {results['study']['pred'][e]:.3f}"
                           for e in epitopes))
    print(f"\n  reference TCRs surviving the ladder, per (query, cohort) pair:")
    print(f"  {'epitope':<12}{'(i)':>8}{'(ii)':>8}{'(iii)':>8}{'(iv)':>8}")
    for e in epitopes:
        print(f"  {e:<12}" + "".join(f"{results[k]['surv'][e]:8.1f}"
                                     for k in ("self", "family", "donor", "study")))

    print(f"\n  per-epitope argmax accuracy:")
    print(f"  {'epitope':<12}" + "".join(f"{k:>10}" for k in
                                         ("(i)", "(ii)", "(iii)", "(iv)", "null")))
    for e in epitopes:
        if np.isnan(results["self"]["per_ep"][e]):
            continue
        print(f"  {e:<12}" + "".join(
            f"{results[k]['per_ep'][e]:10.3f}"
            for k in ("self", "family", "donor", "study", "null")))

    # ---------- MHC restriction held fixed ----------
    # TTDPSFLGRY and LTDEMIAQY are the only HLA-A*01:01 cohorts. V-gene usage
    # differs by restriction element, so "pick the right epitope" can be won by
    # picking the right HLA without any peptide information. Rerunning on the
    # six A*02:01 epitopes alone removes that route entirely.
    a2 = sorted({e for e, h in zip(ep, hla) if h == "A2"})
    print(f"\n=== MHC restriction held fixed: the {len(a2)} HLA-A2 epitopes only ===")
    print("  " + ", ".join(a2))
    hla_res = {}
    for level, label in (("self", "(i)   query only"),
                         ("family", "(ii)  + single-linkage family"),
                         ("study", "(iv)  + donor + publication")):
        hla_res[level] = run(level, panel=a2)
    hla_res["null"] = run("study", permute=True, panel=a2)
    print(f"  chance argmax = {hla_res['self']['chance']:.3f}, "
          f"{hla_res['self']['n_query']} queries")
    print(f"  {'exclusion level':<34}{'argmax acc':>12}{'± sd':>8}"
          f"{'mean rank of true':>20}")
    for level, label in (("self", "(i)   query only"),
                         ("family", "(ii)  + single-linkage family"),
                         ("study", "(iv)  + donor + publication"),
                         ("null", "permuted-label null")):
        r = hla_res[level]
        print(f"  {label:<34}{r['acc']:12.3f}{r['acc_sd']:8.3f}{r['rank']:20.2f}")
    print(f"\n  per-epitope, HLA-A2 panel:")
    print(f"  {'epitope':<12}" + "".join(f"{k:>10}" for k in
                                         ("(i)", "(ii)", "(iv)", "null")))
    for e in a2:
        if np.isnan(hla_res["self"]["per_ep"][e]):
            continue
        print(f"  {e:<12}" + "".join(f"{hla_res[k]['per_ep'][e]:10.3f}"
                                     for k in ("self", "family", "study", "null")))

    # ---------- the symmetric nearest-neighbour excision sweep ----------
    print("\n=== removing the query's m nearest TCRs from EVERY cohort ===")
    print("  Symmetric by construction, so it cannot favour any epitope, and it")
    print("  needs no family radius, donor field or citation to be populated.")
    print(f"  {'m removed per cohort':<24}{'argmax acc':>12}{'null':>8}"
          f"{'ratio':>8}{'mean rank':>11}")
    sweep = []
    for m in (0, 5, 10, 20, 40):
        r = run("study", panel=a2, drop_nn=m)
        rn = run("study", permute=True, panel=a2, drop_nn=m)
        print(f"  {m:<24d}{r['acc']:12.3f}{rn['acc']:8.3f}"
              f"{r['acc']/rn['acc']:8.2f}{r['rank']:11.2f}")
        sweep.append([m, r["acc"], rn["acc"], r["acc"] / rn["acc"], r["rank"]])

    # ---------- is it V-gene usage or CDR3 sequence? ----------
    # Epitope-biased V-gene usage is well documented and would produce a
    # cohort-identification signal that is real but is not CDR3-level
    # recognition. cdr3_only drops the V/J genes from the representation.
    print("\n=== V-gene usage or CDR3 sequence? (A2 panel, level iv, m=20) ===")
    print(f"  {'representation':<16}{'argmax acc':>12}{'null':>8}{'ratio':>8}")
    rep_rows = []
    for rep in ("default", "cdr3_only", "a_sceptr", "b_sceptr"):
        c = np.load(ROOT / f"data/interim/sceptr/{rep}.npz", allow_pickle=True)
        W = c["vectors"]
        W = W / np.linalg.norm(W, axis=1, keepdims=True)
        po = {t: i for i, t in enumerate(c["tcr_id"])}
        Xr = W[[po[t] for t in ids]]
        Dr = 1.0 - Xr @ Xr.T
        np.fill_diagonal(Dr, 0.0)
        r = run("study", panel=a2, drop_nn=20, Dmat=Dr)
        rn = run("study", permute=True, panel=a2, drop_nn=20, Dmat=Dr)
        print(f"  {rep:<16}{r['acc']:12.3f}{rn['acc']:8.3f}"
              f"{r['acc']/rn['acc']:8.2f}")
        rep_rows.append([rep, r["acc"], rn["acc"], r["acc"] / rn["acc"]])

    lines = ["# E88 — The antigen-switch counterfactual: specificity, not revalidation",
             "",
             "Proposed by a collaborator. Everything measured before this is **task A,",
             "revalidation** (is this annotation real?), where D scores 0.860. This is",
             "**task B, specificity** (which peptide is the real target?).",
             "Generated by `scripts/79_antigen_switch_counterfactual.py`.", "",
             "## Why the within-TCR panel does not exist here", "",
             "The competition data holds exactly **one** alpha-beta pair under two",
             "epitopes, negative under both, so no experimentally-tested within-TCR",
             "switch is available. The PRP panel (16 peptides x 21 TCRs) has real",
             "tested negatives but no reference cohort for any of its peptides, so no",
             "density can be computed against them. The eight competition cohorts are",
             "therefore used as the panel, and the 'never tested' contamination is",
             "bounded rather than ignored:", "",
             f"- {len(per_tcr)} distinct alpha-beta pairs in VDJdb are annotated to any",
             f"  of the eight; **{int((per_tcr > 1).sum())} ({xr:.2%})** to more than one.",
             "- Contamination runs one way only: real cross-reactivity makes a true",
             "  signal look *worse*, never better. A positive result is unambiguous.",
             "",
             "## The exclusion ladder", "",
             f"{REPEATS} subsamples of {N_SUBSAMPLE} TCRs per cohort, so a large",
             "reference set cannot win on size. Chance argmax = "
             f"{1/len(epitopes):.3f}; the permuted-label null measures it empirically.",
             "",
             "| exclusion level | argmax accuracy | ± sd | mean rank of true epitope |",
             "|---|---:|---:|---:|"]
    for level, label in (("self", "(i) query only"),
                         ("family", "(ii) + single-linkage family"),
                         ("donor", "(iii) + same donor"),
                         ("null", "permuted-label null (level iii)")):
        r = results[level]
        lines.append(f"| {label} | **{r['acc']:.3f}** | {r['acc_sd']:.3f} "
                     f"| {r['rank']:.2f} |")
    lines += ["", "Share of queries assigned to each epitope at level (iv) -- a",
              "degenerate predictor would concentrate here:", "",
              "| " + " | ".join(e[:4] for e in epitopes) + " |",
              "|" + "---:|" * len(epitopes),
              "| " + " | ".join(f"{results['study']['pred'][e]:.3f}"
                                for e in epitopes) + " |", "",
              "Reference TCRs surviving the ladder, per (query, cohort) pair "
              f"(out of {N_SUBSAMPLE}):", "",
              "| epitope | (i) | (ii) | (iii) | (iv) |", "|---|---:|---:|---:|---:|"]
    for e in epitopes:
        lines.append(f"| {e} | " + " | ".join(
            f"{results[k]['surv'][e]:.1f}"
            for k in ("self", "family", "donor", "study")) + " |")
    lines += ["", "Per-epitope argmax accuracy:", "",
              "| epitope | (i) | (ii) | (iii) | null |", "|---|---:|---:|---:|---:|"]
    for e in epitopes:
        if np.isnan(results["self"]["per_ep"][e]):
            continue
        lines.append(f"| {e} | " + " | ".join(
            f"{results[k]['per_ep'][e]:.3f}"
            for k in ("self", "family", "donor", "null")) + " |")
    lines.append("")

    lines += ["## MHC restriction held fixed", "",
              "TTDPSFLGRY and LTDEMIAQY are the only HLA-A\\*01:01 cohorts, and V-gene",
              "usage differs by restriction element -- so the panel above can be won",
              "by identifying the HLA, with no peptide information at all. Rerun on",
              f"the {len(a2)} HLA-A\\*02:01 epitopes alone "
              f"({hla_res['self']['n_query']} queries, chance "
              f"{hla_res['self']['chance']:.3f}):", "",
              "| exclusion level | argmax accuracy | ± sd | mean rank |",
              "|---|---:|---:|---:|"]
    for level, label in (("self", "(i) query only"),
                         ("family", "(ii) + family"),
                         ("study", "(iv) + donor + publication"),
                         ("null", "permuted-label null")):
        r = hla_res[level]
        lines.append(f"| {label} | **{r['acc']:.3f}** | {r['acc_sd']:.3f} "
                     f"| {r['rank']:.2f} |")
    lines += ["", "| epitope | (i) | (ii) | (iv) | null |", "|---|---:|---:|---:|---:|"]
    for e in a2:
        if np.isnan(hla_res["self"]["per_ep"][e]):
            continue
        lines.append(f"| {e} | " + " | ".join(
            f"{hla_res[k]['per_ep'][e]:.3f}"
            for k in ("self", "family", "study", "null")) + " |")
    lines.append("")

    lines += ["## The symmetric nearest-neighbour excision", "",
              "The ladder above is weak in absolute terms -- it removes only 10-20 of",
              "110 reference TCRs. This removes the query's `m` nearest TCRs from",
              "*every* cohort, which is symmetric by construction and needs no family",
              "radius, donor field or citation to be populated. HLA-A2 panel, level",
              "(iv).", "",
              "| m removed per cohort | argmax accuracy | null | ratio | mean rank |",
              "|---:|---:|---:|---:|---:|"]
    for r in sweep:
        lines.append(f"| {r[0]} | **{r[1]:.3f}** | {r[2]:.3f} | {r[3]:.2f} "
                     f"| {r[4]:.2f} |")
    lines += ["", "## Is it V-gene usage or CDR3 sequence?", "",
              "Epitope-biased V-gene usage is well documented and would produce a",
              "real cohort-identification signal that is nonetheless not CDR3-level",
              "recognition. `cdr3_only` drops the V and J genes from the",
              "representation. A2 panel, level (iv), m=20.", "",
              "| representation | argmax accuracy | null | ratio |",
              "|---|---:|---:|---:|"]
    for r in rep_rows:
        lines.append(f"| {r[0]} | **{r[1]:.3f}** | {r[2]:.3f} | {r[3]:.2f} |")
    lines.append("")

    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    np.save(ROOT / "data/interim/antigen_switch.npy",
            np.array([results], dtype=object), allow_pickle=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
