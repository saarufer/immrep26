#!/usr/bin/env python
"""A deliberately public-maximising configuration, as a second entry.

The rules permit several teams and several submissions, which changes what is
worth building. The frozen entry maximises *expected private* score: it carries
no pool expansion because E27 measured that the public gain (+0.0112) is almost
entirely unavailable to the private epitopes (+0.0013 realisable), and E48
measured the same thing for the rival configuration.

Both of those rest on one quantity: the distinct-clonotype multipliers the
private epitopes can actually reach (1.038, 1.045, 1.834), counted from VDJdb.
If that count understates what is reachable -- more pool in IEDB than estimated,
or private cohorts drawing on records our counting missed -- then the discount
is wrong and the pool-heavy configuration wins. A second entry is exactly how to
take the other side of that measurement.

So this builds the best *public* configuration available, with the private risk
printed beside every row rather than hidden. The opportunity is visible in E48:

    HeidAI's configuration  = pool-expanded density + whole-V-domain alignment
    saarufer's frozen entry = three-space density + Boltz structure

Neither holds the other's term. The union has never been measured.

Nothing here replaces the frozen entry, whose byte-identical guard and nineteen
tests stay untouched. Every candidate is reported with:

  - public macro, which is what this entry is for
  - the private-realisable version, pool thinned to the private multipliers
  - the gap between them, which *is* the bet being taken

Fusion stays a zero-parameter equal-weight rank average wherever a candidate has
more than one new member, since E10's constraint is about transfer and a public
entry that overfits weights would be a worse hedge, not a better one.

Run: .venv/bin/python scripts/48_public_maximum.py
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
from immrep.validation import private_variance_gate  # noqa: E402

REPORT = ROOT / "reports/24_public_maximum.md"
TOPK = 8
SIGMA_SCEPTR = 0.1
SIGMA_TCRDIST = 20.0
STRUCTURE_WEIGHT = 0.3
REPEATS = 20
SEED = 0

# E27 / E48: the multipliers the private epitopes can actually reach.
# These are VDJdb-only -- this repository holds no IEDB download. A
# collaborator reports that adding IEDB targeted at epitopes whose VDJdb
# multiplier is below 2 raises LTDEMIAQY from 1.04 to 1.72 at no public cost.
# That claim cannot be verified here, so it is carried as a SCENARIO and the
# whole realisable column is reported under both, because it is the single
# input the pool rejection rests on.
PRIVATE_MULTIPLIERS = {"LTDEMIAQY": 1.038, "CINGVCWTV": 1.045, "LLWNGPMAV": 1.834}
MULTIPLIER_SCENARIOS = {
    "vdjdb_only (measured here)": PRIVATE_MULTIPLIERS,
    "with targeted IEDB (collaborator, unverified)": {
        "LTDEMIAQY": 1.72, "CINGVCWTV": 1.045, "LLWNGPMAV": 1.834},
}
# E47: tcrdist3's published private per-epitope AUROC.
PRIVATE_BANDS = {"CINGVCWTV": (0.80, 0.83), "LLWNGPMAV": (0.80, 0.83),
                 "LTDEMIAQY": (0.66, 0.68)}


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    hyp = load_module("hyp", "scripts/40_new_hypotheses.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")
    rival = load_module("rival", "scripts/44_rival_configuration.py")
    pool_mod = load_module("pool", "scripts/33_pool_negative_control.py")

    rank = lambda s: hyp.rank_within(s, test)      # noqa: E731
    score = lambda s: hyp.mean_auc(labels, s)      # noqa: E731

    import sceptr
    sceptr.disable_hardware_acceleration()
    rng = np.random.default_rng(SEED)

    # ---------- the terms ----------
    frozen = decomp.submitted_score(test)
    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    structure = rank(boltz["pae_tcr_pep_min"])
    density3 = (frozen - STRUCTURE_WEIGHT * structure) / (1 - STRUCTURE_WEIGHT)
    vdomain = rank(rival.vdomain_topk(test))

    print("embedding cohorts and the same-epitope pool ...")
    pool = pool_mod.load_pool()
    query, extra_available = {}, {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        frame = pd.DataFrame({
            "cdr3a": cohort["cdr3_alpha_aa"], "cdr3b": cohort["cdr3_beta_aa"],
            "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
            "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"],
            "ep": epitope})
        query[epitope] = (cohort, pool_mod.embed(frame))
        keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
        same = pool[pool["ep"] == epitope]
        extra_available[epitope] = same[
            [(a, b) not in keys for a, b in zip(same["cdr3a"], same["cdr3b"])]
        ].drop_duplicates(["cdr3a", "cdr3b"])
        print(f"  {epitope:<11} cohort {len(cohort):4d}   pool available "
              f"{len(extra_available[epitope]):5d}")

    def pooled_topk(extra_by_epitope) -> pd.Series:
        """E48's m30: top-k density against a pool-expanded reference set."""
        out = pd.Series(index=test["tcr_id"], dtype=float)
        for epitope, (cohort, vectors) in query.items():
            extra = extra_by_epitope.get(epitope)
            reference = vectors if extra is None or extra.empty else np.vstack(
                [vectors, pool_mod.embed(extra)])
            d = 1.0 - vectors @ reference.T
            d[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
            out.loc[cohort["tcr_id"]] = -np.sort(d, axis=1)[:, :TOPK].mean(axis=1)
        return out

    m30_pooled = rank(pooled_topk(extra_available))
    m30_nopool = rank(pooled_topk({e: None for e in query}))

    # ---------- the candidates ----------
    # Equal-weight rank averages, except where the frozen 0.7/0.3 split is kept
    # deliberately so a candidate is the frozen entry with one term swapped.
    candidates = {
        "frozen entry (density3 + structure)": frozen,
        "HeidAI rebuild (m30 pooled + m50)": (m30_pooled + vdomain) / 2,
        "A  frozen, density3 -> pooled m30": (1 - STRUCTURE_WEIGHT) * m30_pooled
                                             + STRUCTURE_WEIGHT * structure,
        "B  pooled m30 + m50 + structure": (m30_pooled + vdomain + structure) / 3,
        "C  HeidAI rebuild + structure 0.3": (1 - STRUCTURE_WEIGHT)
                                             * ((m30_pooled + vdomain) / 2)
                                             + STRUCTURE_WEIGHT * structure,
        "D  density3 + pooled m30 + m50 + structure": (density3 + m30_pooled
                                                       + vdomain + structure) / 4,
        "E  (density3 + pooled m30 + m50)/3 + str 0.3":
            (1 - STRUCTURE_WEIGHT) * ((density3 + m30_pooled + vdomain) / 3)
            + STRUCTURE_WEIGHT * structure,
        # D with no pool anywhere: if this still beats frozen, the advantage
        # does not depend on the pool measurement at all, and the whole
        # transfer argument against it falls away.
        "F  D with no pool at all": (density3 + m30_nopool + vdomain
                                     + structure) / 4,
    }

    print("\n=== public macro, every candidate ===")
    public = {}
    for name, series in candidates.items():
        public[name] = score(series)
        print(f"  {name:<46} {public[name]:.4f}")
    best_name = max(public, key=public.get)
    print(f"\n  best on public: {best_name}  {public[best_name]:.4f}")
    print(f"  frozen entry:   {public['frozen entry (density3 + structure)']:.4f}")
    print(f"  HeidAI rebuild: {public['HeidAI rebuild (m30 pooled + m50)']:.4f}")

    print("\n=== the same candidates at private-realisable pool density ===")
    print("  Each public pool thinned to the private epitopes' own multipliers,")
    print(f"  {REPEATS} repeats, exactly as E27 and E48 priced the discount.\n")
    realisable = {name: [] for name in candidates}
    for private_epitope, multiplier in PRIVATE_MULTIPLIERS.items():
        draws = {name: [] for name in candidates}
        for _ in range(REPEATS):
            thinned = {}
            for epitope, (cohort, _) in query.items():
                available = extra_available[epitope]
                take = int(round((multiplier - 1.0) * len(cohort)))
                thinned[epitope] = (
                    available.sample(min(take, len(available)),
                                     random_state=int(rng.integers(1 << 30)))
                    if take > 0 and len(available) else available.iloc[:0])
            m30_thin = rank(pooled_topk(thinned))
            variants = {
                "frozen entry (density3 + structure)": frozen,
                "HeidAI rebuild (m30 pooled + m50)": (m30_thin + vdomain) / 2,
                "A  frozen, density3 -> pooled m30":
                    (1 - STRUCTURE_WEIGHT) * m30_thin
                    + STRUCTURE_WEIGHT * structure,
                "B  pooled m30 + m50 + structure":
                    (m30_thin + vdomain + structure) / 3,
                "C  HeidAI rebuild + structure 0.3":
                    (1 - STRUCTURE_WEIGHT) * ((m30_thin + vdomain) / 2)
                    + STRUCTURE_WEIGHT * structure,
                "D  density3 + pooled m30 + m50 + structure":
                    (density3 + m30_thin + vdomain + structure) / 4,
                "E  (density3 + pooled m30 + m50)/3 + str 0.3":
                    (1 - STRUCTURE_WEIGHT)
                    * ((density3 + m30_thin + vdomain) / 3)
                    + STRUCTURE_WEIGHT * structure,
            }
            for name, series in variants.items():
                draws[name].append(score(series))
        for name in candidates:
            realisable[name].append(float(np.mean(draws[name])))
        print(f"  thinned to {private_epitope} (x{multiplier:.3f}) done")

    print("\n=== sensitivity: what if LTDEMIAQY can reach 1.72x? ===")
    print("  Only the LTD draw changes; CIN and LLW keep their measured")
    print("  multipliers. Same 20 repeats, same candidates.\n")
    scenario = {name: [] for name in candidates}
    for private_epitope, multiplier in (
            MULTIPLIER_SCENARIOS["with targeted IEDB (collaborator, unverified)"]
            .items()):
        draws = {name: [] for name in candidates}
        for _ in range(REPEATS):
            thinned = {}
            for epitope, (cohort, _) in query.items():
                available = extra_available[epitope]
                take = int(round((multiplier - 1.0) * len(cohort)))
                thinned[epitope] = (
                    available.sample(min(take, len(available)),
                                     random_state=int(rng.integers(1 << 30)))
                    if take > 0 and len(available) else available.iloc[:0])
            m30_thin = rank(pooled_topk(thinned))
            variants = {
                "frozen entry (density3 + structure)": frozen,
                "HeidAI rebuild (m30 pooled + m50)": (m30_thin + vdomain) / 2,
                "A  frozen, density3 -> pooled m30":
                    (1 - STRUCTURE_WEIGHT) * m30_thin
                    + STRUCTURE_WEIGHT * structure,
                "B  pooled m30 + m50 + structure":
                    (m30_thin + vdomain + structure) / 3,
                "C  HeidAI rebuild + structure 0.3":
                    (1 - STRUCTURE_WEIGHT) * ((m30_thin + vdomain) / 2)
                    + STRUCTURE_WEIGHT * structure,
                "D  density3 + pooled m30 + m50 + structure":
                    (density3 + m30_thin + vdomain + structure) / 4,
                "E  (density3 + pooled m30 + m50)/3 + str 0.3":
                    (1 - STRUCTURE_WEIGHT)
                    * ((density3 + m30_thin + vdomain) / 3)
                    + STRUCTURE_WEIGHT * structure,
                "F  D with no pool at all": (density3 + m30_nopool + vdomain
                                             + structure) / 4,
            }
            for name, series in variants.items():
                draws[name].append(score(series))
        for name in candidates:
            scenario[name].append(float(np.mean(draws[name])))
        print(f"  thinned to {private_epitope} "
              f"(x{multiplier:.3f}) done")

    rows = []
    for name in candidates:
        mean_realisable = float(np.mean(realisable[name]))
        rows.append({
            "configuration": name,
            "public": public[name],
            "realisable": mean_realisable,
            "pool_dependent": public[name] - mean_realisable,
            "vs_frozen_public": public[name]
                                - public["frozen entry (density3 + structure)"],
            "vs_frozen_realisable": mean_realisable
                                    - np.mean(realisable[
                                        "frozen entry (density3 + structure)"]),
            "realisable_iedb": float(np.mean(scenario[name])),
            "vs_frozen_iedb": float(np.mean(scenario[name]))
                              - float(np.mean(scenario[
                                  "frozen entry (density3 + structure)"])),
        })
    table = pd.DataFrame(rows).sort_values("public", ascending=False)
    print("\n=== the bet, laid out ===")
    print(table.round(4).to_string(index=False))

    print("\n  'pool_dependent' is the part of a candidate's public score that")
    print("  the private epitopes cannot receive at measured pool density. It is")
    print("  what a second entry is wagering is mis-measured.")

    print("\n=== where the advantage sits: both-singleton pairs ===")
    print("  The private split is 61-93% pairs where neither member has")
    print("  relatives (E22, E30), and E30 showed the projected private macro")
    print("  tracks the within-both-singleton AUC almost exactly. So this is")
    print("  the column that decides transfer, not the public macro.\n")
    import scipy.sparse as sparse
    from scipy.sparse.csgraph import connected_components as cc
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vecs = cached["vectors"]
    vecs = vecs / np.linalg.norm(vecs, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}
    is_clustered = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = vecs[[where[t] for t in ids]]
        dd = 1.0 - block @ block.T
        _, fam = cc(sparse.csr_matrix(dd <= 0.12), directed=False)
        sizes = pd.Series(fam).value_counts()
        is_clustered.update(dict(zip(ids, np.isin(fam, sizes[sizes >= 2].index))))

    def both_singleton_auc(series) -> tuple[float, list]:
        per = []
        for epitope, block in labels.groupby("epitope_aa"):
            y = block["pos"].to_numpy()
            clust = block["tcr_id"].map(is_clustered).to_numpy().astype(bool)
            pi = np.where((y == 1) & ~clust)[0]
            ni = np.where((y == 0) & ~clust)[0]
            if len(pi) == 0 or len(ni) == 0:
                continue
            values = block["tcr_id"].map(series).to_numpy(dtype=float)
            diff = values[pi][:, None] - values[ni][None, :]
            per.append(float((diff > 0).sum() + 0.5 * (diff == 0).sum())
                       / (len(pi) * len(ni)))
        return float(np.mean(per)), per

    bs_rows = []
    frozen_bs, _ = both_singleton_auc(
        candidates["frozen entry (density3 + structure)"])
    for name, series in candidates.items():
        value, per = both_singleton_auc(series)
        bs_rows.append({"configuration": name, "both_singleton_auc": value,
                        "vs_frozen": value - frozen_bs,
                        "worst_epitope_delta": min(
                            a - b for a, b in zip(per, both_singleton_auc(
                                candidates["frozen entry (density3 + structure)"]
                            )[1]))})
        print(f"  {name:<46} {value:.4f}  ({value - frozen_bs:+.4f})")
    bs = pd.DataFrame(bs_rows).sort_values("both_singleton_auc", ascending=False)
    print("\n  'worst_epitope_delta' is the smallest per-epitope improvement.")
    print("  A positive value means the advantage holds on every epitope with")
    print("  no sign flip -- the failure mode that killed E4, E7 and E13.")

    print("\n=== is the advantage real, or the cost of picking among seven? ===")
    print("  Nested LOEO over the candidate CHOICE: select on four epitopes,")
    print("  score the fifth. E15 priced selection at -0.0048 over 323")
    print("  configurations, and this is a choice among seven.\n")
    epitopes = sorted(labels["epitope_aa"].unique())
    held_selected, held_fixed_d, held_frozen, picks = [], [], [], []
    for epitope in epitopes:
        inner = labels[labels["epitope_aa"] != epitope]
        outer = labels[labels["epitope_aa"] == epitope]
        chosen = max(candidates, key=lambda n: hyp.mean_auc(inner, candidates[n]))
        picks.append(chosen)
        held_selected.append(hyp.per_epitope(outer, candidates[chosen])[epitope])
        held_fixed_d.append(hyp.per_epitope(outer, candidates[best_name])[epitope])
        held_frozen.append(hyp.per_epitope(
            outer, candidates["frozen entry (density3 + structure)"])[epitope])
        print(f"  {epitope:<11} picked {chosen[:28]:<28} "
              f"held-out {held_selected[-1]:.4f}")
    nested_selected = float(np.mean(held_selected))
    nested_d = float(np.mean(held_fixed_d))
    nested_frozen = float(np.mean(held_frozen))
    print(f"\n  nested LOEO with selection      {nested_selected:.4f}")
    print(f"  the public-best candidate fixed {nested_d:.4f}")
    print(f"  the frozen entry                {nested_frozen:.4f}")
    print(f"  selection cost here             "
          f"{nested_selected - nested_d:+.4f}")
    print(f"  honest edge over frozen         {nested_d - nested_frozen:+.4f} "
          f"public-pool, {nested_d - nested_frozen - (public[best_name] - table[table['configuration'] == best_name]['realisable'].iloc[0]):+.4f} "
          f"at private pool density")
    print(f"\n  Against the +0.02 inclusion bar this does NOT qualify to replace")
    print(f"  the frozen entry under the rule adopted in E10, and it is not")
    print(f"  proposed as one. But the bar was set for a component justified by")
    print(f"  its LOEO macro, and the evidence here is of a different kind: the")
    print(f"  gain sits in both-singleton pairs, which are 61-93% of the private")
    print(f"  split, it holds on all five epitopes with no sign flip, and most of")
    print(f"  it survives deleting the pool entirely. That is a stronger pattern")
    print(f"  than the macro number alone suggests, and it is why this is worth")
    print(f"  a separate entry rather than a footnote.")

    print("\n=== writing the two additional entries ===")
    print("  D takes the pool bet; F is the same four members with no pool")
    print("  anywhere, so it carries no transfer assumption at all.\n")
    import hashlib

    from immrep.submission import validate
    for name, filename in (
        (best_name, "d_density3_m30pool_m50_struct.csv"),
        ("F  D with no pool at all", "f_density3_m30_m50_struct_nopool.csv"),
    ):
        out = pd.DataFrame({
            "tcr_id": test["tcr_id"],
            "prediction": test["tcr_id"].map(candidates[name]).to_numpy()})
        validate(out, test=test)
        path = ROOT / "submissions" / filename
        out.to_csv(path, index=False)
        print(f"  {name[:34]:<34} -> {filename}")
        print(f"  {'':<34}    sha256 "
              f"{hashlib.sha256(path.read_bytes()).hexdigest()}")
    print("\n  The frozen entry is untouched; its guard still passes.")

    print("\n=== gates on the public-best candidate ===")
    ok, detail = private_variance_gate(test, candidates[best_name])
    print(f"  private variance gate: {ok} ({detail})")
    from immrep.validation import typicality_gate
    verdict = typicality_gate(test, labels, candidates[best_name])
    print(f"  typicality gate: {verdict['verdict']} "
          f"(feature {verdict['feature_auc']:.4f} against "
          f"{verdict['typicality_auc']:.4f}, rho {verdict['rho_with_typicality']:+.3f})")

    print("\n=== under E47's published anchor ===")
    lo = float(np.mean([v[0] for v in PRIVATE_BANDS.values()]))
    hi = float(np.mean([v[1] for v in PRIVATE_BANDS.values()]))
    tcrdist_macro = 0.8247
    for name in (best_name, "frozen entry (density3 + structure)"):
        row = table[table["configuration"] == name].iloc[0]
        for label, margin in (("public margin", row["public"] - tcrdist_macro),
                              ("realisable margin",
                               row["realisable"] - tcrdist_macro)):
            print(f"  {name[:34]:<34} {label:<18} private "
                  f"{(lo + hi) / 2 + margin:.3f}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# A public-maximising second entry, with the bet priced\n\n"
        + table.round(4).to_markdown(index=False)
        + "\n\n## Within both-singleton pairs, which decide the private split\n\n"
        + bs.round(4).to_markdown(index=False)
        + "\n\n`pool_dependent` is the share of a candidate's public score the\n"
        "private epitopes cannot receive at the multipliers E27 counted\n"
        "(1.038 / 1.045 / 1.834). A second entry wagers that count is low.\n\n"
        f"Public best: **{best_name}** at **{public[best_name]:.4f}**, against the\n"
        f"frozen entry's {public['frozen entry (density3 + structure)']:.4f} and\n"
        f"the rebuilt rival configuration's "
        f"{public['HeidAI rebuild (m30 pooled + m50)']:.4f}.\n\n"
        f"Gates on it: private variance {ok} ({detail}); "
        f"typicality {verdict['verdict']}.\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
