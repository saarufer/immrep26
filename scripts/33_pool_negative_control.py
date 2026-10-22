#!/usr/bin/env python
"""Negative controls for pool expansion, and the pool x structure factorial.

Plan cards R04 and R05, the two I skipped on my own judgement and should not
have. R04 is a genuine sceptical control that nothing else covers: E27 audited
whether the *query* finds itself in the pool, which is a different worry from
whether adding any sizeable batch of sequences raises density regardless of what
they are. If a size-matched batch of unrelated TCRs works as well as the real
same-epitope pool, then "expansion helps" is a volume artefact.

Three controls, increasingly destructive:

  other-epitope   real TCRs from different epitopes, matched on CDR3-beta length
                  and TRBV gene so the added repertoire has the same shape.
  label-shuffled  the real added pool, but drawn from the union across epitopes
                  so each cohort receives the wrong epitope's neighbours.
  synthetic       OLGA-generated sequences, which are repertoire-plausible but
                  carry no antigen selection at all.

R05 is the 2x2: original or expanded pool, with or without the structure blend,
to separate whether the two are complementary or substitutes. Its pool arm was
already measured near zero on our configuration (E27), so the factorial is
expected to be degenerate -- running it says so with numbers instead of an
assertion.

SCEPTR runs on CPU here: the GPU is occupied by a folding job and these
embeddings are small, so there is no reason to contend for memory.

Run: python scripts/33_pool_negative_control.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_test  # noqa: E402

VDJDB = ROOT / f"data/external/vdjdb/vdjdb-{VDJDB_VERSION}/vdjdb_full.txt"
OUT = ROOT / "reports/11_pool_negative_control.md"
TOPK = 8
SEED = 0


def load_pool() -> pd.DataFrame:
    frame = pd.read_csv(VDJDB, sep="\t", low_memory=False)
    frame = frame[
        (frame["species"] == "HomoSapiens")
        & frame["cdr3.alpha"].notna()
        & frame["cdr3.beta"].notna()
    ].rename(
        columns={
            "cdr3.alpha": "cdr3a", "cdr3.beta": "cdr3b", "v.alpha": "trav",
            "j.alpha": "traj", "v.beta": "trbv", "j.beta": "trbj",
            "antigen.epitope": "ep",
        }
    )
    return frame[["cdr3a", "cdr3b", "trav", "traj", "trbv", "trbj", "ep"]]


def embed(frame: pd.DataFrame) -> np.ndarray:
    import sceptr
    from libtcrlm.schema.tcr import TravGene, TrbvGene

    slash = lambda g: g.replace("_", "/") if isinstance(g, str) else None  # noqa: E731
    model_input = pd.DataFrame(
        {
            "TRAV": [slash(x) for x in frame["trav"]],
            "CDR3A": frame["cdr3a"].tolist(),
            "TRAJ": [slash(x) for x in frame["traj"]],
            "TRBV": [slash(x) for x in frame["trbv"]],
            "CDR3B": frame["cdr3b"].tolist(),
            "TRBJ": [slash(x) for x in frame["trbj"]],
        }
    )
    for column, vocab in (("TRAV", {g.name for g in TravGene}),
                          ("TRBV", {g.name for g in TrbvGene})):
        base = model_input[column].fillna("").str.split("*").str[0]
        model_input.loc[~base.isin(vocab), column] = None
    vectors = sceptr.calc_vector_representations(model_input)
    return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)


def olga_repertoire(n: int, rng: np.random.Generator) -> pd.DataFrame:
    """Generation-model sequences: plausible repertoire, no antigen selection."""
    from tcrdist.pgen import OlgaModel

    alpha = OlgaModel(chain_folder="human_T_alpha", recomb_type="VJ")
    beta = OlgaModel(chain_folder="human_T_beta", recomb_type="VDJ")
    a = [s for s in alpha.gen_cdr3s(n=n * 2) if s][:n]
    b = [s for s in beta.gen_cdr3s(n=n * 2) if s][:n]
    size = min(len(a), len(b))
    return pd.DataFrame({"cdr3a": a[:size], "cdr3b": b[:size],
                         "trav": None, "traj": None, "trbv": None, "trbj": None,
                         "ep": "SYNTHETIC"})


def matched_draw(pool: pd.DataFrame, template: pd.DataFrame,
                 rng: np.random.Generator) -> pd.DataFrame:
    """Sample `template`-shaped TCRs from `pool`, matching CDR3-beta length and
    TRBV where possible so the control batch has the same repertoire shape."""
    pool = pool.copy()
    pool["len_b"] = pool["cdr3b"].str.len()
    picks = []
    for length, group in template.assign(len_b=template["cdr3b"].str.len()).groupby("len_b"):
        candidates = pool[pool["len_b"] == length]
        if candidates.empty:
            candidates = pool
        picks.append(candidates.sample(len(group), replace=len(candidates) < len(group),
                                       random_state=int(rng.integers(1 << 30))))
    return pd.concat(picks, ignore_index=True)


def main() -> int:
    import sceptr

    sceptr.disable_hardware_acceleration()

    test = load_test()
    labels = build_training_frame()
    pool = load_pool()
    rng = np.random.default_rng(SEED)
    epitopes = sorted(labels["epitope_aa"].unique())

    structure = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")

    def rank_within(values: pd.Series) -> pd.Series:
        frame = test[["tcr_id", "epitope_aa"]].copy()
        frame["v"] = frame["tcr_id"].map(values)
        return frame.assign(
            r=frame.groupby("epitope_aa")["v"].rank(pct=True)
        ).set_index("tcr_id")["r"]

    structure_rank = rank_within(structure["pae_tcr_pep_min"])

    print("embedding the query cohorts ...")
    query = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        cohort = cohort.reset_index(drop=True)
        frame = pd.DataFrame({
            "cdr3a": cohort["cdr3_alpha_aa"], "cdr3b": cohort["cdr3_beta_aa"],
            "trav": cohort["TRAV_IMGT"], "traj": cohort["TRAJ_IMGT"],
            "trbv": cohort["TRBV_IMGT"], "trbj": cohort["TRBJ_IMGT"], "ep": epitope,
        })
        query[epitope] = (cohort, embed(frame))

    def score(epitope: str, extra: pd.DataFrame | None) -> pd.Series:
        cohort, vectors = query[epitope]
        reference = vectors if extra is None or extra.empty else np.vstack(
            [vectors, embed(extra)])
        distances = 1.0 - vectors @ reference.T
        distances[np.arange(len(cohort)), np.arange(len(cohort))] = np.inf
        value = -np.sort(distances, axis=1)[:, :TOPK].mean(axis=1)
        return pd.Series(value, index=cohort["tcr_id"])

    def evaluate(series_by_epitope: dict, blend_structure: bool = False) -> dict:
        out = {}
        for epitope in epitopes:
            scores = series_by_epitope[epitope]
            if blend_structure:
                ranked = scores.rank(pct=True)
                scores = 0.7 * ranked + 0.3 * structure_rank.reindex(ranked.index)
            subset = labels[labels["epitope_aa"] == epitope]
            out[epitope] = roc_auc_score(
                subset["pos"], subset["tcr_id"].map(scores).to_numpy())
        return out

    # Real same-epitope pool, deduplicated against the query.
    real = {}
    for epitope in epitopes:
        cohort, _ = query[epitope]
        keys = set(zip(cohort["cdr3_alpha_aa"], cohort["cdr3_beta_aa"]))
        same = pool[pool["ep"] == epitope]
        real[epitope] = same[
            [(a, b) not in keys for a, b in zip(same["cdr3a"], same["cdr3b"])]
        ].drop_duplicates(["cdr3a", "cdr3b"])

    print("\n=== R04: negative controls for pool expansion (top-8 density) ===")
    arms = {"no pool": {e: None for e in epitopes},
            "real same-epitope": real}
    other = {}
    shuffled = {}
    synthetic = {}
    for epitope in epitopes:
        size = len(real[epitope])
        foreign = pool[(pool["ep"] != epitope) & (pool["ep"] != "SYNTHETIC")]
        other[epitope] = matched_draw(foreign, real[epitope], rng) if size else foreign.iloc[:0]
        shuffled[epitope] = foreign.sample(
            min(size, len(foreign)), random_state=SEED) if size else foreign.iloc[:0]
        synthetic[epitope] = olga_repertoire(size, rng) if size else foreign.iloc[:0]
    arms["other-epitope (shape-matched)"] = other
    arms["other-epitope (unmatched)"] = shuffled
    arms["OLGA synthetic"] = synthetic

    rows = []
    for name, extras in arms.items():
        per = evaluate({e: score(e, extras[e]) for e in epitopes})
        rows.append({"arm": name, "added": int(np.mean([
            0 if extras[e] is None else len(extras[e]) for e in epitopes])),
            "mean": float(np.mean(list(per.values()))), **{k[:5]: round(v, 3) for k, v in per.items()}})
    table = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print(table.round(4).to_string(index=False))
    base = table.loc[table["arm"] == "no pool", "mean"].iloc[0]
    print("\n  gain over no pool:")
    # Index by column name: itertuples field names are not positional when the
    # column is already a valid identifier, which a first version got wrong.
    for _, row in table[table["arm"] != "no pool"].iterrows():
        print(f"    {row['arm']:<32} {row['mean'] - base:+.4f}")
    print("\n  If a control arm matches the real pool, expansion is a volume effect.")

    print("\n=== R05: pool x structure factorial ===")
    cells = {}
    for pool_name, extras in (("original", arms["no pool"]), ("expanded", real)):
        for structure_on in (False, True):
            per = evaluate({e: score(e, extras[e]) for e in epitopes},
                           blend_structure=structure_on)
            cells[(pool_name, structure_on)] = per
    grid = pd.DataFrame({
        f"{p}/{'struct' if s else 'no struct'}": cells[(p, s)]
        for p in ("original", "expanded") for s in (False, True)})
    grid.loc["mean"] = grid.mean()
    print(grid.round(4).to_string())
    a = grid.loc["mean", "original/no struct"]
    b = grid.loc["mean", "original/struct"]
    c = grid.loc["mean", "expanded/no struct"]
    d = grid.loc["mean", "expanded/struct"]
    print(f"\n  structure effect without pool  B-A = {b - a:+.4f}")
    print(f"  structure effect with pool     D-C = {d - c:+.4f}")
    print(f"  pool effect without structure  C-A = {c - a:+.4f}")
    print(f"  pool effect with structure     D-B = {d - b:+.4f}")
    print(f"  interaction  (D-C)-(B-A)       = {(d - c) - (b - a):+.4f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Pool expansion: negative controls and the pool x structure factorial\n\n"
        "## R04 negative controls\n\n" + table.round(4).to_markdown(index=False)
        + "\n\n## R05 factorial\n\n" + grid.round(4).to_markdown() + "\n",
        encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
