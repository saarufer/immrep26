#!/usr/bin/env python
"""Four new hypotheses from the PI, tested under the existing discipline.

The predictor was frozen on 2026-10-25. That was my call, not the brief: the
objective is to keep trying until the submission deadline, and freezing is only
right when no hypothesis is left worth testing. Four were.

H3 is already closed and is not in this script. It asked whether convergent
recombination can be told from duplicate sampling of one clone. Its own stop
condition settles it: VDJdb carries **no nucleotide column at all**, and
`meta.subject.id` is 0% non-null on two of the three private epitopes
(LLWNGPMAV, LTDEMIAQY) with `meta.clone.id` 0% on LTDEMIAQY. On the epitopes
that decide the ranking there is nothing to distinguish independent rearrangement
from repeat curation, and record count must not stand in for it.

The three tested here, with why each is not already excluded:

H1  Biophysical substitution distance. E32's card D04 compared identity (0.8007),
    BLOSUM62 (0.7999) and BLOSUM45 (0.8014) and concluded the cost matrix does
    not matter -- but all three are evolutionary or trivial. A physicochemical
    distance is a different axis. Grantham (1974) is used rather than a matrix
    invented here: composition, polarity and side-chain volume, with the
    published coefficients.

H2  Alpha-chain partner diversity. E6 found alpha outranks beta (0.799 vs
    0.751) and E29's D03 tested alpha/beta rank *discordance*, but neither asked
    whether one alpha tolerates many sequence-diverse beta partners. That is a
    bipartite-graph property, not a distance. A beta-centric mirror is computed
    as a control: if beta scores the same, nothing is alpha-specific.

H4  Expression and assembly liabilities. Nothing has tested this. The assay
    re-expresses historical TCRs, so some negatives may fail at expression,
    pairing or folding rather than at recognition. Pre-declared indicators only:
    extra cysteines, N-glycosylation sequons, length deviation, non-canonical
    CDR3 termini, hydrophobic fraction.

Discipline, unchanged: every feature passes E3's viability gate (it must vary
within all three private epitopes) before its LOEO score is believed, fusion is
a zero-parameter rank average, and a component enters the submission only if it
adds more than +0.02 under strict nested LOEO.

Run: python scripts/40_new_hypotheses.py
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

REPORT = ROOT / "reports/17_new_hypotheses.md"
SIGMA = 0.1
THRESHOLD = 0.02

# Grantham 1974: composition, polarity, volume, and the published coefficients.
GRANTHAM = {
    "S": (1.42, 9.2, 32), "R": (0.65, 10.5, 124), "L": (0.00, 4.9, 111),
    "P": (0.39, 8.0, 32.5), "T": (0.71, 8.6, 61), "A": (0.00, 8.1, 31),
    "V": (0.00, 5.9, 84), "G": (0.74, 9.0, 3), "I": (0.00, 5.2, 111),
    "F": (0.00, 5.2, 132), "Y": (0.20, 6.2, 136), "C": (2.75, 5.5, 55),
    "H": (0.58, 10.4, 96), "Q": (0.89, 10.5, 85), "N": (1.33, 11.6, 56),
    "K": (0.33, 11.3, 119), "D": (1.38, 13.0, 54), "E": (0.92, 12.3, 83),
    "M": (0.00, 5.7, 105), "W": (0.13, 5.4, 170),
}
ALPHA, BETA, GAMMA = 1.833, 0.1018, 0.000399
HYDROPHOBIC = set("AILMFWVY")


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def grantham_matrix():
    """A parasail substitution matrix from Grantham distance, as similarity."""
    import parasail

    letters = "ACDEFGHIKLMNPQRSTVWY"
    raw = np.zeros((20, 20))
    for i, a in enumerate(letters):
        for j, b in enumerate(letters):
            ca, pa, va = GRANTHAM[a]
            cb, pb, vb = GRANTHAM[b]
            raw[i, j] = np.sqrt(ALPHA * (ca - cb) ** 2 + BETA * (pa - pb) ** 2
                                + GAMMA * (va - vb) ** 2)
    # Alignment wants similarity, so invert and scale to integers comparable in
    # magnitude to BLOSUM62 (identity ~ +11, most dissimilar ~ -4).
    similarity = np.rint(11 - 15 * raw / raw.max()).astype(int)
    matrix = parasail.matrix_create(letters, 1, 0)
    for i in range(20):
        for j in range(20):
            matrix.set_value(i, j, int(similarity[i, j]))
    return matrix


# --------------------------------------------------------------------------
# shared scoring


def rank_within(values: pd.Series, test: pd.DataFrame) -> pd.Series:
    frame = test[["tcr_id", "epitope_aa"]].copy()
    frame["v"] = frame["tcr_id"].map(values)
    return frame.assign(
        r=frame.groupby("epitope_aa")["v"].rank(pct=True)).set_index("tcr_id")["r"]


def per_epitope(labels: pd.DataFrame, series: pd.Series) -> dict:
    out = {}
    for epitope, block in labels.groupby("epitope_aa"):
        values = block["tcr_id"].map(series).to_numpy()
        out[epitope] = (float(roc_auc_score(block["pos"], values))
                        if len(np.unique(values)) > 1 else np.nan)
    return out


def mean_auc(labels, series) -> float:
    return float(np.nanmean(list(per_epitope(labels, series).values())))


def viability(test: pd.DataFrame, series: pd.Series) -> tuple[bool, str]:
    """E3's rule: a feature that is constant inside a private epitope can score
    well in cross-validation and contribute exactly nothing where it counts."""
    detail = []
    ok = True
    for epitope, cohort in test[test["split"] == "private"].groupby("epitope_aa"):
        values = cohort["tcr_id"].map(series)
        n = values.nunique()
        detail.append(f"{epitope[:5]}:{n}")
        if n <= 1:
            ok = False
    return ok, " ".join(detail)


def nested_blend_gain(test, labels, base_rank: pd.Series,
                      feature: pd.Series) -> float:
    """Strict nested LOEO: the weight is chosen on four epitopes and applied to
    the fifth, so the cost of choosing it is inside the estimate (E15)."""
    weights = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
    ranked = rank_within(feature, test)
    epitopes = sorted(labels["epitope_aa"].unique())
    held = []
    for epitope in epitopes:
        inner = labels[labels["epitope_aa"] != epitope]
        best, best_score = 0.0, -np.inf
        for w in weights:
            blend = (1 - w) * base_rank + w * ranked
            score = mean_auc(inner, blend)
            if score > best_score:
                best, best_score = w, score
        blend = (1 - best) * base_rank + best * ranked
        outer = labels[labels["epitope_aa"] == epitope]
        held.append(roc_auc_score(outer["pos"], outer["tcr_id"].map(blend)))
    return float(np.mean(held))


# --------------------------------------------------------------------------
# H1


def h1_biophysical(test, labels, cards) -> pd.DataFrame:
    import parasail

    print("\n=== H1: biophysical substitution distance (Grantham) ===")
    print("  D04 compared identity 0.8007, BLOSUM62 0.7999, BLOSUM45 0.8014 and")
    print("  concluded the matrix does not matter. All three are evolutionary or")
    print("  trivial; Grantham is composition, polarity and volume.\n")

    schemes = {
        "Grantham (physicochemical)": (grantham_matrix(), 11, 1),
        "BLOSUM62 (D04 reference)": (parasail.blosum62, 11, 1),
        "identity (D04 reference)": (
            parasail.matrix_create("ACDEFGHIKLMNPQRSTVWY", 1, 0), 2, 1),
    }
    rows = []
    for name, (matrix, go, ge) in schemes.items():
        values = pd.Series(index=test["tcr_id"], dtype=float)
        for epitope, cohort in test.groupby("epitope_aa"):
            cohort = cohort.reset_index(drop=True)
            distances = cards.cdr3_distance(cohort, matrix, go, ge).astype(float)
            np.fill_diagonal(distances, np.inf)
            values.loc[cohort["tcr_id"]] = np.exp(-distances / SIGMA).mean(axis=1)
        rows.append({"feature": name, "mean_auc": mean_auc(labels, values),
                     "series": values})
        print(f"  {name:<28} {rows[-1]['mean_auc']:.4f}")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# H2


def h2_partner_diversity(test, labels) -> pd.DataFrame:
    print("\n=== H2: alpha-chain partner diversity ===")
    print("  For each TCR, take its alpha-chain neighbourhood and measure how")
    print("  sequence-diverse those neighbours' beta chains are. A beta-centric")
    print("  mirror is the control: if it scores the same, nothing is")
    print("  alpha-specific.\n")

    cache = {c: np.load(ROOT / f"data/interim/sceptr/{c}.npz", allow_pickle=True)
             for c in ("a_sceptr", "b_sceptr")}
    rows = []
    for lead, partner, label in (("a_sceptr", "b_sceptr", "alpha-led"),
                                 ("b_sceptr", "a_sceptr", "beta-led (control)")):
        for radius in (0.10, 0.20):
            values = pd.Series(index=test["tcr_id"], dtype=float)
            counts = pd.Series(index=test["tcr_id"], dtype=float)
            for epitope, cohort in test.groupby("epitope_aa"):
                ids = cache[lead][f"ID__{epitope}"]
                order = {t: i for i, t in enumerate(ids)}
                pick = [order[t] for t in cohort["tcr_id"]]
                dl = cache[lead][f"D__{epitope}"][np.ix_(pick, pick)]
                dp = cache[partner][f"D__{epitope}"][np.ix_(pick, pick)]
                near = dl <= radius
                np.fill_diagonal(near, False)
                diversity = np.zeros(len(pick))
                for i in range(len(pick)):
                    idx = np.where(near[i])[0]
                    if len(idx) >= 2:
                        sub = dp[np.ix_(idx, idx)]
                        diversity[i] = sub[np.triu_indices(len(idx), 1)].mean()
                values.loc[cohort["tcr_id"]] = diversity
                counts.loc[cohort["tcr_id"]] = near.sum(axis=1)
            name = f"{label} partner diversity r={radius}"
            rows.append({"feature": name, "mean_auc": mean_auc(labels, values),
                         "series": values})
            print(f"  {name:<36} {rows[-1]['mean_auc']:.4f}"
                  f"   (median neighbours {counts.median():.0f})")
            if label.startswith("alpha"):
                rows.append({"feature": f"alpha neighbour count r={radius}",
                             "mean_auc": mean_auc(labels, counts),
                             "series": counts})
                print(f"  {'  ^ neighbour count alone':<36} "
                      f"{rows[-1]['mean_auc']:.4f}   <- the trivial confound")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# H4


def h4_liabilities(test, labels) -> pd.DataFrame:
    print("\n=== H4: expression and assembly liabilities ===")
    print("  The assay re-expresses historical TCRs, so a negative may fail at")
    print("  expression, pairing or folding rather than at recognition.")
    print("  Pre-declared indicators only, no fitted classifier.\n")

    def sequons(seq: str) -> int:
        return sum(1 for i in range(len(seq) - 2)
                   if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST")

    modal = {c: test[c].str.len().mode().iloc[0]
             for c in ("cdr3_alpha_aa", "cdr3_beta_aa")}
    frame = pd.DataFrame(index=test["tcr_id"])
    a, b = test["cdr3_alpha_aa"], test["cdr3_beta_aa"]
    da, db = test["full_seq_alpha_aa"], test["full_seq_beta_aa"]

    frame["extra_cys"] = (a.str[1:].str.count("C").values
                          + b.str[1:].str.count("C").values)
    frame["nglyc_cdr3"] = [sequons(x) + sequons(y) for x, y in zip(a, b)]
    frame["nglyc_domain"] = [sequons(x) + sequons(y) for x, y in zip(da, db)]
    frame["len_deviation"] = (
        (a.str.len() - modal["cdr3_alpha_aa"]).abs().values
        + (b.str.len() - modal["cdr3_beta_aa"]).abs().values)
    frame["noncanonical_ends"] = (
        (~a.str[-1].isin(list("FW"))).astype(int).values
        + (~b.str[-1].isin(list("FW"))).astype(int).values
        + (~a.str[0].eq("C")).astype(int).values
        + (~b.str[0].eq("C")).astype(int).values)
    frame["hydrophobic_frac"] = [
        (sum(c in HYDROPHOBIC for c in x) + sum(c in HYDROPHOBIC for c in y))
        / (len(x) + len(y)) for x, y in zip(a, b)]
    frame["liability_sum"] = (
        frame["extra_cys"] + frame["nglyc_cdr3"] + frame["noncanonical_ends"]
        + (frame["len_deviation"] > 4).astype(int))

    rows = []
    for column in frame.columns:
        # A liability should make validation *less* likely, so the predictor is
        # the negated indicator; report both so a sign flip is visible.
        auc = mean_auc(labels, -frame[column])
        rows.append({"feature": f"H4 {column} (negated)", "mean_auc": auc,
                     "series": -frame[column]})
        print(f"  {column:<20} AUC {auc:.4f}   "
              f"prevalence {np.mean(frame[column] > 0):.1%}")
    return pd.DataFrame(rows)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    cards = load_module("cards", "scripts/27_remaining_sequence_cards.py")
    decomp = load_module("decomp", "scripts/20_auc_error_decomposition.py")

    print("scoring the submitted blend as the baseline ...")
    base = decomp.submitted_score(test)
    base_rank = rank_within(base, test)
    baseline = mean_auc(labels, base_rank)
    print(f"  submitted blend: {baseline:.4f}")

    results = pd.concat([h1_biophysical(test, labels, cards),
                         h2_partner_diversity(test, labels),
                         h4_liabilities(test, labels)], ignore_index=True)

    print("\n=== E3 viability gate, then nested-LOEO blend gain ===")
    print("  A feature constant inside a private epitope is discarded before")
    print("  its score is believed, however well it does publicly.\n")
    rows = []
    for record in results.to_dict("records"):
        ok, detail = viability(test, record["series"])
        gain = np.nan
        if ok:
            gain = nested_blend_gain(test, labels, base_rank,
                                     record["series"]) - baseline
        rho = np.nanmean([
            spearmanr(cohort["tcr_id"].map(record["series"]),
                      cohort["tcr_id"].map(base_rank)).statistic
            for _, cohort in test.groupby("epitope_aa")])
        rows.append({"feature": record["feature"], "auc": record["mean_auc"],
                     "private_levels": detail, "viable": ok,
                     "rho_with_blend": rho, "nested_gain": gain})
        flag = "" if ok else "  DISCARDED: constant on a private epitope"
        print(f"  {record['feature']:<36} auc {record['mean_auc']:.4f}  "
              f"rho {rho:+.3f}  gain {gain:+.4f}{flag}"
              if ok else
              f"  {record['feature']:<36} auc {record['mean_auc']:.4f}{flag}")

    table = pd.DataFrame(rows).sort_values("nested_gain", ascending=False)
    best = table["nested_gain"].max()
    print(f"\n  best nested gain {best:+.4f} against a threshold of +{THRESHOLD}")
    print("  " + ("A component clears the bar -- revisit the frozen submission."
                  if best > THRESHOLD else
                  "Nothing clears the bar. The frozen submission stands."))

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Four new hypotheses\n\n"
        f"Submitted blend baseline {baseline:.4f}; inclusion threshold "
        f"+{THRESHOLD} under strict nested LOEO.\n\n"
        "H3 is closed by its own stop condition: VDJdb has no nucleotide column, "
        "and subject identifiers are absent on two of three private epitopes.\n\n"
        + table.round(4).to_markdown(index=False) + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
