#!/usr/bin/env python
"""Does E30's private projection actually transfer? A leave-one-epitope-out test.

E30 revised the private estimate from LOEO's 0.846 to 0.80 by decomposing AUC
into E22's pair classes and recombining with the private epitopes' class
composition. A reviewer identified a problem with that, and the problem is real.

**The composition is not label-free.** A pair class such as `pos_clu_neg_sing`
needs to know which member is positive. Writing a = P(clustered | positive),
b = P(clustered | negative) and pi = P(positive), the shares are

    pos_clu_neg_sing = a (1 - b)        both_singleton   = (1 - a)(1 - b)
    pos_sing_neg_clu = (1 - a) b        clustered-both   = a b

while the only thing observable without labels is the overall clustered rate

    q = pi a + (1 - pi) b

One equation, two unknowns. So E30 did not measure the private composition; it
must have assumed something that maps q to (a, b), and the result is driven by
that assumption, because E20 established a >> b on the public side.

**And there is a second unobservable the review did not name: pi.** Validation
rates run 28.7% to 66.3% across the public epitopes, so pi is not pinned either.
Splitting the both-clustered mass into same- and different-cluster pairs needs a
third, s. The projection therefore has three transferred parameters, not zero.

**A third problem, found while reading the code rather than the review.** E22
clusters only the *labelled* rows (`scripts/20`, which groups `labels`), while
E29's family rates — the table E30 reasoned from — cluster the *whole* cohort
(`scripts/25`, which groups `test`). E30 combined error rates measured on one
basis with compositions measured on the other. Everything here is recomputed on
the whole-cohort basis, which is the only one available for private epitopes.

**What this script does.** It rebuilds the projection as an explicit function of
(q, pi, lambda, s) and then tests the whole pipeline out of sample: for each
public epitope, estimate every transferred quantity from the other four, predict
the held-out epitope's AUC from its label-free q alone, and compare. Three
variants separate where any error comes from:

    full      error rates and composition both transferred from four epitopes
    oracle-composition   true composition of the held-out epitope, transferred errors
    oracle-errors        true error rates of the held-out epitope, transferred composition

A projection that cannot predict a held-out public epitope has no standing to
predict a private one. Finally the private estimate is reported as a sensitivity
surface over the unobservables rather than as a point.

Run: python scripts/37_private_projection_validation.py
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
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/14_private_projection_validation.md"
CLUSTER_RADIUS = 0.12
CLASSES = ("pos_clu_neg_sing", "pos_sing_neg_clu", "both_singleton",
           "same_cluster", "diff_cluster")


def submitted_score(test: pd.DataFrame) -> pd.Series:
    """The scoring function actually submitted, reused from scripts/20."""
    spec = importlib.util.spec_from_file_location(
        "decomp", ROOT / "scripts/20_auc_error_decomposition.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["decomp"] = module
    spec.loader.exec_module(module)
    return module.submitted_score(test)


def cluster_cohort(ids: np.ndarray, normed: np.ndarray,
                   position: dict) -> tuple[np.ndarray, np.ndarray]:
    """Single-linkage clustering over the WHOLE cohort, which is the only basis
    a private epitope can supply."""
    block = normed[[position[t] for t in ids]]
    distances = 1.0 - block @ block.T
    _, cluster = connected_components(
        sp.csr_matrix(distances <= CLUSTER_RADIUS), directed=False)
    sizes = pd.Series(cluster).value_counts()
    clustered = np.isin(cluster, sizes[sizes >= 2].index)
    return cluster, clustered


def observed(test, labels, scores, normed, position) -> tuple[pd.DataFrame, dict]:
    """Per-epitope truth: parameters, class shares, class error rates, AUC."""
    rows, errors = [], {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        cluster, clustered = cluster_cohort(ids, normed, position)
        order = {t: i for i, t in enumerate(ids)}
        record = {"epitope": epitope,
                  "split": cohort["split"].iloc[0],
                  "q_cohort": float(clustered.mean())}

        subset = labels[labels["epitope_aa"] == epitope]
        if subset.empty:
            rows.append(record)
            continue

        pick = np.array([order[t] for t in subset["tcr_id"]])
        y = subset["pos"].to_numpy()
        s = subset["tcr_id"].map(scores).to_numpy()
        c = clustered[pick]
        cid = cluster[pick]

        pos_idx, neg_idx = np.where(y == 1)[0], np.where(y == 0)[0]
        pc, nc = c[pos_idx][:, None], c[neg_idx][None, :]
        same = (cid[pos_idx][:, None] == cid[neg_idx][None, :]) & pc & nc
        klass = np.full((len(pos_idx), len(neg_idx)), "", dtype=object)
        klass[pc & nc & same] = "same_cluster"
        klass[pc & nc & ~same] = "diff_cluster"
        klass[pc & ~nc] = "pos_clu_neg_sing"
        klass[~pc & nc] = "pos_sing_neg_clu"
        klass[~pc & ~nc] = "both_singleton"

        sp_s, sn_s = s[pos_idx][:, None], s[neg_idx][None, :]
        wrong = (sp_s < sn_s).astype(float) + 0.5 * (sp_s == sn_s)
        total = wrong.size

        a = float(c[pos_idx].mean())
        b = float(c[neg_idx].mean())
        both = (pc & nc)
        record.update({
            "auc": float(roc_auc_score(y, s)),
            "pi": float(y.mean()), "a": a, "b": b,
            "lam": a / b if b > 0 else np.nan,
            "s_same": float(same.sum() / both.sum()) if both.sum() else 0.0,
            "q_labelled": float(c.mean()), "n_pairs": int(total),
        })
        per_class = {}
        for name in CLASSES:
            mask = klass == name
            record[f"share_{name}"] = mask.sum() / total
            per_class[name] = (float(wrong[mask].sum()), int(mask.sum()))
            record[f"err_{name}"] = (wrong[mask].mean() if mask.any() else np.nan)
        errors[epitope] = per_class
        rows.append(record)
    return pd.DataFrame(rows), errors


def pooled_errors(errors: dict, exclude: str | None) -> dict:
    """Class error rates pooled over epitopes, weighted by pair count."""
    out = {}
    for name in CLASSES:
        num = den = 0.0
        for epitope, per_class in errors.items():
            if epitope == exclude:
                continue
            w, n = per_class[name]
            num += w
            den += n
        out[name] = num / den if den else np.nan
    return out


def compose(q: float, pi: float, lam: float, s_same: float) -> dict:
    """Class shares implied by the observable q and the three transferred
    parameters. Solves q = pi*a + (1-pi)*b under a = lam*b."""
    denominator = pi * lam + (1.0 - pi)
    b = q / denominator
    a = lam * b
    a, b = min(a, 1.0), min(b, 1.0)
    return {
        "pos_clu_neg_sing": a * (1 - b),
        "pos_sing_neg_clu": (1 - a) * b,
        "both_singleton": (1 - a) * (1 - b),
        "same_cluster": a * b * s_same,
        "diff_cluster": a * b * (1 - s_same),
    }


def predict(shares: dict, errs: dict) -> float:
    return 1.0 - sum(shares[c] * errs[c] for c in CLASSES)


def main() -> int:
    test = load_test()
    labels = build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    print("scoring with the submitted model ...")
    scores = submitted_score(test)

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    normed = cached["vectors"] / np.linalg.norm(
        cached["vectors"], axis=1, keepdims=True)

    frame, errors = observed(test, labels, scores, normed, position)
    public = frame[frame["split"] == "public"].copy()
    private = frame[frame["split"] == "private"].copy()

    print("\n=== the identifiability problem, in the public data ===")
    print("  q is observable without labels; a and b are not, and nor is pi.")
    print(public[["epitope", "q_cohort", "q_labelled", "pi", "a", "b", "lam",
                  "s_same"]].round(3).to_string(index=False))
    print(f"\n  lambda = a/b spans {public['lam'].min():.2f} to "
          f"{public['lam'].max():.2f}; pi spans {public['pi'].min():.3f} to "
          f"{public['pi'].max():.3f}.")
    print("  Neither is observable on a private epitope, so the composition")
    print("  there is assumed, not measured.")

    print("\n=== the reconstruction E30 cited is an identity, not a test ===")
    identity = []
    for row in public.to_dict("records"):
        rebuilt = 1.0 - sum(row[f"share_{c}"] * row[f"err_{c}"] for c in CLASSES)
        identity.append({"epitope": row["epitope"], "actual": row["auc"],
                         "rebuilt": rebuilt, "residual": rebuilt - row["auc"]})
    identity = pd.DataFrame(identity)
    print(identity.round(5).to_string(index=False))
    print(f"  max |residual| {identity['residual'].abs().max():.2e} -- this uses each")
    print("  epitope's own error rates and own shares, so it can only return that")
    print("  epitope's own AUC. It checks the arithmetic, not transferability.")

    print("\n=== P0: leave-one-epitope-out transfer test ===")
    print("  Everything transferred is estimated on the other four epitopes;")
    print("  the held-out epitope contributes only its label-free q.\n")
    results = []
    for row in public.to_dict("records"):
        epitope = row["epitope"]
        others = public[public["epitope"] != epitope]
        errs = pooled_errors(errors, exclude=epitope)
        lam = float(others["lam"].median())
        pi = float(others["pi"].median())
        s_same = float(others["s_same"].median())

        transferred = compose(row["q_cohort"], pi, lam, s_same)
        truth = {c: row[f"share_{c}"] for c in CLASSES}
        own_errs = {c: row[f"err_{c}"] for c in CLASSES}

        results.append({
            "epitope": epitope,
            "actual": row["auc"],
            "full": predict(transferred, errs),
            "oracle_composition": predict(truth, errs),
            "oracle_errors": predict(transferred, own_errs),
        })
    results = pd.DataFrame(results)
    for column in ("full", "oracle_composition", "oracle_errors"):
        results[f"d_{column}"] = results[column] - results["actual"]
    print(results.round(4).to_string(index=False))

    summary = {}
    for column in ("full", "oracle_composition", "oracle_errors"):
        d = results[f"d_{column}"]
        summary[column] = {"mean_error": d.mean(), "mae": d.abs().mean(),
                           "max_abs": d.abs().max()}
    summary = pd.DataFrame(summary).T
    print("\n  prediction error of each variant (predicted minus actual):")
    print(summary.round(4).to_string())
    print("\n  If `full` has a large MAE, the projection cannot predict a held-out")
    print("  public epitope and therefore cannot be quoted as a private point.")
    print("  Comparing the two oracles says whether the damage is in the error")
    print("  rates or in the composition.")

    print("\n=== private projection as a sensitivity surface ===")
    errs_all = pooled_errors(errors, exclude=None)
    lam_grid = [public["lam"].min(), public["lam"].median(), public["lam"].max()]
    pi_grid = [0.30, 0.40, 0.50, 0.60, 0.70]
    s_same = float(public["s_same"].median())
    surface = []
    for lam in lam_grid:
        for pi in pi_grid:
            macro = np.mean([
                predict(compose(r["q_cohort"], pi, lam, s_same), errs_all)
                for r in private.to_dict("records")])
            surface.append({"lambda": lam, "pi": pi, "private_macro": macro})
    surface = pd.DataFrame(surface)
    grid = surface.pivot(index="lambda", columns="pi", values="private_macro")
    print(grid.round(3).to_string())
    print(f"\n  spread across the grid: {surface['private_macro'].min():.3f} to "
          f"{surface['private_macro'].max():.3f}")

    # The honest interval has to carry the out-of-sample error too.
    bias = summary.loc["full", "mean_error"]
    mae = summary.loc["full", "mae"]
    centre = float(surface["private_macro"].median())
    print(f"\n  grid median {centre:.3f}; LOEO transfer bias {bias:+.4f}, "
          f"MAE {mae:.4f}")
    print(f"  bias-corrected centre {centre - bias:.3f}, and any interval must")
    print(f"  also carry the +/-{mae:.3f} the method misses by on held-out public")
    print("  epitopes, on top of the grid spread above.")

    print("\n=== correcting a factual claim in E30 ===")
    lowest_public = public.loc[public["q_cohort"].idxmin()]
    highest_private = private.loc[private["q_cohort"].idxmax()]
    print(f"  lowest public clustered rate:  {lowest_public['epitope']} "
          f"{lowest_public['q_cohort']:.3f}")
    print(f"  highest private clustered rate: {highest_private['epitope']} "
          f"{highest_private['q_cohort']:.3f}")
    print("  E30 said all five public epitopes are family-richer than the private")
    print("  three. That is false: "
          f"{highest_private['epitope']} sits above {lowest_public['epitope']}.")
    print("  Two of three private epitopes are below every public one; the third")
    print("  is inside the public range.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Does E30's private projection transfer?\n\n"
        "## Observable and unobservable quantities\n\n"
        + public[["epitope", "q_cohort", "pi", "a", "b", "lam", "s_same"]]
        .round(3).to_markdown(index=False)
        + "\n\n`q` is label-free. `pi`, `a`, `b` are not, and the class "
          "composition needs `a` and `b` separately.\n\n"
          "## Leave-one-epitope-out transfer test\n\n"
        + results.round(4).to_markdown(index=False)
        + "\n\n" + summary.round(4).to_markdown()
        + "\n\n## Private sensitivity surface\n\n"
        + grid.round(3).to_markdown() + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
