#!/usr/bin/env python
"""E83 -- the contact fingerprint as a density space, not as a confidence score.

Every structural feature tried so far (E14, E19, E28, E36) is a **scalar
summarising one interface**, and E36 settled that those carry no
peptide-specific information: folding a TCR against anchor-matched decoys moved
ipTM by +0.0061 against a decoy-choice noise of 0.036.

This asks a different question. The 60-dimensional contact fingerprint from
`scripts/24` is indexed **by peptide position** — per position, heavy-atom
contact counts and minimum distances from TCR alpha and beta plus the CDR3
share. Within an epitope the peptide is identical, so position is a shared
coordinate and two TCRs can be compared by *where they touch the antigen*.

So: not "how confident is this interface" but "**do TCRs that engage the
peptide the same way validate the same way**". E36's negative does not cover
it, because that result is about the TCR-versus-peptide axis while this is
about the TCR-versus-TCR axis in a structural coordinate.

D already carries `pae_tcr_pep_min` at weight 0.3. The test is whether a
structural *neighbourhood* adds anything to, or beats, that structural
*scalar*.

Run: .venv/bin/python scripts/74_structural_neighbourhood.py
"""
from __future__ import annotations
import importlib.util, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate, typicality_gate  # noqa: E402

REPORT = ROOT / "reports/49_structural_neighbourhood.md"
K = 8


def main() -> int:
    spec = importlib.util.spec_from_file_location(
        "hyp", ROOT / "scripts/40_new_hypotheses.py")
    hyp = importlib.util.module_from_spec(spec)
    sys.modules["hyp"] = hyp
    spec.loader.exec_module(hyp)
    test, labels = load_test(), build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    cached = np.load(ROOT / "data/interim/contact_fingerprints.npz", allow_pickle=True)
    raw = cached["fingerprint"].astype(float)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}

    print("=== two scalings, because the 60 columns are on different units ===")
    print("  contact counts, capped minimum distances and CDR3 shares live on")
    print("  different scales, so an unscaled cosine would be dominated by the")
    print("  counts. Both the raw and the per-column standardised version are")
    print("  reported rather than one being chosen.\n")
    variants = {"raw cosine": raw,
                "column-standardised": (raw - raw.mean(0)) / (raw.std(0) + 1e-9)}
    series = {}
    for name, matrix in variants.items():
        normed = matrix / (np.linalg.norm(matrix, axis=1, keepdims=True) + 1e-9)
        values = {}
        for _, cohort in test.groupby("epitope_aa"):
            ids = [t for t in cohort["tcr_id"] if t in where]
            block = normed[[where[t] for t in ids]]
            d = 1.0 - block @ block.T
            np.fill_diagonal(d, np.inf)
            for i, tcr in enumerate(ids):
                values[tcr] = -float(np.sort(d[i])[:K].mean())
        series[name] = pd.Series(values)

    boltz = pd.read_csv(ROOT / "data/interim/boltz_features.csv").set_index("tcr_id")
    series["pae_tcr_pep_min (D's term)"] = boltz["pae_tcr_pep_min"]
    print(f"{'feature':<30}{'macro':>8}   per epitope")
    rows = []
    for name, values in series.items():
        per = {e: roc_auc_score(b["pos"], b["tcr_id"].map(values))
               for e, b in labels.groupby("epitope_aa")}
        rows.append({"feature": name, "macro": float(np.mean(list(per.values()))),
                     **{e[:3]: v for e, v in per.items()}})
        print(f"  {name:<28}{rows[-1]['macro']:>8.4f}   "
              + "  ".join(f"{e[:3]} {v:.3f}" for e, v in per.items()))
    table = pd.DataFrame(rows)

    print("\n=== does it add to configuration D? ===")
    submission = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv")
    d_score = submission.set_index("tcr_id")["prediction"]
    base = hyp.mean_auc(labels, d_score)
    best_name = max((r["macro"], r["feature"]) for r in rows
                    if "pae" not in r["feature"])[1]
    best = series[best_name]
    ranked = hyp.rank_within(best, test)
    rho = float(np.mean([spearmanr(b["tcr_id"].map(d_score),
                                   b["tcr_id"].map(ranked)).statistic
                         for _, b in test.groupby("epitope_aa")]))
    oracle = max((hyp.mean_auc(labels, (1 - w) * d_score + w * ranked), w)
                 for w in np.round(np.arange(0, 0.55, 0.05), 2))
    nested = hyp.nested_blend_gain(test, labels, d_score, best)
    ok, detail = private_variance_gate(test, best)
    gate = typicality_gate(test, labels, best)
    print(f"  best structural neighbourhood: {best_name}")
    print(f"  D {base:.4f}; within-epitope rho with D {rho:+.3f}")
    print(f"  oracle single global weight {oracle[0]:.4f} (w={oracle[1]:.2f}, "
          f"{oracle[0] - base:+.4f})")
    print(f"  nested LOEO {nested:.4f} ({nested - base:+.4f}) against +0.02")
    print(f"  private variance gate {ok} ({detail})")
    print(f"  typicality gate: {gate['verdict']} "
          f"({gate['feature_auc']:.4f} vs {gate['typicality_auc']:.4f}, "
          f"rho {gate['rho_with_typicality']:+.3f})")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("# E83 - contact fingerprint as a density space\n\n"
                      + table.round(4).to_markdown(index=False)
                      + f"\n\nBest: {best_name}. Into D: oracle {oracle[0]:.4f} "
                      f"at w={oracle[1]:.2f}, nested {nested:.4f}, rho with D "
                      f"{rho:+.3f}.\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
