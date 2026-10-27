#!/usr/bin/env python
"""E79 -- weighting a neighbour's contribution by its annotation quality.

Different from what is closed: E4 scored `vdjdb_score` as a feature OF THE
QUERY and E33's S05 did the same for provenance. This weights THE NEIGHBOUR,
using only VDJdb metadata and no validation labels.

Run: .venv/bin/python scripts/70_reference_quality_weighting.py
"""
from __future__ import annotations
import importlib.util, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immrep.constants import VDJDB_VERSION  # noqa: E402
from immrep.data import build_training_frame, load_test  # noqa: E402
from immrep.validation import private_variance_gate  # noqa: E402

REPORT = ROOT / "reports/45_reference_quality.md"
K = 8


def main() -> int:
    spec = importlib.util.spec_from_file_location(
        "hyp", ROOT / "scripts/40_new_hypotheses.py")
    hyp = importlib.util.module_from_spec(spec)
    sys.modules["hyp"] = hyp
    spec.loader.exec_module(hyp)
    test, labels = load_test(), build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    full = pd.read_csv(ROOT / f"data/external/vdjdb/vdjdb-{VDJDB_VERSION}/vdjdb_full.txt",
                       sep="\t", low_memory=False)
    full = full[(full["species"] == "HomoSapiens")
                & full["cdr3.alpha"].notna() & full["cdr3.beta"].notna()]
    grouped = full.groupby(["cdr3.alpha", "cdr3.beta", "antigen.epitope"])
    quality = pd.DataFrame({"n_refs": grouped["reference.id"].nunique(),
                            "n_meth": grouped["method.identification"].nunique(),
                            "score": grouped["vdjdb.score"].max()})
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}

    print("=== field variability, where E3's gate fires ===")
    variability = []
    for epitope, cohort in test.groupby("epitope_aa"):
        key = cohort.set_index("tcr_id")[["cdr3_alpha_aa", "cdr3_beta_aa"]]
        rows = [quality.loc[(key.loc[t, "cdr3_alpha_aa"],
                             key.loc[t, "cdr3_beta_aa"], epitope)]
                for t in cohort["tcr_id"]
                if (key.loc[t, "cdr3_alpha_aa"], key.loc[t, "cdr3_beta_aa"],
                    epitope) in quality.index]
        block = pd.DataFrame(rows)
        variability.append({"epitope": epitope, "split": cohort["split"].iloc[0],
                            **{c: int(block[c].nunique()) for c in quality.columns}})
        print(f"  {epitope:<11} {cohort['split'].iloc[0]:<8} "
              + "  ".join(f"{c} {block[c].nunique()}" for c in quality.columns))
    variability = pd.DataFrame(variability)

    def build(column):
        out = {}
        for epitope, cohort in test.groupby("epitope_aa"):
            ids = cohort["tcr_id"].to_numpy()
            block = vectors[[where[t] for t in ids]]
            d = 1.0 - block @ block.T
            np.fill_diagonal(d, np.inf)
            key = cohort.set_index("tcr_id")[["cdr3_alpha_aa", "cdr3_beta_aa"]]
            if column is None:
                w = np.ones(len(ids))
            else:
                w = np.array([
                    float(quality.loc[(key.loc[t, "cdr3_alpha_aa"],
                                       key.loc[t, "cdr3_beta_aa"], epitope), column])
                    if (key.loc[t, "cdr3_alpha_aa"], key.loc[t, "cdr3_beta_aa"],
                        epitope) in quality.index else 0.0 for t in ids])
                # Floored at 1 so a score-0 neighbour is down-weighted, not erased:
                # an unfloored weight made np.average raise on all-zero windows.
                w = 1.0 + w
            w = w / w.mean()
            for i, tcr in enumerate(ids):
                order = np.argsort(d[i])[:K]
                out[tcr] = -float(np.average(d[i][order], weights=w[order]))
        return pd.Series(out)

    print("\n=== weighted density ===")
    base = build(None)
    baseline = hyp.mean_auc(labels, base)
    print(f"  unweighted {baseline:.4f}")
    results = [{"weighting": "unweighted", "macro": baseline, "delta": 0.0}]
    for column in ("n_refs", "score", "n_meth"):
        series = build(column)
        value = hyp.mean_auc(labels, series)
        results.append({"weighting": column, "macro": value,
                        "delta": value - baseline})
        print(f"  by {column:<8} {value:.4f}  {value - baseline:+.4f}")
    table = pd.DataFrame(results)
    best = build("n_refs")
    submission = pd.read_csv(ROOT / "submissions/d_density3_m30pool_m50_struct.csv")
    d_score = submission.set_index("tcr_id")["prediction"]
    start = hyp.mean_auc(labels, d_score)
    ranked = hyp.rank_within(best, test)
    oracle = max((hyp.mean_auc(labels, (1 - w) * d_score + w * ranked), w)
                 for w in np.round(np.arange(0, 0.55, 0.05), 2))
    ok, detail = private_variance_gate(test, best)
    print(f"\n  into D: {start:.4f} -> oracle {oracle[0]:.4f} (w={oracle[1]:.2f}, "
          f"{oracle[0] - start:+.4f})")
    print(f"  private variance gate {ok} ({detail}) -- but note the weight "
          "itself is constant on LTDEMIAQY;")
    print("  the weighted density inherits its variance from the distances.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("# E79 - reference annotation quality weighting\n\n"
                      + variability.to_markdown(index=False) + "\n\n"
                      + table.round(4).to_markdown(index=False)
                      + f"\n\nInto D: oracle {oracle[0]:.4f} at w={oracle[1]:.2f}.\n",
                      encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
