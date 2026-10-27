#!/usr/bin/env python
"""E78 -- neighbourhood density against a V/J, length and Pgen matched background.

Proposed externally (as "E74" in a collaborator's plan; renumbered, E74 is
taken). The question: a TCR has many reference sequences nearby, but **given
its V gene, CDR3 length and generation background, how unusual is that?**

E8 ran the coarse version -- cohort density divided by density against the
other epitopes' TCRs -- and it lost: own-cohort 0.8396, difference 0.8357, log
ratio 0.8262, with background density *alone* scoring 0.60 in the positive
direction. Its settled mechanism was that crowded regions of TCR space are
where genuinely public clonotypes live, so normalising subtracts signal.

This is the refinement E8 did not run: the background is matched per query on
**TRBV stem x CDR3-beta length x paired-Pgen tercile**, so "crowded" is
measured against sequences of the same germline and generation class rather
than against the whole foreign pool. If E8's mechanism is the whole story the
sign will not change; if the coarse background was simply too coarse, it might.

Prediction recorded before running, because a prediction is not a measurement:
the same failure, because E77 independently found that local redundancy *is*
the signal -- de-duplicating the top-8 by family costs 0.0033 and the effective
family count scores 0.2957 on its own.

Run: .venv/bin/python scripts/68_matched_background_enrichment.py
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
from immrep.validation import private_variance_gate  # noqa: E402

REPORT = ROOT / "reports/43_matched_background.md"
TOPK = 8
MIN_BACKGROUND = 30


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

    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    vectors = cached["vectors"]
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    where = {t: i for i, t in enumerate(cached["tcr_id"])}

    pgen = pd.read_csv(ROOT / "data/interim/pgen.csv").set_index("tcr_id")
    frame = test.copy()
    frame["trbv"] = frame["TRBV_IMGT"].astype(str).str.split("*").str[0]
    frame["len_b"] = frame["cdr3_beta_aa"].astype(str).str.len()
    frame["pgen"] = frame["tcr_id"].map(pgen["log10_pgen_beta"])
    frame["pgen_tercile"] = pd.qcut(frame["pgen"], 3, labels=False,
                                    duplicates="drop")
    frame["len_bin"] = (frame["len_b"] // 2)
    index_of = {t: i for i, t in enumerate(frame["tcr_id"])}
    matrix = vectors[[where[t] for t in frame["tcr_id"]]]
    epitope = frame["epitope_aa"].to_numpy()

    def topk_density(row: int, candidates: np.ndarray) -> float:
        if len(candidates) < 2:
            return np.nan
        distances = 1.0 - matrix[candidates] @ matrix[row]
        distances = distances[candidates != row]
        k = min(TOPK, len(distances))
        return float(-np.sort(distances)[:k].mean())

    print("=== building the matched background per query ===")
    print("  background = TCRs of OTHER epitopes sharing TRBV stem, CDR3-beta")
    print(f"  length bin and Pgen tercile; relaxed if fewer than "
          f"{MIN_BACKGROUND}.\n")
    own, back, matched_on, pool_size = {}, {}, {}, {}
    keys = list(zip(frame["trbv"], frame["len_bin"], frame["pgen_tercile"]))
    for i, tcr in enumerate(frame["tcr_id"]):
        same = np.where(epitope == epitope[i])[0]
        own[tcr] = topk_density(i, same)
        foreign = epitope != epitope[i]
        for level, mask in (
            ("trbv+len+pgen", foreign
             & np.array([k == keys[i] for k in keys])),
            ("trbv+len", foreign & np.array([k[:2] == keys[i][:2] for k in keys])),
            ("len+pgen", foreign
             & np.array([(k[1], k[2]) == (keys[i][1], keys[i][2]) for k in keys])),
            ("len only", foreign & np.array([k[1] == keys[i][1] for k in keys])),
            # Final fallback: the whole foreign pool, which is E8's original
            # background. Without it a few hundred queries get no background
            # at all and the feature is NaN for them.
            ("any foreign (E8's background)", foreign),
        ):
            candidates = np.where(mask)[0]
            if len(candidates) >= MIN_BACKGROUND:
                back[tcr] = topk_density(i, candidates)
                matched_on[tcr] = level
                pool_size[tcr] = len(candidates)
                break
        else:
            back[tcr] = np.nan
            matched_on[tcr] = "none"
            pool_size[tcr] = 0
    levels = pd.Series(matched_on).value_counts()
    print("  match level achieved: " + ", ".join(
        f"{k} {v}" for k, v in levels.items()))
    print(f"  median background pool {np.median(list(pool_size.values())):.0f}")

    own_s = pd.Series(own)
    back_s = pd.Series(back)
    enrich = own_s - back_s
    rank = lambda s: hyp.rank_within(s, test)          # noqa: E731

    print("\n=== does the matched enrichment beat raw density? ===")
    rows = []
    for name, series in (("own-cohort top-8 (raw)", own_s),
                         ("matched background alone", back_s),
                         ("enrichment = own - background", enrich)):
        per = {e: roc_auc_score(b["pos"], b["tcr_id"].map(series))
               for e, b in labels.groupby("epitope_aa")
               if b["tcr_id"].map(series).notna().all()}
        rows.append({"feature": name, "macro": float(np.mean(list(per.values()))),
                     **{e[:3]: v for e, v in per.items()}})
        print(f"  {name:<32} macro {rows[-1]['macro']:.4f}   "
              + "  ".join(f"{e[:3]} {v:.3f}" for e, v in per.items()))
    table = pd.DataFrame(rows)

    print("\n=== and does it add anything to configuration D? ===")
    submission = pd.read_csv(
        ROOT / "submissions/d_density3_m30pool_m50_struct.csv")
    d_score = submission.set_index("tcr_id")["prediction"]
    base = hyp.mean_auc(labels, d_score)
    enrich_rank = rank(enrich)
    rho = float(np.mean([
        spearmanr(b["tcr_id"].map(d_score), b["tcr_id"].map(enrich_rank)
                  ).statistic for _, b in test.groupby("epitope_aa")]))
    print(f"  D as submitted {base:.4f}; enrichment's within-epitope rho with "
          f"D is {rho:+.3f}")
    grid = np.round(np.arange(0.0, 0.55, 0.05), 2)
    best = max((hyp.mean_auc(labels, (1 - w) * d_score + w * enrich_rank), w)
               for w in grid)
    nested = hyp.nested_blend_gain(test, labels, d_score, enrich)
    print(f"  oracle single global weight {best[0]:.4f} (w={best[1]:.2f}, "
          f"{best[0] - base:+.4f})")
    print(f"  nested LOEO blend gain      {nested:.4f} ({nested - base:+.4f}) "
          f"against the +0.02 bar")
    ok, detail = private_variance_gate(test, enrich)
    print(f"  private variance gate: {ok} ({detail})")

    print("\n=== verdict ===")
    if best[0] - base <= 0.02:
        print("  The prediction holds. A matched background does not rescue the")
        print("  design: an oracle allowed to pick the weight gains "
              f"{best[0] - base:+.4f},")
        print("  below the +0.02 bar, so no fitted version can clear it either.")
        print("  E8's mechanism survives the refinement -- crowding is signal,")
        print("  and dividing it out removes signal, however finely the")
        print("  background is matched.")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E78 - matched-background neighbourhood enrichment\n\n"
        "Background matched per query on TRBV stem x CDR3-beta length bin x "
        f"Pgen tercile, relaxed when fewer than {MIN_BACKGROUND} candidates.\n\n"
        "## Match levels\n\n" + levels.to_frame("n").to_markdown()
        + "\n\n## Feature performance\n\n"
        + table.round(4).to_markdown(index=False)
        + f"\n\n- D as submitted {base:.4f}\n"
        f"- enrichment's within-epitope rho with D: {rho:+.3f}\n"
        f"- oracle global weight: {best[0]:.4f} (w={best[1]:.2f})\n"
        f"- nested LOEO: {nested:.4f}\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
