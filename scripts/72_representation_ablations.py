#!/usr/bin/env python
"""E81 -- twelve representations, including three trained on corrupted data.

SCEPTR ships fourteen variants and E6 tested seven. The untested ones include
published ablation controls that answer how much of the representation's value
is learned antigen biology: `synthetic_data` (trained on OLGA sequences, no
antigen selection), `shuffled_data` (alpha-beta pairing randomised) and
`dropout_noise_only`.

**All twelve use the same reduction (top-8)** so the comparison between them is
exact. E6's published numbers used each representation's own best reduction,
so they are not directly comparable row by row and are not mixed into the table
here; `default` appears in both and reads 0.8414 at top-8 against E6's 0.8396
at kernel sigma=0.1.

Gene names need standardising first: every variant rejects `TRAV14_DV4*01`.

Run: .venv/bin/python scripts/72_representation_ablations.py
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/47_representation_ablations.md"
K = 8
NAMES = ["default", "small", "large", "blosum", "left_aligned", "average_pooling",
         "tiny", "dropout_noise_only", "shuffled_data", "synthetic_data",
         "mlm_only", "finetuned"]


def main() -> int:
    import sceptr
    from sceptr import variant as V
    from libtcrlm.schema.tcr import TravGene, TrbvGene
    sceptr.disable_hardware_acceleration()
    test, labels = load_test(), build_training_frame()[["tcr_id", "epitope_aa", "pos"]]
    slash = lambda g: g.replace("_", "/") if isinstance(g, str) else None  # noqa: E731
    frame = pd.DataFrame({
        "TRAV": [slash(x) for x in test["TRAV_IMGT"]],
        "CDR3A": test["cdr3_alpha_aa"].tolist(),
        "TRAJ": [slash(x) for x in test["TRAJ_IMGT"]],
        "TRBV": [slash(x) for x in test["TRBV_IMGT"]],
        "CDR3B": test["cdr3_beta_aa"].tolist(),
        "TRBJ": [slash(x) for x in test["TRBJ_IMGT"]]})
    for column, vocab in (("TRAV", {g.name for g in TravGene}),
                          ("TRBV", {g.name for g in TrbvGene})):
        base = frame[column].fillna("").str.split("*").str[0]
        frame.loc[~base.isin(vocab), column] = None

    per = {}
    for name in NAMES:
        embedding = getattr(V, name)().calc_vector_representations(frame)
        embedding = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)
        position = {t: i for i, t in enumerate(test["tcr_id"])}
        top = {}
        for _, cohort in test.groupby("epitope_aa"):
            ids = cohort["tcr_id"].to_numpy()
            block = embedding[[position[t] for t in ids]]
            d = 1.0 - block @ block.T
            np.fill_diagonal(d, np.inf)
            for i, tcr in enumerate(ids):
                top[tcr] = -np.sort(d[i])[:K].mean()
        series = pd.Series(top)
        per[name] = {e: roc_auc_score(b["pos"], b["tcr_id"].map(series))
                     for e, b in labels.groupby("epitope_aa")}
    epitopes = sorted(labels["epitope_aa"].unique())
    table = pd.DataFrame([{"variant": n, **{e[:3]: per[n][e] for e in epitopes},
                           "macro": float(np.mean([per[n][e] for e in epitopes]))}
                          for n in NAMES]).sort_values("macro", ascending=False)
    print(table.round(4).to_string(index=False))

    print("\n=== nested: pick the best variant on four epitopes, use it on the fifth ===")
    held, picks = [], []
    for epitope in epitopes:
        inner = [e for e in epitopes if e != epitope]
        best = max(NAMES, key=lambda n: np.mean([per[n][i] for i in inner]))
        picks.append(best)
        held.append(per[best][epitope])
        print(f"  held out {epitope:<11} picked {best:<20} scored {per[best][epitope]:.3f}")
    nested = float(np.mean(held))
    fixed_best = float(np.mean([per['blosum'][e] for e in epitopes]))
    fixed_default = float(np.mean([per['default'][e] for e in epitopes]))
    print(f"\n  nested selection       {nested:.4f}")
    print(f"  fixed blosum (post hoc) {fixed_best:.4f}")
    print(f"  fixed default (D uses)  {fixed_default:.4f}")
    print(f"  selection cost          {nested - fixed_best:+.4f}")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("# E81 - twelve representations at a common top-8 reduction\n\n"
                      + table.round(4).to_markdown(index=False)
                      + f"\n\nNested selection {nested:.4f}; fixed blosum "
                      f"{fixed_best:.4f}; fixed default {fixed_default:.4f}; "
                      f"selection cost {nested - fixed_best:+.4f}.\n"
                      f"Folds picked: {picks}\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
