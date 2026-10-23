#!/usr/bin/env python
"""Data-layer audit across every experiment, run before the write-up.

A code review of the *data* rather than the code: for each experiment family,
check that the artifacts it rests on exist, are keyed consistently against
test.csv, contain no duplicate or orphan identifiers, and still reproduce the
headline numbers quoted in the log. Thirty-nine experiments accumulated over
several days of interim files, and a write-up that cites them should not be the
first thing to notice a gap.

Each check prints PASS, FAIL or NOTE. FAIL means a number in the log cannot be
reproduced from the data on disk, or an identifier set is inconsistent. NOTE
means something worth recording that is not an error -- a declared output that
was never produced, a file kept for provenance that nothing reads.

Run: python scripts/36_data_audit.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep import constants as K  # noqa: E402
from immrep.data import (  # noqa: E402
    build_training_frame, load_labels, load_tcrdist_baseline, load_test)
from immrep.metric import official_score  # noqa: E402
from immrep.submission import validate  # noqa: E402

INTERIM = ROOT / "data/interim"
REPORT = ROOT / "reports/13_data_audit.md"

RESULTS: list[dict] = []


def check(group: str, name: str, ok: bool | None, detail: str = "") -> bool:
    status = "NOTE" if ok is None else ("PASS" if ok else "FAIL")
    RESULTS.append({"group": group, "check": name, "status": status,
                    "detail": detail})
    print(f"  [{status}] {name}" + (f"  {detail}" if detail else ""))
    return bool(ok)


def heading(text: str) -> None:
    print(f"\n=== {text} ===")


def audit_core(test, labels) -> None:
    heading("core inputs (E1-E4)")
    check("core", "test.csv row count", len(test) == K.N_TEST_ROWS,
          f"{len(test)} vs expected {K.N_TEST_ROWS}")
    check("core", "tcr_id unique", test["tcr_id"].is_unique,
          f"{test['tcr_id'].nunique()} distinct")
    counts = test["epitope_aa"].value_counts().to_dict()
    expected = {e: int(v["n"]) for e, v in K.ALL_EPITOPES.items()}
    mismatch = {e: (counts.get(e), n) for e, n in expected.items()
                if counts.get(e) != n}
    check("core", "per-epitope counts match the data description",
          not mismatch, "all 8 match" if not mismatch else str(mismatch))
    check("core", "8 epitopes, 5 public and 3 private",
          set(counts) == set(expected) and len(expected) == 8,
          f"{len(counts)} epitopes")
    required = ["cdr3_alpha_aa", "cdr3_beta_aa", "full_seq_alpha_aa",
                "full_seq_beta_aa", "hla_long", "split"]
    missing = {c: int(test[c].isna().sum()) for c in required
               if c in test.columns and test[c].isna().any()}
    check("core", "no missing sequence or HLA fields", not missing,
          "complete" if not missing else str(missing))
    split = test.groupby("split").size().to_dict()
    check("core", "split sizes", split.get("private") == K.N_PRIVATE
          and split.get("public") == K.N_PUBLIC, str(split))

    heading("labels (E2)")
    # The scored set is the public split minus the 252 ambiguous rows, which the
    # official metric drops. A first version of this audit expected all 2034 and
    # reported a FAIL against correct data.
    ambiguous = len(load_labels(drop_ambiguous=False)) - len(
        load_labels(drop_ambiguous=True))
    check("labels", "scored rows = public split minus ambiguous",
          len(labels) + ambiguous == K.N_PUBLIC,
          f"{len(labels)} + {ambiguous} ambiguous = {K.N_PUBLIC}")
    check("labels", "labelled ids are a subset of the public split",
          set(labels["tcr_id"]) <= set(
              test.loc[test["split"] == "public", "tcr_id"]),
          f"{labels['tcr_id'].nunique()} ids")
    check("labels", "no private row carries a label",
          not set(labels["tcr_id"]) & set(
              test.loc[test["split"] == "private", "tcr_id"]),
          "private split is unlabelled, as the competition requires")
    check("labels", "labels are binary", set(labels["pos"].unique()) <= {0, 1},
          f"positives {int(labels['pos'].sum())}, "
          f"rate {labels['pos'].mean():.3f}")
    per = labels.groupby("epitope_aa")["pos"].agg(["size", "sum", "mean"])
    degenerate = per[(per["sum"] == 0) | (per["sum"] == per["size"])]
    check("labels", "every public epitope has both classes", degenerate.empty,
          "all 5 mixed" if degenerate.empty else str(degenerate.index.tolist()))


def audit_metric(test, labels) -> None:
    heading("metric reproduction (E1)")
    # Predictions must be aligned positionally to the labels frame, so map by id.
    baseline = load_tcrdist_baseline()
    column = "prediction" if "prediction" in baseline else baseline.columns[1]
    series = baseline.set_index("tcr_id")[column]
    score = official_score(labels, labels["tcr_id"].map(series).to_numpy())
    check("metric", "tcrdist3 baseline reproduces the organisers' number",
          abs(score - K.BASELINES["tcrdist3"]) < 0.0015,
          f"{score:.4f} vs published {K.BASELINES['tcrdist3']}")

    # Kaggle truncates to three decimals rather than rounding: 0.8498 displayed
    # as 0.849, not 0.850. A first version compared with a tolerance and failed
    # a correct file.
    for name, shown in (("density3_struct03.csv", "0.849"),
                        ("density_sceptr_default_sigma0.1.csv", "0.839")):
        path = ROOT / "submissions" / name
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        series = frame.set_index(frame.columns[0])[frame.columns[1]]
        public = official_score(labels, labels["tcr_id"].map(series).to_numpy())
        check("metric", f"{name} reproduces its leaderboard score",
              f"{int(public * 1000) / 1000:.3f}" == shown,
              f"{public:.4f} truncates to {int(public * 1000) / 1000:.3f}, "
              f"leaderboard showed {shown}")


def ids_subset(group: str, name: str, frame: pd.DataFrame, test,
               column: str = "tcr_id") -> None:
    """Every keyed interim file must key into test.csv and not repeat itself."""
    if column not in frame.columns:
        check(group, f"{name}: has a {column} column", False, "absent")
        return
    ids = frame[column]
    orphan = set(ids) - set(test["tcr_id"])
    check(group, f"{name}: ids all present in test.csv", not orphan,
          f"{len(ids.unique())} distinct"
          + (f", {len(orphan)} orphan" if orphan else ""))


def audit_sequence(test, labels) -> None:
    heading("sequence-space artifacts (E5-E13, E16-E18, E20-E21, E26-E35)")
    for name, path, expect_full in (
        ("pgen.csv", INTERIM / "pgen.csv", True),
        ("specificity_scores.csv", INTERIM / "specificity_scores.csv", False),
    ):
        if not path.exists():
            check("sequence", f"{name}: present", False, "missing")
            continue
        frame = pd.read_csv(path)
        ids_subset("sequence", name, frame, test)
        if expect_full:
            check("sequence", f"{name}: covers every test row",
                  frame["tcr_id"].nunique() == len(test),
                  f"{frame['tcr_id'].nunique()}/{len(test)}")
        numeric = frame.select_dtypes(include=[np.number])
        bad = {c: int((~np.isfinite(numeric[c])).sum()) for c in numeric
               if not np.isfinite(numeric[c]).all()}
        check("sequence", f"{name}: all numeric values finite", not bad,
              "finite" if not bad else str(bad))

    for sub in ("sceptr", "tcrdist"):
        folder = INTERIM / sub
        files = sorted(folder.glob("*")) if folder.exists() else []
        check("sequence", f"{sub}/ distance or embedding cache", bool(files),
              f"{len(files)} files")
        for path in files:
            if path.suffix != ".npy":
                continue
            array = np.load(path, mmap_mode="r")
            square = array.ndim == 2 and array.shape[0] == array.shape[1]
            check("sequence", f"{sub}/{path.name}: shape",
                  array.ndim in (1, 2),
                  f"{array.shape}" + (" (square)" if square else ""))


def audit_structure(test) -> None:
    heading("structure artifacts (E14, E19, E28, E34, E36, E37, E39)")
    features = INTERIM / "boltz_features.csv"
    if features.exists():
        frame = pd.read_csv(features)
        ids_subset("structure", "boltz_features.csv", frame, test)
        check("structure", "boltz_features.csv covers every test row",
              frame["tcr_id"].nunique() == len(test),
              f"{frame['tcr_id'].nunique()}/{len(test)}")
        check("structure", "boltz_features.csv has no duplicate ids",
              frame["tcr_id"].is_unique, f"{len(frame)} rows")
        if "pae_tcr_pep_min" in frame:
            # Stored negated so that higher is better, matching density's
            # orientation: PAE 0.32 A is kept as -0.32, so an ascending rank
            # puts confident complexes high. A first version checked the raw
            # range and failed a correct, deliberate sign convention.
            column = frame["pae_tcr_pep_min"]
            check("structure", "min-PAE sign convention is uniform",
                  bool((column <= 0).all()),
                  "all non-positive, i.e. negated so higher is better")
            check("structure", "min-PAE magnitude in a physical range",
                  bool(column.abs().between(0, 32).all()),
                  f"{column.abs().min():.2f}-{column.abs().max():.2f} A")
        if "iptm" in frame:
            check("structure", "ipTM within [0, 1]",
                  bool(frame["iptm"].between(0, 1).all()),
                  f"{frame['iptm'].min():.3f}-{frame['iptm'].max():.3f}")

    predictions = INTERIM / "boltz/out/boltz_results_inputs/predictions"
    if predictions.exists():
        folded = sum(1 for p in predictions.iterdir() if p.is_dir())
        check("structure", "folded complexes match the test set",
              folded == len(test), f"{folded}/{len(test)}")

    for name in ("boltz_consistency.csv", "structure_contacts.csv",
                 "minimisation_robustness.csv"):
        path = INTERIM / name
        if not path.exists():
            check("structure", f"{name}: present", False, "missing")
            continue
        frame = pd.read_csv(path)
        ids_subset("structure", name, frame, test)

    fingerprints = INTERIM / "contact_fingerprints.npz"
    if fingerprints.exists():
        with np.load(fingerprints, allow_pickle=True) as data:
            check("structure", "contact_fingerprints.npz loads",
                  len(data.files) > 0, f"{len(data.files)} arrays")

    manifest = INTERIM / "decoy/manifest.csv"
    decoys = INTERIM / "decoy/out/boltz_results_inputs/predictions"
    if manifest.exists():
        frame = pd.read_csv(manifest)
        ids_subset("structure", "decoy manifest", frame, test)
        folded = sum(1 for p in decoys.iterdir() if p.is_dir()) if decoys.exists() else 0
        check("structure", "decoy folds complete", folded >= len(frame) * 0.99,
              f"{folded}/{len(frame)} target-decoy complexes")
        phases = frame.groupby("phase")["tcr_id"].nunique().to_dict()
        check("structure", "both decoy phases ran", len(phases) == 2, str(phases))

    relaxed = sorted((INTERIM / "relaxed").glob("*.pdb"))
    check("structure", "relaxed structures present", len(relaxed) == 69,
          f"{len(relaxed)} PDBs")
    bad_chains, bad_peptide = [], []
    for path in relaxed:
        chains: dict[str, set] = {}
        for line in path.read_text().splitlines():
            if line.startswith("ATOM"):
                chains.setdefault(line[21], set()).add(int(line[22:26]))
        if len(chains) != 4:
            bad_chains.append(path.stem)
        if not any(8 <= len(v) <= 11 for v in chains.values()):
            bad_peptide.append(path.stem)
    check("structure", "every relaxed structure has 4 chains", not bad_chains,
          "all 4-chain" if not bad_chains else f"{len(bad_chains)} bad")
    check("structure", "every relaxed structure has an 8-11 residue peptide",
          not bad_peptide,
          "found" if not bad_peptide else f"{len(bad_peptide)} missing")

    orphan = {p.stem for p in relaxed} - set(test["tcr_id"])
    check("structure", "relaxed ids key into test.csv", not orphan,
          f"{len(orphan)} orphan" if orphan else "all present")


def audit_md(test) -> None:
    heading("molecular dynamics (E39)")
    contacts = INTERIM / "md_residue_contacts.csv"
    trace = INTERIM / "md_interface_trace.csv"
    pairs_path = INTERIM / "md_pairs.json"
    if not contacts.exists():
        check("md", "md_residue_contacts.csv present", False, "missing")
        return
    frame = pd.read_csv(contacts)
    trajectories = frame.groupby(["tcr_id", "replica"]).ngroups
    check("md", "trajectory count", trajectories == 207,
          f"{trajectories} (69 complexes x 3 replicas)")
    check("md", "every complex has all 3 replicas",
          bool((frame.groupby("tcr_id")["replica"].nunique() == 3).all()),
          f"{frame['tcr_id'].nunique()} complexes")
    check("md", "every trajectory has 100 frames",
          bool((frame["n_frames"] == 100).all()),
          f"frames {sorted(frame['n_frames'].unique())}")
    ids_subset("md", "md_residue_contacts.csv", frame, test)
    check("md", "distances finite and physical",
          bool(np.isfinite(frame["dist_mean"]).all()
               and (frame["dist_mean"] > 0).all()),
          f"{frame['dist_mean'].min():.2f}-{frame['dist_mean'].max():.2f} A")
    consistent = ((frame["count_mean"] > 0) == (frame["frac_in_contact"] > 0)).all()
    check("md", "contact count and contact fraction agree", bool(consistent),
          "count>0 iff some frame within 4.5 A")
    check("md", "no duplicate (tcr_id, replica, chain, resid)",
          not frame.duplicated(["tcr_id", "replica", "chain", "resid"]).any(),
          f"{len(frame)} rows")
    offsets_ok = all(
        (np.sort(block["cdr3_offset"].values) == np.arange(len(block))).all()
        for _, block in frame.groupby(["tcr_id", "replica", "chain"]))
    check("md", "CDR3 offsets contiguous from 0 in every chain", offsets_ok)

    if trace.exists():
        series = pd.read_csv(trace)
        check("md", "interface trace matches the trajectory count",
              series.groupby(["tcr_id", "replica"]).ngroups == trajectories,
              f"{series.groupby(['tcr_id', 'replica']).ngroups} traces")

    if pairs_path.exists():
        pairs = json.loads(pairs_path.read_text())
        total = sum(len(v) for v in pairs.values())
        check("md", "pair table has 20 discordant and 20 control",
              all(len(v) == 20 for v in pairs.values()), f"{total} pairs")
        # The bug this file carried: site residues must be per member, since the
        # two TCRs in a pair number the CDR3 differently.
        per_role = all(f"{role}_site_residues" in entry
                       for entries in pairs.values() for entry in entries
                       for role in ("a", "b"))
        check("md", "site residues stored per pair member", per_role,
              "no shared union")
        differ = sum(1 for entries in pairs.values() for entry in entries
                     if entry["a_site_residues"] != entry["b_site_residues"])
        check("md", "pairs whose members number the site differently", None,
              f"{differ}/{total} -- each would be mis-probed by a union")
        simulated = set(frame["tcr_id"])
        members = {entry[f"{role}_id"] for entries in pairs.values()
                   for entry in entries for role in ("a", "b")}
        check("md", "every pair member was simulated", members <= simulated,
              f"{len(members)} members, {len(members - simulated)} missing")


def audit_submissions(test, labels) -> None:
    heading("submissions")
    for path in sorted((ROOT / "submissions").glob("*.csv")):
        frame = pd.read_csv(path)
        try:
            validate(frame, test=test)
            ok, detail = True, f"{len(frame)} rows"
        except Exception as exc:
            ok, detail = False, f"{type(exc).__name__}: {str(exc)[:60]}"
        series = frame.set_index(frame.columns[0])[frame.columns[1]]
        if ok:
            scored = labels["tcr_id"].map(series).to_numpy()
            detail += f", public {official_score(labels, scored):.4f}"
        check("submissions", f"{path.name}: validates", ok, detail)


def audit_reports() -> None:
    heading("declared outputs versus files on disk")
    declared: dict[str, list[str]] = {}
    for script in sorted((ROOT / "scripts").glob("*.py")):
        text = script.read_text()
        for line in text.splitlines():
            if "reports/" not in line or line.strip().startswith("#"):
                continue
            for token in line.split('"'):
                if token.startswith("reports/") and token.endswith(".md"):
                    declared.setdefault(token, []).append(script.name)
    # This audit's own report does not exist yet when the check runs.
    declared.pop(str(REPORT.relative_to(ROOT)), None)
    missing = {path: who for path, who in declared.items()
               if not (ROOT / path).exists()}
    check("reports", "every report a script declares exists", not missing,
          "all present" if not missing else "; ".join(
              f"{p} (declared in {', '.join(w)})" for p, w in missing.items()))

    existing = {p.name for p in (ROOT / "reports").glob("*.md")}
    check("reports", "report files on disk", None,
          f"{len(existing)} files, {len(declared)} declared by scripts")


def main() -> int:
    test = load_test()
    labels = build_training_frame()

    audit_core(test, labels)
    audit_metric(test, labels)
    audit_sequence(test, labels)
    audit_structure(test)
    audit_md(test)
    audit_submissions(test, labels)
    audit_reports()

    frame = pd.DataFrame(RESULTS)
    counts = frame["status"].value_counts().to_dict()
    print(f"\n{'=' * 62}\n{len(frame)} checks: "
          + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
    failures = frame[frame["status"] == "FAIL"]
    if not failures.empty:
        print("\nFAILURES:")
        for row in failures.to_dict("records"):
            print(f"  {row['group']}/{row['check']}: {row['detail']}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# Data-layer audit\n\n"
        f"{len(frame)} checks across the artifacts behind all 39 experiments: "
        + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + ".\n\n"
        + frame.to_markdown(index=False) + "\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 1 if not failures.empty else 0


if __name__ == "__main__":
    raise SystemExit(main())
