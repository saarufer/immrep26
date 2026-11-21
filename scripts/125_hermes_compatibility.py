#!/usr/bin/env python
"""E148 -- HERMES interface compatibility: does this TCR's environment support this peptide?

An external review's third proposal, and the only route that passed all four
admission questions without being executed. HERMES (Visani et al.,
StatPhysBio/hermes, PMC12557503) is a 3D rotationally-equivariant network that
predicts amino-acid propensity at a masked focal residue from its local atomic
neighbourhood, pre-trained on ~10k CASP12 ProteinNet chains -- the general
protein universe, not binding data.

**The score**, as the review proposed it:

    dLL(t) = sum over peptide sites i of
             [ log P(aa_i | full TCR-pMHC environment)
             - log P(aa_i | same coordinates, TCR chains removed) ]

It asks whether the local environment *this TCR* provides makes the peptide's
own amino-acid composition more likely. Within an epitope cohort the peptide
sequence is constant, so everything that varies is the TCR and the environment
it creates -- which is exactly the within-epitope contrast the metric scores.

Against the admission standard:

  Q1  new knowledge source. General-protein-structure amino-acid statistics are
      distinct from neighbourhood density (SCEPTR/tcrdist), from folding
      confidence (PAE/ipTM, measured at -0.0193 as a residual), from polarity
      (polar_frac) and from docking quality (pDockQ, E145).
  Q2  needs only structures, and we hold PDBs for all 2,636 rows including all
      602 private ones.
  Q3  clean. Trained on masked amino-acid classification over general protein
      chains; it has never seen a TCR-pMHC binding label, so it is independent
      of the VDJdb annotation process this task audits.
  Q4  information.

**Fidelity note.** HERMES's `bp` models take hydrogens as a trained input
channel (`channels = C,N,O,S,H,SASA,charge`,
`get_physicochemical_info_for_hydrogens = True`), placed by `reduce` in the
released preprocessing. `reduce` 4.17 was built from source
(github.com/rlabduke/reduce) so this run matches the training distribution; a
pdbfixer-only substitute was used for timing and produced no recorded score.
Hydrogen placement is seeded (seed=0) because it is otherwise random.

**Pre-registered, before the run.** Residual-on-D **0.52-0.57**, nested
**-0.002 to +0.002**. Reasoning on record: four independent channels have now
landed at residual 0.5406-0.5609 and converted to nothing, and the closest
analogue -- E145's pDockQ, also structure-derived, also at 100% private
coverage, with the project's best leak-free survival (stratified 0.5717,
binomial p = 0.021) -- gave nested +0.0001.

Controls computed alongside, as the review asked: whether dLL is driven by
interface contact count, burial, or structure confidence rather than by
chemistry.

Run: .venv/bin/python scripts/125_hermes_compatibility.py
"""

from __future__ import annotations

import glob
import os
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
PRED = Path(os.environ.get(
    "BOLTZ_PRED", ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"))

# E152: this was hardcoded to a session-scoped scratchpad, which made the script
# unrunnable once that directory was collected. Point HERMES_MODELS at the
# `hermes_bp_000` directory of a `StatPhysBio/hermes` checkout.
MODELS = Path(os.environ.get(
    "HERMES_MODELS", ROOT / "data/external/hermes/trained_models/hermes_bp_000"))

OUT = ROOT / "data/interim/hermes_compat.csv"
# E152: the per-site 20-way distributions, which the first run computed and
# discarded. One pair of inferences per complex scores every single-residue
# peptide variant by lookup, so this is the cache any substitution-scanning work
# needs. The arrays hold LOG-PROBABILITIES, not raw logits: each row is
# log_softmax applied once to the seed-averaged logits. Keys are
# `<tcr_id>|logprob_full`, `<tcr_id>|logprob_pmhc`, `<tcr_id>|target_idx`,
# `<tcr_id>|site_ids`, `<tcr_id>|native_seq`, plus a single `aa_order`.
OUT_LOGITS = ROOT / "data/interim/hermes_logits.npz"
PROVENANCE = ROOT / "data/interim/hermes_provenance.json"

KEEP_PMHC = ("M", "P")
SEED = 0
BATCH = 64            # predict_from_pdbfile's batch size; was an unnamed literal


def strip_to_pmhc(src: Path, dst: Path) -> None:
    """Delete the TCR chains, leaving every pMHC coordinate untouched.

    CONECT records are dropped wholesale: after removing atoms they may point at
    serial numbers that no longer exist, and nothing downstream reads them.
    """
    keep = []
    for l in src.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            if len(l) > 21 and l[21] in KEEP_PMHC:
                keep.append(l)
        elif not l.startswith("CONECT"):
            keep.append(l)
    dst.write_text("\n".join(keep) + "\n")


def site_ids(res: dict) -> list[tuple]:
    """Ordered (resname, chain, resnum, insertion code) of the scored residues.

    The pdb-name field is dropped: the stripped copy is written under a different
    filename by construction. Channels that legitimately change when the TCR is
    removed -- SASA above all -- are not compared.
    """
    ids = res.get("res_ids")
    if ids is None:
        return []
    out = []
    for r in np.asarray(ids):
        f = [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in r]
        out.append((f[0], f[2], f[3], f[4]))
    return out


def logprob_of_target(res: dict) -> tuple:
    """Ensemble at logits level, then log_softmax, then pick the observed aa.

    NOTE the ensembling convention: logits are averaged across the 10 seeds and
    the softmax is taken once, which is not the same as averaging probabilities.
    E152 checks this against the authors' released scoring.

    Returns the summed and mean target log-probability, the site count, the
    observed-residue indices, and the full (n_sites, 20) log-probability matrix.
    """
    from scipy.special import log_softmax
    lg = np.asarray(res["logits"])            # (n_models, n_sites, 20)
    tg = np.asarray(res["targets"]).astype(int)
    lp = log_softmax(lg.mean(axis=0), axis=-1)   # (n_sites, 20)
    per_site = lp[np.arange(len(tg)), tg]
    return float(per_site.sum()), float(per_site.mean()), len(tg), tg, lp


def main() -> int:
    from hermes.inference.inference_hermes import (load_hermes_models,
                                                   predict_from_pdbfile)
    from hermes.utils.protein_naming import ind_to_ol_size
    dirs = sorted(glob.glob(str(MODELS / "*seed=*")))
    assert dirs, f"no model seeds under {MODELS}"
    models, hp, _ = load_hermes_models(dirs)
    reduce_ver = os.popen("reduce -version 2>&1 | head -1").read().strip()
    print(f"=== {MODELS.name}, {len(dirs)}-member ensemble, reduce "
          f"{reduce_ver} ===")

    # E152: the first run printed the reduce version and recorded nothing. A
    # structural score is only interpretable with the weights and the hydrogen
    # placer that produced it, so both are written next to the output.
    import hashlib
    import json
    digests = {}
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "*.pt"))):
            h = hashlib.sha256(Path(f).read_bytes()).hexdigest()
            digests[os.path.relpath(f, MODELS)] = h
    PROVENANCE.write_text(json.dumps({
        "models_dir": str(MODELS), "n_seeds": len(dirs),
        "reduce_version": reduce_ver, "hydrogen_seed": SEED,
        "ensembling": "mean of logits across seeds, then one log_softmax",
        "batch_size": BATCH, "boltz_predictions": str(PRED),
        "checkpoint_sha256": digests,
    }, indent=2))
    print(f"  provenance -> {PROVENANCE.relative_to(ROOT)} "
          f"({len(digests)} checkpoint file(s) hashed)")

    done = {}
    if OUT.exists():
        prev = pd.read_csv(OUT)
        # E152: the first version keyed on tcr_id alone, so a row that failed
        # once was skipped forever on every later resume. Only finite rows count
        # as done.
        ok = prev[prev["dll"].notna()]
        done = dict(zip(ok["tcr_id"], ok.to_dict("records")))
        n_retry = len(prev) - len(ok)
        print(f"  resuming: {len(done)} rows done"
              + (f", {n_retry} failed row(s) will be retried" if n_retry else ""))

    pdbs = sorted(glob.glob(str(PRED / "*" / "*_model_0.pdb")))
    print(f"  {len(pdbs)} complexes to score\n")
    tmp = Path(tempfile.mkdtemp())
    rec = list(done.values())
    # HERMES orders its 20 outputs by residue SIZE, not alphabetically. Any
    # consumer of this cache must read the axis through `aa_order`.
    logits: dict[str, np.ndarray] = {
        "aa_order": np.array("".join(ind_to_ol_size[i] for i in range(20)))}
    if OUT_LOGITS.exists():
        with np.load(OUT_LOGITS, allow_pickle=False) as z:
            logits.update({k: z[k] for k in z.files})
        n_cached = sum(1 for k in logits if k.endswith("|logprob_full"))
        print(f"  {n_cached} complexes already have cached log-probabilities")
    t0 = time.time()
    for i, p in enumerate(pdbs):
        p = Path(p)
        tid = p.parent.name
        if tid in done:
            continue
        try:
            pm = tmp / "pmhc.pdb"
            strip_to_pmhc(p, pm)
            rf = predict_from_pdbfile(str(p), models, hp, BATCH,
                                      regions={"whole_chains": ["P"]}, seed=SEED)
            rp = predict_from_pdbfile(str(pm), models, hp, BATCH,
                                      regions={"whole_chains": ["P"]}, seed=SEED)
            kf, kp = list(rf.keys())[0], list(rp.keys())[0]
            sf, mf, nf, tf, lpf = logprob_of_target(rf[kf])
            sp, mp, np_, tp, lpp = logprob_of_target(rp[kp])
            if nf != np_:
                raise ValueError(f"site count differs {nf} vs {np_}")
            if not np.array_equal(tf, tp):
                raise ValueError("observed peptide residues differ between "
                                 "the full and stripped complex")
            rec.append(dict(tcr_id=tid, ll_full=sf, ll_pmhc=sp,
                            dll=sf - sp, dll_mean=mf - mp, n_sites=nf))
            # (n_sites, 20) log-probabilities for both conditions, plus the
            # observed residue indices. Every single-residue peptide variant is
            # then a lookup rather than another fold.
            idf, idp = site_ids(rf[kf]), site_ids(rp[kp])
            if idf != idp:
                raise ValueError("site identity differs between the full and "
                                 "stripped complex")
            logits[f"{tid}|logprob_full"] = lpf.astype(np.float32)
            logits[f"{tid}|logprob_pmhc"] = lpp.astype(np.float32)
            logits[f"{tid}|target_idx"] = tf.astype(np.int16)
            logits[f"{tid}|site_ids"] = np.array(idf, dtype="U8")
            logits[f"{tid}|native_seq"] = np.array(
                "".join(ind_to_ol_size[int(i)] for i in tf))
        except Exception as e:
            rec.append(dict(tcr_id=tid, ll_full=np.nan, ll_pmhc=np.nan,
                            dll=np.nan, dll_mean=np.nan, n_sites=0))
            print(f"  [{tid}] failed: {type(e).__name__}: {e}")
        if (len(rec)) % 100 == 0:
            el = time.time() - t0
            pd.DataFrame(rec).to_csv(OUT, index=False)
            np.savez_compressed(OUT_LOGITS, **logits)
            n_new = len(rec) - len(done)
            rate = el / max(n_new, 1)
            print(f"  {len(rec)}/{len(pdbs)}  {rate:.2f}s/complex  "
                  f"eta {(len(pdbs)-len(rec))*rate/3600:.1f}h", flush=True)
    pd.DataFrame(rec).to_csv(OUT, index=False)
    np.savez_compressed(OUT_LOGITS, **logits)
    D = pd.DataFrame(rec)
    print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(D)} rows, "
          f"{int(D.dll.notna().sum())} finite)")
    n_cached = sum(1 for k in logits if k.endswith("|logprob_full"))
    print(f"wrote {OUT_LOGITS.relative_to(ROOT)}  "
          f"({n_cached} complexes x 2 conditions x (n_sites, 20) log-probabilities, "
          f"aa_order={logits['aa_order']})")
    print(D[["ll_full", "ll_pmhc", "dll", "dll_mean"]].describe().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
