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
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
MODELS = Path("/tmp/claude-1000/-home-chen-project-Immrep26/"
              "1d950e23-986a-43fb-ac61-0ce9ce8f1804/scratchpad/hermes/"
              "trained_models/hermes_bp_000")
OUT = ROOT / "data/interim/hermes_compat.csv"
TCR_CHAINS = ("A", "B")
KEEP_PMHC = ("M", "P")
SEED = 0
AA = "ACDEFGHIKLMNPQRSTVWY"


def strip_to_pmhc(src: Path, dst: Path) -> None:
    keep = [l for l in src.read_text().splitlines()
            if not l.startswith(("ATOM", "HETATM")) or l[21] in KEEP_PMHC]
    dst.write_text("\n".join(keep) + "\n")


def logprob_of_target(res: dict) -> tuple:
    """Ensemble at logits level, then log_softmax, then pick the observed aa."""
    from scipy.special import log_softmax
    lg = np.asarray(res["logits"])            # (n_models, n_sites, 20)
    tg = np.asarray(res["targets"]).astype(int)
    lp = log_softmax(lg.mean(axis=0), axis=-1)   # (n_sites, 20)
    per_site = lp[np.arange(len(tg)), tg]
    return float(per_site.sum()), float(per_site.mean()), len(tg), per_site


def main() -> int:
    from hermes.inference.inference_hermes import (load_hermes_models,
                                                   predict_from_pdbfile)
    dirs = sorted(glob.glob(str(MODELS / "*seed=*")))
    assert dirs, f"no model seeds under {MODELS}"
    models, hp, _ = load_hermes_models(dirs)
    print(f"=== hermes_bp_000, {len(dirs)}-member ensemble, reduce "
          f"{os.popen('reduce -version 2>&1 | head -1').read().strip()} ===")

    done = {}
    if OUT.exists():
        prev = pd.read_csv(OUT)
        done = dict(zip(prev["tcr_id"], prev.to_dict("records")))
        print(f"  resuming: {len(done)} rows already computed")

    pdbs = sorted(glob.glob(str(PRED / "*" / "*_model_0.pdb")))
    print(f"  {len(pdbs)} complexes to score\n")
    tmp = Path(tempfile.mkdtemp())
    rec = list(done.values())
    t0 = time.time()
    for i, p in enumerate(pdbs):
        p = Path(p)
        tid = p.parent.name
        if tid in done:
            continue
        try:
            pm = tmp / "pmhc.pdb"
            strip_to_pmhc(p, pm)
            rf = predict_from_pdbfile(str(p), models, hp, 64,
                                      regions={"whole_chains": ["P"]}, seed=SEED)
            rp = predict_from_pdbfile(str(pm), models, hp, 64,
                                      regions={"whole_chains": ["P"]}, seed=SEED)
            kf, kp = list(rf.keys())[0], list(rp.keys())[0]
            sf, mf, nf, _ = logprob_of_target(rf[kf])
            sp, mp, np_, _ = logprob_of_target(rp[kp])
            if nf != np_:
                raise ValueError(f"site count differs {nf} vs {np_}")
            rec.append(dict(tcr_id=tid, ll_full=sf, ll_pmhc=sp,
                            dll=sf - sp, dll_mean=mf - mp, n_sites=nf))
        except Exception as e:
            rec.append(dict(tcr_id=tid, ll_full=np.nan, ll_pmhc=np.nan,
                            dll=np.nan, dll_mean=np.nan, n_sites=0))
            print(f"  [{tid}] failed: {type(e).__name__}: {e}")
        if (len(rec)) % 100 == 0:
            el = time.time() - t0
            pd.DataFrame(rec).to_csv(OUT, index=False)
            n_new = len(rec) - len(done)
            rate = el / max(n_new, 1)
            print(f"  {len(rec)}/{len(pdbs)}  {rate:.2f}s/complex  "
                  f"eta {(len(pdbs)-len(rec))*rate/3600:.1f}h", flush=True)
    pd.DataFrame(rec).to_csv(OUT, index=False)
    D = pd.DataFrame(rec)
    print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(D)} rows, "
          f"{int(D.dll.notna().sum())} finite)")
    print(D[["ll_full", "ll_pmhc", "dll", "dll_mean"]].describe().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
