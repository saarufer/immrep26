# Reproducing this project

The predictor is frozen as of 2026-10-25. The claim this document has to
support is narrow and checkable: **`make submission` rebuilds the entered file
byte-for-byte**, and `make check` verifies that along with the data invariants
and the 57-check data audit.

```
make check          # guards, audit, rebuild + sha256 verify   (~5 min, CPU)
```

Everything else is opt-in, because the full pipeline is roughly 28 GPU-hours
and a day of CPU, and none of it is needed to reproduce the submission.

---

## 1. What you need before anything runs

### The licensed data, which is not in this repository

TCRvdb functional-validation labels are **request-only** from
`tcrvdb.nkiai.nl`. The licence does not permit redistribution in whole or in
part, so the file is gitignored and was never committed. Work using it must
cite `doi:10.1016/j.immuni.2026.09.007`.

Without it, `make test` still runs — the tests that need labels skip rather
than fail — but no score can be recomputed.

Everything else is fetched by `bash scripts/00_fetch_external.sh`: VDJdb
releases, IMGT/HLA allele sequences, and the PDB entries used in E34.

### Two Python environments, deliberately

```bash
uv venv --python 3.12 .venv
VIRTUAL_ENV=$PWD/.venv uv pip install -e ".[dev]"
```

That covers everything except two things that cannot live in the same
environment:

**Boltz-2** pins `numpy<2`, which conflicts with pyarrow. This is why every
interim file in this project is CSV rather than parquet. Install it into the
same venv only when you intend to fold:

```bash
VIRTUAL_ENV=$PWD/.venv uv pip install ".[fold]"
```

**OpenMM needs conda-forge, not PyPI.** The PyPI wheel ships an empty
`plugins/` directory — no CUDA plugin at all — so it silently falls back to
CPU and the MD stage becomes infeasible. This cost two wrong diagnoses before
the real cause was found. The working environment is separate:

```bash
micromamba create -y -p ./.mmenv-openmm -c conda-forge \
  python=3.12 "openmm>=8.2" cuda-version=12.8 pandas scipy
```

Verify it before trusting it:

```bash
PYTHONNOUSERSITE=1 micromamba run -p ./.mmenv-openmm python -c \
  "from openmm import Platform; print([Platform.getPlatform(i).getName() \
   for i in range(Platform.getNumPlatforms())])"
# must print [... 'CUDA'], not just ['Reference', 'CPU']
```

`PYTHONNOUSERSITE=1` is not optional: without it `~/.local` leaks a different
pandas into the conda environment and the import fails.

### Hardware this was run on

RTX 5090 (32 GB), 6 CPU threads left free for dataloaders. The GPU stages were
run **one process at a time**. Concurrency was measured and does not help:
four concurrent OpenMM processes total 1057 ns/day against 1040 for one, so
the implicit-solvent system already saturates the device.

---

## 2. Stage costs

| target | what it does | cost | needed for the submission? |
|---|---|---|---|
| `make test` | metric, data invariants, frozen-file guards | seconds | — |
| `make audit` | 57 data-layer checks | ~1 min | — |
| `make submission` | rebuild + sha256 verify | ~3 min | **yes** |
| `make qc` | id alignment, per-complex provenance | ~2 min | recommended |
| `make verify` | the analyses the write-up quotes | ~20 min | no |
| `make analysis` | every CPU analysis stage | hours | no |
| `make fold-main` | Boltz-2 over 2,636 complexes | 7.4 GPU-h | no* |
| `make fold-multi` | 5-sample re-fold (E19) | 13.2 GPU-h | no |
| `make fold-decoys` | 2,706 decoy folds (E36) | 7 GPU-h | no |
| `make relax` | 69 complexes minimised (E37) | ~1 CPU-h | no |
| `make md` | 414 ns across 207 trajectories (E39) | ~10 h | no |

\* `make submission` needs three caches. `data/interim/boltz_features.csv`
**is committed** — 1.8 MB of `tcr_id` plus 34 numeric descriptors, no sequences
and no labels — precisely so the 7.4 GPU-hours are not a prerequisite. The
other two are not committed because they are 29 MB and cheap:

```bash
.venv/bin/python scripts/05_sceptr_embeddings.py   # ~3 min
.venv/bin/python scripts/03_tcrdist_matrices.py    # ~5 min
```

Run those two on a fresh clone and `make submission` reproduces the frozen
file byte-for-byte without any GPU work.

---

## 3. Layout

```
src/immrep/        metric, data loading, features, validation, submission guard
scripts/           numbered stages; the number is the order they were run in
tests/             guards on the frozen predictor and the metric
reports/           00_experiment_log.md is the record; the rest are generated
submissions/       density3_struct03.csv is the frozen entry
data/external/     fetched or licensed inputs — not committed
data/interim/      derived caches
```

`reports/00_experiment_log.md` is the primary record: 42 experiments in the
order they were run, each with what it asked and what it settled, including
the ones that were later corrected or withdrawn. The compendium
(`reports/IMMREP26_experiment_compendium.{md,html}` and the `_zh` edition) is
a reorganised read of the same material.

---

## 4. Things that will bite you

These are recorded because each one cost real time here.

**Gene name format.** `test.csv` writes `TRAV14_DV4*01`; tcrdist3 and SCEPTR
want `TRAV14/DV4*01`. Convert `_` to `/`. After conversion, zero genes are
unmatched.

**SCEPTR rejects `TRBV23-1*01`** — an ORF, absent from its vocabulary. Exactly
one row is affected, in the private CINGVCWTV cohort. Set that chain's V gene
to `None` so the model falls back to CDR3 only, rather than dropping the row.

**OpenMM renames chains on write.** The input `M` and `P` become `C` and `D`,
so any analysis that keys on chain letter silently finds nothing. Identify the
peptide by residue count (8–11) instead.

**Density reference sets.** Density must be computed over the **whole cohort**,
not just the labelled rows. Computing it on labelled rows alone is a slip that
reverses the operator ordering in E6's table.

**Clustering bases are not interchangeable.** `scripts/20` clusters labelled
rows; `scripts/25` clusters the whole cohort. Mixing error rates from one with
compositions from the other is the mistake E40 found in E30.

**Do not compute results in an unsaved heredoc.** This cost the project five
times, most seriously in E30, whose withdrawn 0.80 estimate could not be
reconstructed because the code did not exist. If a number reaches a report, the
code that produced it belongs in `scripts/`.

**`pgrep -f <pattern>` matches its own command line.** A wait loop written as
`until ! pgrep -f 'my_script'` never exits, because the loop's own shell
contains the pattern. This produced a monitor that would never have fired.
