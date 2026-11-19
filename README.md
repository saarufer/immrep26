# IMMREP26 — predicting which reported TCR–pMHC pairs validate

Kaggle competition [`immrep26`](https://www.kaggle.com/competitions/immrep26),
closed 2026-11-02. **Final private score 0.819** (submission `D`), against the
organisers' tcrdist3 baseline at 0.768 and a best-in-field 0.846 among entries
that did not recover the labels.

Every TCR in the test set is already *reported* in VDJdb as recognising its
annotated peptide–HLA. The task is to rank them by whether that annotation
holds up in the Schumacher lab's standardised functional screen, where roughly
44% of re-tested pairs showed no reactivity. So this is not specificity
prediction against random decoys — it is **predicting annotation reliability
among historically reported positives**.

## What the competition measured

`test.csv` holds 2,636 paired αβ TCRs over 8 peptides in one file, split by
Kaggle into two leaderboards:

| | peptides | rows | role |
|---|---|---:|---|
| public | GILGFVFTL, YLQPRTFLL, TTDPSFLGRY, NLVPMVATV, GLCTLVAML | 2,034 (1,782 scored) | labels published via TCRvdb → **this is the training set** |
| private | CINGVCWTV, LLWNGPMAV, LTDEMIAQY | 602 (243 / 229 / 130) | decided the final ranking |

The metric is **mean per-peptide ROC AUC with ambiguous rows dropped**,
reverse-engineered from the organisers' baseline file and then confirmed
end-to-end: a submission predicted locally at 0.8396 scored 0.839. Only the
ranking *within* a peptide is read, and each peptide carries equal weight
regardless of size. `src/immrep/metric.py` documents the candidates ruled out.

Two consequences, both of which the final board confirmed:

- **A near-perfect public score is a label join, not skill.** 19 teams reached
  ≥0.984 on the public board. Public standing was dropped from consideration
  after E2 and never optimised.
- **Differences below ~0.047 between unrelated methods are not resolvable**
  on 602 rows and three cohorts. Measured, not asserted:
  `scripts/129_leaderboard_resolution.py` → `reports/108_*`.

## Results

`reports/107_final_results.md` is the post-mortem. Entries from this repository:

| submission | public | private |
|---|---:|---:|
| `C1_consensus_S4_FHB.csv` — rank consensus with a collaborator's model | 0.882 | 0.811 |
| `S4_D_stack9.csv` — D plus nine refinements | 0.871 | 0.810 |
| **`d_density3_m30pool_m50_struct.csv` — D, four equal-weight terms** | 0.860 | **0.819** |
| organisers' tcrdist3 baseline | — | 0.768 |

**The public ordering reversed on the private peptides.** The two more
elaborate entries bought public score and lost private score against the
plainest member of the family. For C1 the swing is quantified: its expected
margin over D on a private-shaped bootstrap was +0.022 ± 0.006, and the
realised margin was −0.008, about five paired standard deviations away. The
pre-registered private estimate for D, recorded in its submission description
on 2026-10-26, was 0.798–0.802 against a realised 0.819.

Counting by team — collapsing multi-account clusters, removing the seven
label-recovery entries — this is **19th of 107**, inside a 35-team tie that
spans the top third of the clean field.

## The method

Validated TCRs sit in tight, label-pure clusters and non-validating ones are
isolated, so every method here is really answering *is this TCR part of a
convergent, antigen-driven cluster*. `reports/02_signal_landscape.md` is the
main evidence document; that signal is one-dimensional and it saturates.

**D** is the equal-weight average of four within-peptide percentile ranks,
with **no fitted parameters** in the aggregation:

| component | definition |
|---|---|
| `density3` | exponential kernel density in default-SCEPTR (bw 0.1), small-SCEPTR (0.1) and paired α+β tcrdist3 (20) spaces; self excluded |
| `pooled_support` | negative mean of the eight smallest default-SCEPTR cosine distances in a same-peptide reference pool (cohort + paired human VDJdb 2026-06-03) |
| `full_sequence_support` | negative mean of the eight smallest α+β full-V-domain alignment distances within the cohort |
| `structure_support` | negative minimum TCR–peptide PAE from Boltz-2 complexes — model confidence, not affinity |

No binding labels, no peptide supervision, and **no VDJdb-trained binding
predictor** — ruled out in E121 as circular for a task auditing VDJdb (the
strongest public candidate overlapped the test set 100% and carried 872
inverted labels).

A single SCEPTR density alone is 0.8396 of D's 0.8601, i.e. **94.3% of the
above-chance signal**; the ladder is 0.8396 → +structure 0.8498 → +V-domain
0.8557 → +pool 0.8601. The premise fails where cohort convergence is absent:
density alone reaches 0.927 on YLQPRTFLL and **0.677 on GLCTLVAML**, and
nothing tested across 150 experiments recovered that cohort.

## What does not help

Each with a mechanism in `reports/00_experiment_log.md`: VDJdb annotation
metadata (actively harmful — `vdjdb_score` scores 0.485 and flips sign between
peptides), OLGA generation probability, the choice of representation,
background-density normalisation, cluster-shape features, supervised fusion,
graph centrality and label propagation, an unseen-peptide specificity model
(0.481, below chance), multi-sample docking consistency (+0.001 for 13
GPU-hours), cluster-conditional blend weights, interface polarity beyond one
structure, HERMES compatibility, pDockQ/ipSAE, expression-framework residues,
and the spectral-physicochemical correction head of E150.

The error was decomposed rather than only reduced
(`reports/06_auc_error_decomposition.md`): of the remaining 0.15, **62% comes
from pairs where neither TCR is clustered** and **21% from a class covering
4.5% of pairs** — isolated positive against clustered negative — where the
model is inverted and errs 67% of the time. Reweighting cannot reach it:
helping that class costs ~13× more elsewhere, because scores are per TCR.

Five peptides cannot support a fitted combiner — logistic regression on top of
density loses 0.006, gradient boosting 0.021 — so **any fusion here is a
zero-parameter rank average**, and a component enters only if it adds >0.02 in
leave-one-peptide-out. Enforced in `scripts/12_final_fusion.py`. Having no
fitted parameters is not the same as having no selection exposure: the
architecture, features and blend weight were chosen in light of labelled data,
and `scripts/13_strict_loeo.py` measures that exposure at 0.0048 over 323
candidate configurations.

## A note on the benchmark

The top seven private entries scored 0.953–0.999, then the field drops to
0.846 and continues smoothly down. A 0.999 over 602 rows and three unseen
peptides is not attainable by prediction.

No leak exists inside the released files — `tcr_id` plate coordinates, all 17
metadata columns, missingness patterns and cross-peptide clonotype overlap
were all checked and are noise. But NCBI BioProject
[PRJNA1464244](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1464244),
registered **2026-05-10**, holds ONT sequencing of the pooled VDJdb-10 screens
for **all eight peptides**, including `CINGVCWTV_CD69_pos_rep_1–3`,
`LLWNGPMAV_CD69_pos_rep_1–3` and `LTDEMIAQY_CD69_pos_rep_1–3` with matched
baselines. With `test.csv`'s full-length sequences to map against and the five
released TCRvdb peptides to calibrate a threshold, the private labels are
derivable from public data.

This **corrects E2**, which concluded the private split was sound after
checking only the TCRvdb repository and its changelog, not the sequence
archives. No reconstruction was attempted here and nothing is attributed to
any team. Details in `reports/107_final_results.md` §4.

## Layout

```
src/immrep/      metric, data loading, features, LOEO harness, submission QC
scripts/         numbered, each runnable standalone, 00 → 129
reports/         findings, regenerated by the scripts that produced them
data/            all gitignored; see data/external/README.md for provenance
```

Validation is leave-one-peptide-out throughout. Random splits leak badly here:
near-identical clonotypes sit on both sides of a random split. The 18
pseudo-peptide folds used late in the project split by (peptide, source study)
and aggregate with equal weight per peptide; their permutation null has
sd 0.0045.

## Reproducing

```bash
uv venv --python 3.12 .venv && source .venv/bin/activate
uv pip install torch --index-url https://download.pytorch.org/whl/cu128
uv pip install numpy pandas scipy scikit-learn lightgbm tcrdist3 sceptr kaggle biopython
./scripts/00_fetch_external.sh            # VDJdb; TCRvdb must be requested
kaggle competitions download -c immrep26 -p data/raw && unzip data/raw/immrep26.zip -d data/raw
python scripts/02_confirm_metric.py       # reproduces the metric derivation
python scripts/05_sceptr_embeddings.py
python scripts/13_strict_loeo.py          # honest estimate, incl. selection cost
python scripts/18_final_submission.py     # builds the submitted file
python scripts/21_submission_qc.py        # id-alignment assertions
python scripts/129_leaderboard_resolution.py   # what the board can resolve
```

Structure prediction (`scripts/08`–`10`) needs `boltz`, which pins `numpy<2`
and so conflicts with `pyarrow`; interim files are CSV to avoid the clash.
Folding all 2,636 complexes takes ~10 h on an RTX 5090 at 14 s each.
`data/interim/boltz_features.csv` is the one committed cache, so
`make submission` works from a fresh clone without a GPU.

## Data terms

TCRvdb is request-only from [tcrvdb.nkiai.nl](https://tcrvdb.nkiai.nl) and its
terms **prohibit redistribution**, so it is gitignored and never committed.
Work using it must cite [doi:10.1016/j.immuni.2026.09.007](https://doi.org/10.1016/j.immuni.2026.09.007).
VDJdb carries its own CC-BY-NC-SA terms. Reports contain aggregate statistics
derived from TCRvdb, not records.
