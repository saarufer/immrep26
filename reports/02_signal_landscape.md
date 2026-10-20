# What predicts functional validation, and what does not

Status: 2026-10-20. All numbers are mean per-epitope ROC AUC over the five
public epitopes (1,782 scored rows), the official metric. Reference points:
organisers' tcrdist3 baseline **0.8247**, AlphaFold3 min-PAE baseline **0.787**.

The standard error of the public estimate is about **0.011**, so differences
below roughly 0.03 are not separable. On the private side only three epitopes
contribute, each with equal weight, giving a standard error near **0.023**.

## Confirmed facts

**The metric is mean per-epitope ROC AUC with ambiguous rows dropped.**
Scoring the organisers' own baseline file under every plausible metric matches
its 0.824 public score at exactly one: 0.824728. Pooled AUC gives 0.8577,
pooled AP 0.8602, mean per-epitope AP 0.8338, n-weighted per-epitope AUC
0.8466. Only within-epitope ranking is ever read, so cross-epitope calibration
is irrelevant — which is fortunate, because validation rates run from 28.7%
(NLVPMVATV) to 66.3% (GILGFVFTL).

**The public rows are the training set.** The five published epitopes' labels
come from TCRvdb, so 2,034 of the 2,636 rows in `test.csv` are labelled and the
remaining 602 across three unseen epitopes decide the ranking. The public
leaderboard measures a label join, not skill.

**Kaggle's "~68%" public split is the scored fraction.** 1,782 / 2,636 = 67.6%,
the 2,034 public rows minus 252 ambiguous. This reconciles the two figures in
the competition text and confirms ambiguous rows carry no weight.

## The signal has one dimension

Validated TCRs sit in tight, highly pure clusters; non-validating ones are
isolated. Median paired tcrdist to the nearest validated neighbour, by label:

| epitope | positives | negatives | singletons (d>90) | singleton pos-rate vs cohort |
|---|---:|---:|---:|---|
| GILGFVFTL | 24 | 169 | 34% | 0.31 vs 0.66 |
| YLQPRTFLL | 42 | 198 | 58% | 0.10 vs 0.41 |
| NLVPMVATV | 24 | 212 | 70% | 0.13 vs 0.29 |
| TTDPSFLGRY | 108 | 183 | 83% | 0.39 vs 0.50 |
| GLCTLVAML | 24 | 144 | 37% | 0.48 vs 0.66 |

Clusters at d<=90 with >=3 members have median label purity 0.89-1.00. So the
question every method is really answering is "does this TCR belong to a
convergent, antigen-driven cluster".

Every neighbourhood method extracts that one dimension equally well, and
saturates:

| representation | best density reduction | mean AUC | GLCTLVAML |
|---|---|---:|---:|
| SCEPTR `default` | kernel sigma=0.1 | **0.8396** | 0.713 |
| SCEPTR `small` | kernel sigma=0.1 | 0.8341 | 0.702 |
| tcrdist3 alpha+beta | kernel sigma=20 | 0.8329 | 0.699 |
| SCEPTR `large` | kernel sigma=0.1 | 0.8260 | 0.676 |
| SCEPTR `finetuned` | top-8 mean | 0.8206 | 0.715 |
| SCEPTR `a_sceptr` (alpha only) | top-11 mean | 0.8087 | **0.729** |
| SCEPTR `cdr3_only` | kernel sigma=0.1 | 0.8007 | 0.654 |
| SCEPTR `b_sceptr` (beta only) | kernel sigma=0.2 | 0.7519 | 0.583 |

A learned paired-chain representation beats a hand-built distance by 0.007,
inside one standard error. The metric space does not matter.

Two things in that table are worth keeping:

- **The alpha chain carries more of this signal than the beta chain**, in both
  tcrdist space (0.807 vs 0.746) and SCEPTR space (0.809 vs 0.752). That
  inverts the usual expectation that beta dominates specificity.
- **GLCTLVAML is the weak spot everywhere** (0.58-0.73), and it is the one
  epitope where alpha-only beats the paired model. Its singletons have a
  validation rate of 0.48 against a cohort rate of 0.66, so cluster membership
  is nearly uninformative there; density scores its singletons at AUC 0.482.
  Fixing one epitope from 0.69 to 0.85 moves the mean by 0.03, the largest
  single-epitope headroom available.

## Dead ends, with reasons

**Annotation metadata cannot rank the private epitopes.** No VDJdb provenance
field has usable within-epitope variance across all three held-back epitopes:
`vdjdb_score` is constant 0 for LLWNGPMAV and LTDEMIAQY, `method.verification`
is null for all 602 private rows, and LTDEMIAQY is 97.8% a single
`reference.id`. Details in `01_metadata_viability.md`.

Worse, where it does vary it points the wrong way. `vdjdb_score` scores
**0.455** pooled — below chance — and its direction flips between epitopes:
0.581 on GLCTLVAML but 0.399 on NLVPMVATV and 0.443 on YLQPRTFLL. The 158
score-3 NLVPMVATV entries validate at 15.8%. A model that learns "higher
confidence means more reliable" on one epitope is actively harmed on the next.
VDJdb's own confidence score anti-correlates with functional validation.

**Generation probability adds nothing.** OLGA Pgen scores 0.417 for
`-log10 Pgen(paired)`, i.e. *higher* Pgen means more likely to validate — the
opposite of the "false positives are typical repertoire draws" hypothesis. It
is not orthogonal either: Spearman correlation with density reaches -0.63 on
GILGFVFTL and -0.59 on TTDPSFLGRY, because dense convergent clusters are made
of easy-to-generate sequences. Within singletons its AUC is 0.50-0.53. The one
exception is GLCTLVAML singletons at 0.683 against density's 0.482, on n=62,
which is not enough to act on.

**Background normalisation hurts.** Dividing cohort density by density against
the other epitopes' TCRs — a TCRNET/ALICE-style enrichment — drops the score
from 0.8396 to 0.8262 (log ratio) or 0.8357 (difference). Background density
alone scores 0.60 *positive*: TCRs sitting in crowded regions of TCR space are
more likely to validate, not less, because crowded regions are where genuinely
public clonotypes live. Normalising subtracts signal.

## Where the remaining headroom has to come from

Density answers one question and answers it to saturation. Everything above
0.84 must come from information that is not a neighbourhood count:

1. **Epitope-conditioned specificity.** Nothing tried so far reads the epitope
   sequence at all. A model scoring `P(TCR binds this peptide-HLA)`, trained on
   external data with the target epitope held out, is neighbour-free and so can
   rank the singletons that are 34-83% of every cohort.
2. **Structure.** The organisers' AF3 min-PAE baseline reaches 0.787 alone,
   below density, but it is neighbour-free and therefore plausibly
   complementary. Whether density and structure fuse to beat either is
   untested and is the obvious experiment; cluster membership is a principled
   gate for such a fusion.
3. **Ensembling across representations**, which is cheap but should be expected
   to deliver the small gain that near-identical, highly correlated scores
   usually do.

## Later additions to the dead-end list

**Supervised fusion loses to variance.** Under leave-one-epitope-out, the best
single unsupervised feature scores 0.8362. Fitting on top of it makes things
worse, not better: logistic regression on density features alone 0.8359, on
density plus cluster shape 0.8299, a depth-2 gradient booster on both 0.8188.
Cluster shape on its own reaches only 0.6846. Five folds cannot support a
fitted combiner, which fixes a hard design constraint — **any fusion must be
rank-averaging with no fitted weights.** Full table in `03_loeo_fusion.md`.

Cluster shape is worth recording as a real but non-additive signal. Restricted
to rows that sit in a cluster of two or more, per-position CDR3-beta entropy
within the cluster scores 0.773 against density's 0.768, and internal
tightness scores 0.275 — i.e. **tighter clusters validate less**. That is the
convergent-recombination signature: a genuine antigen-driven motif family is
loose and diverse because many independent recombinations converged on it,
while a knot of near-duplicates that entered VDJdb together is tight. Density
scores both as "dense". The diagnostic that prompted this: NLVPMVATV's 23
high-density non-validating TCRs cluster tightly with each other (median
cosine 0.148) but sit far from that epitope's validated TCRs (0.312), whereas
on the other epitopes the high-density mistakes sit right on top of the
positives (GILGFVFTL 0.050, GLCTLVAML 0.088). Once imputed across the full
cohort, though, the signal collapses to 0.5553 — it is largely a restatement
of cluster membership.

**VDJdb enrichment cannot help the private side.** The cohorts could be padded
with same-epitope paired TCRs that were not screened, sharpening the density
estimate: 3,005 extra for GILGFVFTL and 9,881 for NLVPMVATV. But the private
epitopes are already exhausted — 11 extra for CINGVCWTV, 5 for LTDEMIAQY, 191
for LLWNGPMAV. Enrichment would inflate public scores up to 28-fold in
reference-set size while leaving two of the three deciding epitopes untouched,
destroying the correspondence between LOEO and private performance. The
organisers appear to have picked epitopes where VDJdb coverage is spent.

**Graph refinement destroys the signal, monotonically.** Degree centrality —
plain density — is the best of its family at 0.8362: PageRank 0.8329, one round
of label propagation 0.8297, three rounds 0.8095, ten rounds 0.7417,
eigenvector centrality 0.7486. More propagation is always worse, and
eigenvector centrality collapses onto the largest cluster and stops ranking the
others. The signal is strictly local.

**Rank-averaging across representations adds almost nothing**, as expected of
highly correlated scores: the best three-way average (SCEPTR default + small +
tcrdist3) reaches 0.8399 against 0.8362 for the best single, and that +0.0037
is the maximum over 200-plus combinations so the unbiased gain is nearer
+0.002. It is free and reduces sensitivity to the choice of representation, so
it is a reasonable default, but it is not progress.

## Current standing

Confirmed end to end: the submitted rank-normalised SCEPTR density scored
**0.839** on the public leaderboard against a local estimate of 0.8396. The
evaluation harness reproduces the leaderboard exactly, so no further submission
slots need to be spent measuring anything on the public epitopes — which
matters, because the daily submission limit is 1.

Seven independent negative results now point the same way: the extractable
signal in this task is local neighbourhood density, and nothing layered on top
of it survives extrapolation to a new epitope. The one hypothesis still open is
structure, which is neighbour-free by construction. All 2,636 TCR-peptide-MHC
complexes are being folded with Boltz-2 at 14 s each.

## Epitope-conditioned specificity: tried, and it is worse than chance

The one remaining non-structural idea was to read the peptide. A LightGBM
classifier over SCEPTR TCR embeddings concatenated with a BLOSUM62 peptide
encoding was trained on 21,112 VDJdb paired TCR-epitope pairs across 119
epitopes (capped at 1,500 per epitope so the two 30k-TCR epitopes could not
dominate), with negatives being the same TCRs paired with wrong epitopes so
that TCR identity alone carries no label. The target epitope was removed from
training entirely -- it has to be, since the test TCRs are themselves VDJdb
entries annotated to it, and a model that had seen them would call every one
positive.

Held-out per-epitope AUC: NLVPMVATV 0.553, GILGFVFTL 0.533, GLCTLVAML 0.531,
TTDPSFLGRY 0.507, YLQPRTFLL 0.283. **Mean 0.481 -- below chance.** Rank-averaged
with density it drags the mean from 0.8415 down to 0.7293.

Four epitopes sit at no-signal and YLQPRTFLL is strongly inverted. The likely
reading of that inversion: the model learns how *generic* a TCR looks against
VDJdb as a whole, and on YLQPRTFLL the generic-looking TCRs are precisely the
spurious annotations. That is the same failure shape as `vdjdb_score` -- a
quantity whose sign is not stable across epitopes. Predicting binding to a
peptide the model has never seen remains unsolved, and it contributes nothing
here. Cost: 78 min wall, 1,153 CPU-min.

That closes every non-structural hypothesis. Eight negative results, all
consistent with the signal being local neighbourhood density and nothing else.

Current best: rank-average of SCEPTR `default` and `small` kernel density,
**0.8415** mean per-epitope AUC (GILGFVFTL 0.870, GLCTLVAML 0.717, NLVPMVATV
0.846, TTDPSFLGRY 0.846, YLQPRTFLL 0.928), against the 0.8247 baseline.

## Structure: real, small, and not orthogonal

All 2,636 complexes folded (Boltz-2, single-sequence, 14 s each, 7.4 h).

As a standalone predictor the best structure feature reaches 0.759, below the
organisers' AlphaFold3 min-PAE baseline of 0.787 and well below density. But
the informative part is *which* features work. They all do, at about the same
level, including ones that cannot carry TCR-peptide specificity:

| feature | mean AUC |
|---|---:|
| `pae_a_mhc_min` (TCR alpha-MHC) | 0.759 |
| `pae_tcr_pep_min` (the organisers' statistic) | 0.759 |
| `complex_plddt` (global fold quality) | 0.694 |
| `pae_ab_min` (TCR alpha-TCR beta) | 0.656 |
| `iptm_mhc_pep` (**MHC-peptide**) | 0.647 |

`iptm_mhc_pep` is the diagnostic: the MHC-peptide interface involves no TCR at
all and is identical for every complex in a cohort, yet it still scores 0.647.
Structure confidence here is largely reporting *how canonical the TCR is* --
typical sequences fold confidently everywhere through the shared pair
representation -- and canonical means convergent means likely to validate.
Within-epitope Spearman correlation with density confirms it: `confidence_score`
+0.541, `pae_tcr_pep_min` +0.528, and even `iptm_mhc_pep` +0.262.

So structure is not the orthogonal signal it was hoped to be. It does still add
a little. Under nested LOEO with the blend weight chosen only on training
epitopes, blending `pae_tcr_pep_min` at weight 0.2-0.3 lifts 0.8392 to 0.8461,
an honest **+0.0070**, with the weight stably selected at 0.2-0.3 in every fold
and a nested-versus-oracle cost of only 0.0024.

Four of five epitopes improve (GLCTLVAML +0.026, NLVPMVATV +0.014, YLQPRTFLL
+0.006, GILGFVFTL +0.000). The single regression is TTDPSFLGRY at -0.011, and
it is the only 10-mer in the competition as well as the only public A*01:01.
Peptide length changes MHC-I groove geometry -- 10-mers bulge -- which is an
a-priori reason to expect it to behave differently, and **all three private
epitopes are 9-mers**. Restricted to the four 9-mers the gain is +0.0115.
LTDEMIAQY, the private A*01:01, has a structure-confidence distribution
(mean -1.40, sd 1.31) that looks like NLVPMVATV and YLQPRTFLL, where structure
helped, not like TTDPSFLGRY. One caution the other way: LLWNGPMAV's structure
confidence has sd 0.40 against 1.0-2.0 everywhere else, so structure has little
to rank with on that epitope.

## Selecting a configuration is itself overfitting

Every number quoted above for a *chosen* configuration is optimistic, because
sigma, which representations to average, and the structure weight were picked
by looking at public-epitope results. Density is unsupervised, so no label
touches the score -- but selection on five epitopes inflates an estimate just
as fitting does. `scripts/13_strict_loeo.py` separates the three quantities
over 323 candidate configurations:

| | mean per-epitope AUC |
|---|---:|
| oracle: best config scored on all five epitopes | 0.8508 |
| **nested LOEO: config re-chosen inside each fold** | **0.8460** |
| best a priori config, plain LOEO | 0.8508 |

Cross-validated selection is **0.0048 worse** than fixing the configuration in
advance. The TTDPSFLGRY fold shows the mechanism: on the other four epitopes,
`pae_a_mhc_min@0.3` looked best, and it returned 0.8328 on the held-out epitope
where the fixed configuration gave 0.8504. With five folds the selection step
has worse signal-to-noise than the thing being selected.

**Therefore: fix the configuration a priori and make as few decisions as
possible.** The defensible estimate of private performance is the nested
number, **0.846**, since the a priori configuration was chosen by someone who
had already seen the public results.
