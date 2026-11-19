# IMMREP26 — Experiment Compendium

**Validation of public TCR–pMHC data.** Every experiment run on this project, with
what it was asking, how it was measured, and what it settled.

Compiled 2026-10-23, after the data-layer audit in `reports/13_data_audit.md`
(57 checks, all passing). Source of record: `reports/00_experiment_log.md`.

---

## How to read the numbers

The competition metric is the **unweighted mean of per-epitope ROC AUC, with
ambiguous rows dropped** — reverse-engineered in E1 and confirmed three times
against the leaderboard. Unless stated, every figure below is that quantity over
the five public epitopes.

| reference point | value |
|---|---:|
| organisers' tcrdist3 baseline | **0.8247** |
| organisers' AlphaFold3 min-PAE baseline | 0.787 |
| standard error of a five-epitope estimate | ≈ 0.011 |
| standard error on the three private epitopes | ≈ 0.023 |
| pre-registered inclusion threshold for a new component | **> +0.02 in LOEO** |

The dataset is 2,636 rows in `test.csv`: 2,034 public (of which 1,782 carry a
usable label and 252 are ambiguous) and 602 private across three epitopes. Only
the private 602 decide the final ranking.

### The submitted model

`submissions/density3_struct03.csv`, built by `scripts/18_final_submission.py`:
a rank average of three within-epitope kernel densities — SCEPTR `default` and
`small` at cosine σ=0.1, tcrdist3 α+β at σ=20 — blended 0.7/0.3 with Boltz-2
TCR–peptide minimum PAE. **Fully unsupervised, zero fitted parameters.** Public
0.8498 (leaderboard 0.849).

### The honest private estimate

**0.846**, with a documented but unquantified downside risk. E30 argued for
0.80; **E40 withdrew that number** after it failed a leave-one-epitope-out
transfer test at MAE 0.083 — nearly twice the bias it claimed to detect. The
distribution shift E30 found is real; its size is not measurable with five
epitopes. Treat 0.67–0.91 as a scenario range, never a confidence interval.

---

## Thematic index

| theme | experiments |
|---|---|
| Foundations: metric, labels, metadata | E1, E2, E3, E4 |
| What the signal actually is | E5, E6, E9, E20, E32 (D10) |
| Rejected sequence approaches | E7, E8, E10, E11, E13, E16, E17, E32 |
| Candidate-pool expansion | E12, E27, E33, E35, E38 |
| Structure | E14, E19, E28, E34, E36, E37, E39 |
| Validation methodology and error anatomy | E15, E18, E21, E22, E23, E24, E29, **E30** |
| Mechanism of the irreducible error | E20, E31, E37, E39 |
| Independent cross-check | E25, E26 |

---

# Foundations

## E1 — Is the competition's metric knowable without the Evaluation page?

**Goal.** Determine the exact scoring function, so local measurement can replace
submissions (the daily limit is 1).

**Method.** Score the organisers' own `tcrdist_predictions.csv`, whose public
score is 0.824, against the published TCRvdb labels under every plausible metric.
`scripts/02_confirm_metric.py`.

**Result.** Exactly one candidate matches.

| candidate | value |
|---|---:|
| **mean per-epitope ROC AUC, ambiguous dropped** | **0.824728** |
| mean per-epitope average precision | 0.8338 |
| n-weighted mean per-epitope ROC AUC | 0.8466 |
| pooled ROC AUC | 0.8577 |
| pooled average precision | 0.8602 |
| mean per-epitope AUC, ambiguous as positive | 0.7908 |
| mean per-epitope AUC, ambiguous as negative | 0.7769 |

**Conclusion.** Only within-epitope ranking is read, so the 28.7–66.3% spread in
validation rate across epitopes needs no calibration. Each epitope carries equal
weight regardless of size, so the 130-row LTDEMIAQY counts as much as the 243-row
CINGVCWTV. Confirmed end to end twice: a submission predicted locally at 0.8396
scored 0.839, and 0.8498 scored 0.849 — Kaggle truncates rather than rounds.
**The local harness is the leaderboard.** This also resolves the competition
text's internal contradiction: the "~68%" public split is 1,782/2,636 = 67.6%,
i.e. 2,034 public rows minus 252 ambiguous.

## E2 — Can the public labels simply be recovered?

**Goal.** Test whether the public leaderboard measures skill at all.

**Method.** Join TCRvdb to `test.csv` and score the recovered label.

**Result.** `tcr_id` is shared **verbatim**: 2,034/2,034 rows match exactly, no
fuzzy matching needed. Submitting the recovered label scores **1.000 on every
public epitope**. Simulating 5% match failure still gives 0.9988.

**Conclusion.** The public leaderboard measures a `merge`, not skill — the 1.000
and 0.997 standings are reproducible in one line. Public standing was dropped
from consideration entirely after this. **Follow-up: do the private labels leak
too? No.** TCRvdb's changelog lists exactly which epitopes each release
contains: YLQ+GLC (2025-05-01), plus TTD, GIL and NLV (2026-10-01). CINGVCWTV,
LLWNGPMAV and LTDEMIAQY have never been released and are absent from the
downloaded file. **The private leaderboard is sound**, so modelling decides the
ranking.

## E3 — Can annotation provenance rank the private epitopes at all?

**Goal.** Establish, before trusting any cross-validation score, whether
provenance metadata even varies where it must.

**Method.** For every VDJdb metadata field, count distinct levels *within* each
epitope cohort; flag a field unusable if it is constant, or if its dominant level
holds >95% of rows, in any private epitope.
`scripts/01_audit_metadata_viability.py`.

**Result.** **Not one field survives.** `vdjdb_score` is constant 0 for
LLWNGPMAV and LTDEMIAQY; `method.verification` is null for all 602 private rows;
LTDEMIAQY is 97.8% a single `reference.id` and a single `method.identification`.

**Conclusion.** A provenance feature can look strong in cross-validation on the
public epitopes — where GILGFVFTL and NLVPMVATV span `vdjdb_score` 0–3 — and
contribute literally zero on the private side. **This became a standing rule:
every new feature is checked for within-epitope variance on all three private
epitopes before its LOEO score is believed.**

## E4 — Metadata-only baseline: how much does evidence quality explain?

**Goal.** Quantify what VDJdb's own confidence score is worth.

**Method.** Score `vdjdb_score` directly, pooled and within epitope.

**Result.** Pooled AUC **0.455** — below chance — and the direction *flips
between epitopes*: 0.581 on GLCTLVAML, 0.399 on NLVPMVATV, 0.443 on YLQPRTFLL,
undefined on TTDPSFLGRY (constant). The 158 score-3 NLVPMVATV entries validate
at 15.8%, against 23.5% for score-0 in the same cohort.

**Conclusion.** VDJdb's confidence score **anti-correlates** with functional
validation. A model that learns "higher confidence means more reliable" on one
epitope is actively harmed on the next. Metadata is kept as a documented
ablation, never as a submission component.

---

# What the signal actually is

## E5 — What is the signal, structurally?

**Goal.** Characterise what separates validated from non-validating TCRs.

**Method.** Within each epitope, compare median distance from each TCR to its
nearest validated neighbour, split by label; cluster at tcrdist ≤ 90 and measure
label purity and singleton rate.

**Result.**

| epitope | pos → nearest validated | neg → nearest validated | singletons | singleton validation rate vs cohort |
|---|---:|---:|---:|---|
| GILGFVFTL | 24 | 169 | 34% | 0.31 vs 0.66 |
| YLQPRTFLL | 42 | 198 | 58% | 0.10 vs 0.41 |
| NLVPMVATV | 24 | 212 | 70% | 0.13 vs 0.29 |
| TTDPSFLGRY | 108 | 183 | 83% | 0.39 vs 0.50 |
| GLCTLVAML | 24 | 144 | 37% | 0.48 vs 0.66 |

Clusters of ≥3 members have median label purity 0.89–1.00.

**Conclusion.** Validated TCRs sit in tight, label-pure clusters; non-validating
ones are isolated. **Every method in this competition is really answering "is
this TCR part of a convergent, antigen-driven cluster".** GLCTLVAML is the weak
epitope everywhere because being a singleton barely predicts anything there
(0.48 vs a 0.66 cohort rate), and density scores its singletons at AUC 0.482 —
no information at all.

## E6 — Does the representation matter?

**Goal.** Find out whether a better metric space buys anything.

**Method.** Compute the same density reductions (kernel at several bandwidths,
top-k mean, neighbour counts) over tcrdist3 and seven SCEPTR variants.
`scripts/03`, `scripts/05`.

**Result.**

| representation | best reduction | mean |
|---|---|---:|
| SCEPTR `default` | kernel σ=0.1 | **0.8396** |
| SCEPTR `small` | kernel σ=0.1 | 0.8354 |
| tcrdist3 α+β | kernel σ=20 | 0.8329 |
| SCEPTR `large` | kernel σ=0.1 | 0.8260 |
| SCEPTR `a_sceptr` (α only) | top-11 | 0.7986 |
| SCEPTR `cdr3_only` | kernel σ=0.1 | 0.8007 |
| SCEPTR `b_sceptr` (β only) | kernel σ=0.2 | 0.7505 |

**Conclusion.** A learned paired-chain representation beats a hand-built distance
by 0.007, inside one standard error. **The metric space does not matter.** Two
durable facts fell out: the **α chain carries more of this signal than β** (0.799
vs 0.751), inverting the usual expectation, and α-only is the only space that
helps GLCTLVAML.

> **Caveat added by E27.** This table computed density using only the *labelled*
> rows as the reference set, where inference uses the whole cohort. On the full
> cohort top-8 reads 0.8414 against kernel's 0.8396, reversing this table's
> ordering. The conclusion is unaffected and no submission was touched, but the
> ordering here should be read with that caveat.

## E9 — Cluster shape: real motif family or co-annotated artifact?

**Goal.** Separate convergent recombination from a knot of near-duplicates that
entered VDJdb together.

**Method.** Prompted by a diagnostic: NLVPMVATV's 23 high-density
non-validating TCRs cluster tightly with *each other* (median cosine 0.148) but
sit far from that epitope's validated TCRs (0.312), whereas elsewhere the
high-density mistakes sit right on the positives (GILGFVFTL 0.050, GLCTLVAML
0.088). So: measure cluster size, internal tightness, V/J purity, per-position
CDR3 entropy, centrality.

**Result.** Restricted to rows inside a cluster, per-position CDR3β entropy
scores **0.773** against density's 0.768, and internal tightness scores
**0.275** — meaning **tighter clusters validate less**. Imputed across the full
cohort, entropy collapses to 0.5553.

**Conclusion.** A real signal with a clean mechanism — convergent recombination
produces loose, diverse motif families, while a co-annotated knot is tight, and
plain density scores both as "dense" — but it is largely a restatement of cluster
membership and adds nothing on top.

## E20 — Four questions about what density is actually doing

**Goal.** Interrogate the mechanism rather than the score.

**(a) Does the method lean on density more than on the labels?** It uses **no
labels at all.** Density and structure are both unsupervised; all 1,782 labels
serve only to validate. That is a strength — there is nothing to overfit, which
is why the local estimate reproduced the leaderboard exactly — and an unresolved
gap: 1,782 labels exist and the predictor extracts nothing from them. Every
supervised route tried (E10, E13, E16) lost to variance or scored below chance.

**(b) Are clustered TCRs merely label-pure, or genuinely positive-leaning?**
Genuinely positive-leaning. Of 45 clusters with ≥3 members (SCEPTR cosine ≤0.12):

| cluster validation rate | clusters | TCRs |
|---|---:|---:|
| 100% (pure positive) | 26 | 205 |
| 80–100% | 8 | 136 |
| 60–80% | 8 | 31 |
| 20–60% | **0** | **0** |
| 0–20% | 0 | 0 |
| 0% (pure negative) | 3 | 13 |

The 20–60% band is **empty**, and only 13 of 385 clustered TCRs (3.4%) sit in a
pure-negative cluster. Cluster membership is positive evidence in its own right,
not a statement about homogeneity. The ceiling is equally clear: only 487 of
1,782 TCRs are in any cluster, and **about 54% of positives are singletons**.

**(c) Is anything available for singletons?** Density remains the best thing
there — mean AUC 0.758 within singletons, against structure 0.722, multi-sample
confidence 0.729, consistency 0.692, Pgen 0.497. One exception, and it is the
weak epitope: on GLCTLVAML's 81 singletons density collapses to **0.551** while
every structure-derived feature beats it (confidence 0.676, consistency 0.668,
min-PAE 0.659).

**(d) Near-neighbour label discordance.** Some negatives are single-residue
variants of validated positives, sitting inside otherwise-pure positive clusters:
34 of 872 negatives (**3.9%**), and relabelling them would move 0.8362 → 0.8602.

**Conclusion, with two cautions raised in review and both correct.** The
relabelling arithmetic is **selection-affected** — these rows were picked
*because* they sit close to positives, the same property that makes the model
rank them high. And the discordance has two possible causes this analysis cannot
separate: assay dropout, or genuine fine specificity. So 0.024 is not a
demonstrated ceiling. The honest description is **near-positive discordant
clonotypes**: negatives indistinguishable from validated positives by sequence,
of unresolved cause. E22 later relocated their damage and E31 settled the cause.

**(e) Why is GLCTLVAML weak, and does that threaten the private side?** Density
fails there because **everything is dense** — 55.2% of its cohort has a
neighbour inside the clustering radius, against 2.5–14.6% on the private three.
None of the private epitopes share this pathology. *(E21 later corrected the
inference drawn from this: cohort geometry does not forecast AUC.)*

## E32 (card D10) — What the model actually is

**Goal.** State plainly what the predictor reduces to.

**Method.** Bin by nearest-neighbour distance and score within each bin.

**Result.**

| nn1 band | n | validation rate | **within-band AUC** |
|---|---:|---:|---:|
| <0.06 | 310 | 0.823 | 0.791 |
| 0.06–0.12 | 189 | 0.921 | **0.553** |
| 0.12–0.25 | 221 | 0.805 | **0.523** |
| 0.25–0.5 | 526 | 0.432 | 0.690 |
| >0.5 | 534 | 0.139 | 0.659 |

**Conclusion.** Within-band discrimination is poor everywhere (0.52–0.79) while
the validation rate falls 0.92 → 0.14 across bands. **The model is essentially a
monotone function of nearest-neighbour distance and barely discriminates inside
any narrow band of it.** And the private epitopes sit in the far bands — 86.8% of
CINGVCWTV, 75.1% of LLWNGPMAV, 80.0% of LTDEMIAQY beyond 0.25 — which is an
independent, simpler confirmation of E30's reprojection.

---

# Rejected sequence approaches

## E7 — Generation probability

**Goal.** Test a neighbour-free alternative: are spurious annotations just
typical repertoire draws?

**Method.** OLGA Pgen per chain and paired. `scripts/04_pgen.py`.

**Result.** `-log10 Pgen(paired)` scores **0.417** — i.e. *higher* Pgen means
more likely to validate, the opposite of the hypothesis. Not orthogonal either:
Spearman with density reaches −0.63 on GILGFVFTL. Within singletons, 0.50–0.53.

**Conclusion.** Dense convergent clusters are built from easy-to-generate
sequences, so Pgen mostly re-encodes density with an unstable sign. Rejected.

## E8 — Background normalisation (TCRNET/ALICE-style)

**Goal.** Separate antigen-driven convergence from "this region of TCR space is
just crowded".

**Method.** Divide cohort density by density against the other epitopes' TCRs.

**Result.** Own-cohort density 0.8396; difference 0.8357; log ratio 0.8262.
Background density **alone** scores 0.60 *positive*.

**Conclusion.** Crowded regions of TCR space are where genuinely public
clonotypes live, so normalising subtracts signal rather than noise. Rejected.

## E10 — Fitted fusion

**Goal.** Ask whether a learned combiner beats a fixed rank average.

**Method.** Logistic regression and gradient boosting over all density and
cluster-shape features, LOEO. `scripts/07_loeo_fusion.py`.

**Result.** Best single unsupervised feature 0.8362; logistic regression on
density only 0.8359; density + cluster shape 0.8299; depth-2 booster 0.8188;
cluster shape alone 0.6846.

**Conclusion.** Five folds cannot support a fitted combiner. **Hard constraint
adopted: any fusion is a zero-parameter rank average, and a component enters only
if it adds > +0.02 in LOEO.** This threshold governed every later decision.

## E11 — Graph refinement

**Goal.** Test whether global graph structure beats local degree.

**Method.** Replace degree centrality (plain density) with eigenvector
centrality, PageRank, and 1–10 rounds of label propagation.

**Result.** Degree 0.8362; PageRank(0.85) 0.8329; propagation ×1 0.8297, ×3
0.8095, ×10 0.7417; eigenvector centrality 0.7486.

**Conclusion.** Monotonically worse with more propagation. **The signal is
strictly local**; eigenvector centrality collapses onto the largest cluster and
stops ranking the others. Rejected.

## E13 — Epitope-conditioned specificity model

**Goal.** Predict binding from TCR and peptide jointly, generalising to an
unseen epitope.

**Method.** LightGBM over SCEPTR TCR embeddings concatenated with a BLOSUM62
peptide encoding, trained on 21,112 VDJdb paired pairs across 119 epitopes;
negatives are the same TCRs paired with wrong epitopes so TCR identity carries no
label. **The target epitope is removed from training entirely** — it must be,
since the test TCRs are themselves VDJdb entries annotated to it.
`scripts/11_specificity_model.py`.

**Result.** NLVPMVATV 0.553, GILGFVFTL 0.533, GLCTLVAML 0.531, TTDPSFLGRY 0.507,
YLQPRTFLL **0.283**. Mean **0.481 — below chance**. Rank-averaged with density it
drags 0.8415 down to 0.7293. Cost: 78 min wall, 1,153 CPU-min.

**Conclusion.** Predicting binding to a peptide never seen in training remains
unsolved. What the model learns is *how generic a TCR looks against VDJdb*, and
on YLQPRTFLL the generic-looking TCRs are precisely the spurious annotations —
the same failure shape as `vdjdb_score` in E4. Two independent failures with one
mechanism: **a quantity whose sign is not stable across epitopes is worse than no
feature.**

## E16 — Learning the metric that density lives in

**Goal.** Replace a hand-chosen bandwidth and isotropic distance with learned
ones.

**Method.** Learn the kernel bandwidth alone (1 parameter) and a diagonal metric
over the 64 SCEPTR dimensions (64 parameters, shrunk toward identity) by gradient
descent on a within-epitope pairwise ranking surrogate. Strictly nested.
`scripts/14_metric_learning.py`.

**Result.** Plain cosine σ=0.1 **0.8396**; scale-only **0.8199**; diagonal
**0.8144**. Learned bandwidths came out at 0.20–0.26 in every fold, against a
true optimum of 0.07.

**Diagnosis, in three steps, because the first two were wrong.** First guess,
scale coupling: standardising density within epitope before the loss did not help
— σ still converged to 0.20–0.26. Second guess, loss shape: a bounded sigmoid was
*worse*, and both surrogates are **negatively correlated with AUC** across a σ
grid (−0.39 and −0.51). What is actually going on, from a temperature sweep:

| temperature | loss argmin (σ) | corr(−loss, AUC) |
|---|---|---:|
| 1.0 | 0.30 | −0.51 |
| 0.2 | 0.20 | −0.37 |
| 0.05 | 0.15 | −0.12 |
| 0.02 | 0.10 | **+0.03** |

The surrogate only aligns with AUC as it stops being smooth — as temperature → 0
it becomes exactly the 1−AUC step function, by which point gradients have
vanished.

**Conclusion.** The score distribution is **heavy-tailed**: a minority of TCRs
sit in dense clusters at extreme density values. AUC gives full credit for
ranking those at the very top; any smooth margin surrogate saturates there and is
dominated by the ambiguous bulk. So the surrogate is least sensitive exactly
where the metric gets most of its signal. This is a **structural mismatch, not an
implementation flaw.** What the sweep did settle: per-epitope optimal bandwidth
varies five-fold (0.03 to 0.15).

## E17 — Label-free adaptive bandwidth

**Goal.** Exploit E16's real bandwidth heterogeneity without using labels, so it
transfers to private epitopes unchanged.

**Method.** Set σ from each cohort's own distance geometry. Ten rules tested.

**Result.**

| rule | mean |
|---|---:|
| fixed 0.07 | **0.8397** |
| fixed 0.10 | 0.8396 |
| median pairwise / 12 | 0.8385 |
| median nn1 | 0.8054 |
| q10 pairwise | 0.7927 |

**Conclusion.** The rules that *work* barely adapt (`median pairwise/8` picks
0.091–0.115 everywhere — fixed σ in disguise); the rules that genuinely vary
adapt in the **wrong direction**. Mechanism: singleton-dominated cohorts have
large nearest-neighbour distances and need *small* σ so the few genuinely close
pairs dominate, while any bandwidth-scaling heuristic gives sparse cohorts *large*
σ and washes the signal out. Bandwidth heterogeneity is real but **not predictable
from cohort geometry.** Fixed σ≈0.07–0.10 stands.

## E32 — Six remaining sequence cards, all rejected

**Goal.** Close out the plan's remaining sequence ideas.
`scripts/27_remaining_sequence_cards.py`, all judged by strict nested LOEO
against the submitted blend at 0.8498.

| card | result | verdict |
|---|---|---|
| D04 substitution cost on CDR3 | identity 0.8007, BLOSUM62 0.7999, BLOSUM45 0.8014 | the cost matrix does not matter |
| D05 multi-scale neighbourhood | three-scale average 0.8378 vs single σ=0.1 at 0.8396 | worse than one scale |
| S01 residual learnable information | best extra feature +0.0002, `iptm` **−0.0193** | third independent confirmation that the labels add nothing |
| S02 singleton rescue model | global **−0.0186**; within-singleton 0.7704 vs the frozen score's 0.793 | a fitted singleton model is worse than the frozen score even on singletons |
| S07 abstention layer | close-neighbour group AUC 0.720, isolated group 0.778 | distance-based confidence does not map to ranking reliability |

**Conclusion.** S07 is worth stating in the inverse: the group with close
neighbours has a *lower* AUC than the isolated group, because its validation rate
is already 0.86 and there is little left to order. **A distance-based confidence
layer would flag where the prior is extreme, not where the ranking is
unreliable.** D10's result appears above under "What the signal actually is".

---

# Candidate-pool expansion

## E12 — Can the density estimate be sharpened with unscreened VDJdb TCRs?

**Goal.** Decide whether to pad the reference set with same-epitope TCRs not in
the test cohort.

**Method.** Count available extras per epitope.

**Result.** Public: 3,005 extra for GILGFVFTL, 9,881 for NLVPMVATV (28×), 459 for
GLCTLVAML. Private: **11 for CINGVCWTV, 5 for LTDEMIAQY**, 191 for LLWNGPMAV.

**Conclusion at the time.** Rejected on design grounds **without being run** —
enrichment would inflate public scores while leaving two of three deciding
epitopes untouched. *E25 later identified this as a procedural gap, and E27 ran
it properly; the verdict survived but the reasoning was incomplete.*

## E27 — Pool expansion, audited and quantified

**Goal.** Replace E12's assertion with measurement.

**Method and results, in three parts.** `scripts/23_pool_selfhit_audit.py`.

**(a) The query finds itself.** Every test TCR came from VDJdb, so expanding with
VDJdb risks counting the query as its own neighbour. 96.4% of queries find an
exact self-hit, and the count is **not constant** — GILGFVFTL averages 2.19, so
repeatedly curated clonotypes collect extra zero-distance "neighbours". Curation
multiplicity tracks validation (1 record → 49.0% validation, 4+ records →
**93.5%**), making it a label proxy masquerading as density. Confirmed directly:
raw expansion **loses** 0.0024 of mean AUC; query-removed expansion gains 0.0048.

**(b) The aggregation operator decides whether expansion helps at all.**

| operator | no pool | pooled | delta |
|---|---:|---:|---:|
| top-8 | 0.8414 | **0.8526** | **+0.0112** |
| top-3 | 0.8334 | 0.8476 | +0.0142 |
| top-16 | 0.8377 | 0.8488 | +0.0111 |
| kernel σ=0.1 | 0.8396 | 0.8432 | +0.0036 |

NLVPMVATV, whose pool is 94.8% a single dextramer study, goes **−0.0354 under
kernel-mean and +0.0131 under top-8**: a mean over the whole pool is drowned by
one study's records, while top-k never looks at them.

**(c) What the private side can realise: +0.0013 on top-8, +0.0004 on the
submitted blend.** Downsampling each public pool to the private epitopes' own
distinct-clonotype multipliers, 30 repeats: **88% of the public gain cannot
transfer**, because two of three private epitopes can add only ~4% more distinct
clonotypes.

**Three reasons not to quote the point value alone** — all three found by
self-review after the first estimate was already written down:

1. It was first computed on a model not submitted. On the submitted
   configuration the macro estimate falls from +0.0013 to **+0.0004**, naive 95%
   interval [−0.0022, +0.0030], negative on two of three private epitopes.
2. The mean hides a **26-fold spread** (GLCTLVAML +0.0158 to NLVPMVATV −0.0076,
   sd five times the mean). Averaging across public epitopes assumes the gain
   depends only on the multiplier — exactly the exchangeability E21 disproved.
   **Five epitopes cannot estimate this quantity.**
3. Random subsampling is optimistic. LTDEMIAQY's 5 available extras come from a
   pool that is 97.8% one reference; drawing only from each pool's dominant study
   gives **−0.0002**.

**Conclusion.** Not "the gain is small" but **"on our configuration there is no
measurable gain"**. E12's verdict stands; its reasoning was incomplete — the
problem is not only that few extra TCRs exist, but that the ones that do are
duplicates from single studies.

## E33 — Pool and provenance cards, with two results that correct each other

**Goal.** Test whether the pool gain is evidence or volume.

**(a) I02 — expansion rescues positives preferentially.** When the pool gives an
isolated TCR a neighbour inside the clustering radius, that TCR is about **three
times** more likely to be a validated positive (GILGFVFTL 3.2×, pooled ~3×). So
expansion carries genuine evidence — but outside GILGFVFTL the absolute counts
are negligible (two TCRs rescued each on three epitopes).

**(b) R02 — and yet throttling kills the gain.**

| pool policy | mean |
|---|---:|
| no pool | 0.8414 |
| all records | **0.8526** |
| one share per reference | 0.8504 |
| one per subject | 0.8428 |

Deduplicating by subject nearly erases the benefit.

**(c) R08 — how much "clustered" just means "easy to generate" is
epitope-specific.** Matching on Pgen terciles × CDR3β length, the
clustering-to-validation gap goes from 0.400 → **0.067** on GILGFVFTL but
0.658 → **0.711** on NLVPMVATV.

**(d) S05 — provenance weak supervision, dead for the third time.** Records per
clonotype 0.534, distinct references 0.557, assay classes 0.514 — and each takes
only one or two distinct values within every private epitope.

**Conclusion.** The gain has both an **evidential** component (I02) and a
**volume** component (R02), which is why raw-record expansion looked better than
it deserved. And the degree of Pgen confounding is itself epitope-specific, which
is one more reason a single global correction (E8) could only hurt.

## E35 — Reference-database version, source and leakage

**Goal.** Check whether the pool conclusions depend on which VDJdb snapshot was
used, and whether the pool leaks labels.

**Method.** Three releases — 2024-06-13 (32,183 paired human rows), 2025-02-21
(35,672), 2026-06-03 (86,221) — with the pool measurement rerun on each.

**(a) Leakage audit, which had not been done and should have been.** Searching
every human reference for the Immunity 2026 DOI, the 2025 bioRxiv preprint,
"TCRvdb", the Schumacher lab repo and "T-RAP" returns **zero records each**. The
pool is independent of the labels.

**(b) The screening set most resembles a ~2025-02 snapshot.** Median absolute
difference between VDJdb clonotype count and test-set size is 4 for 2025-02-21,
against 10 and 25 for the others; CINGVCWTV matches exactly, 243 against 243.

**(c) The pool gain is not resolvably version-dependent.** Gains of +0.0052,
+0.0099, +0.0112 across the three releases look monotone, but **no pairwise
difference reaches two standard errors**, and the spread is carried almost
entirely by GLCTLVAML. A 2.7-fold increase in database size produces no
resolvable change — which is itself the result. This also falsifies a claim made
in narration before measuring.

**(d) The durable finding is a count, not an AUC, so noise does not apply.**

| epitope | 2024-06 | 2025-02 | 2026-06 |
|---|---:|---:|---:|
| CINGVCWTV | 1.01 | 1.00 | 1.05 |
| LTDEMIAQY | 1.04 | 1.02 | 1.04 |
| LLWNGPMAV | 1.08 | 1.03 | **1.83** |

**Conclusion.** CINGVCWTV and LTDEMIAQY gained essentially no new paired TCRs in
eighteen months, so **their exhaustion is a durable property of the epitope**,
making E27's and E33's conclusion version-robust for two of the three deciding
epitopes.

## E38 — Negative controls and the pool × structure factorial

**Goal.** Both cards were originally skipped on my own judgement. R04 asks the
sceptical question nothing else covered: does adding *any* sizeable batch of
sequences raise density regardless of what they are?

**Method.** Three size-matched controls, all adding 2,684 sequences — the same
count as the real pool: real TCRs from other epitopes matched on CDR3β length and
TRBV; the same drawn unmatched; and OLGA-generated synthetic repertoire.
`scripts/33_pool_negative_control.py`.

**Result.**

| reference added | mean | gain |
|---|---:|---:|
| none | 0.8414 | — |
| **real same-epitope** | **0.8526** | **+0.0112** |
| other-epitope, shape-matched | 0.8333 | **−0.0081** |
| other-epitope, unmatched | 0.8350 | −0.0064 |
| OLGA synthetic repertoire | 0.8381 | −0.0063 |

**Conclusion.** **All three controls actively hurt.** Adding the same number of
irrelevant sequences does not merely fail to help, it degrades the density
estimate. So the expansion gain is **antigen-specific information**, which
together with I02 and R02 completes the picture: a real mechanism, carried by
multiple records from the same subjects, that the private epitopes cannot access.
R05's factorial was degenerate as expected, confirming with numbers what had been
asserted.

---

# Structure

## E14 — Structure prediction at scale

**Goal.** Test whether a predicted interface supplies signal orthogonal to
sequence neighbourhood.

**Method.** Fold all 2,636 TCR–peptide–MHC complexes with Boltz-2 (four chains:
TCRα V, TCRβ V, allele-specific HLA α1α2, peptide; ~416 residues;
single-sequence; 14 s each, 7.4 h total). Extract interface PAE and chain-pair
ipTM. `scripts/08`–`10`.

**Result.** Best standalone structure feature 0.759, below the organisers' AF3
baseline. The informative part is *which* features work — all of them, equally,
including ones that cannot carry specificity:

| feature | mean |
|---|---:|
| `pae_a_mhc_min` (TCRα–MHC) | 0.759 |
| `pae_tcr_pep_min` (organisers' statistic) | 0.759 |
| `complex_plddt` (global fold quality) | 0.694 |
| `pae_ab_min` (TCRα–TCRβ) | 0.656 |
| `iptm_mhc_pep` (**MHC–peptide**) | 0.647 |

`iptm_mhc_pep` is the diagnostic: that interface contains **no TCR** and is
identical for every complex in a cohort, yet scores 0.647.

**Conclusion.** Structure confidence here largely reports **how canonical the TCR
is** — typical sequences fold confidently everywhere through the shared pair
representation — so it is not the orthogonal signal hoped for. It still adds a
little: blended at weight 0.2–0.3, nested LOEO gives **+0.0070** (0.8392 →
0.8461), with the weight stably selected in every fold. Four of five epitopes
improve; the sole regression is TTDPSFLGRY, **the only 10-mer**, and all three
private epitopes are 9-mers, so restricted to the four 9-mers the gain is
+0.0115.

## E19 — Docking consistency across diffusion samples

**Goal.** Find a structural quantity immune to the canonicality confound.

**Method.** Re-fold all 2,636 complexes with 5 diffusion samples (13.2 GPU-hours).
Superpose each pair of samples on the MHC and peptide Cα, then take TCR Cα RMSD.
`scripts/16`, `scripts/17`.

**The hypothesis, and why it was wrong.** Consistency was expected to be immune to
the confound because all five samples share identical sequence input. **That
argument was wrong.** The five samples share more than the input: they share the
single already-computed trunk representation. The diffusion head samples a
distribution conditioned on it, so a confident trunk yields a narrow conditional
and therefore consistent samples. **Consistency is downstream of trunk
confidence.**

**Result.** Docking RMSD spans 0.32–42.7 Å, so the dynamic range is real, but
within-epitope Spearman against density is **+0.530** — indistinguishable from the
+0.528 of the confidence features it was meant to improve on. Blended, the best
gain is **+0.0010** against a +0.02 threshold.

**Conclusion.** 13.2 GPU-hours for +0.001. **Rejected.** The pre-registered rule
did its job. One thing worth keeping: consistency scores 0.775 on GLCTLVAML, the
highest of any single feature on the project's weakest epitope — but weighting it
heavily enough to exploit that drops the overall mean to 0.8128.

## E28 — Contact fingerprints and docking geometry

**Goal.** Move from confidence scores to actual interface geometry.

**Method.** All four cards extracted in one pass over the 2,636 structures,
`scripts/24_structure_contacts.py`. Peptide position is the common coordinate —
the peptide is identical within an epitope — so no IMGT numbering is needed for
the one thing that must be exactly right.

**Result.** **Geometry QC passes:** peptide–MHC contacts median 334, minimum 260,
so every peptide sits in the groove; CDR3 supplies 74–76% of all TCR–peptide
contacts, the biologically expected share. One complex,
`1_6_B11_22_YLQPRTFLL`, has **no TCR–peptide contact at all** — a genuine docking
failure in a public epitope.

Contact-fingerprint density scores **0.7476**, correlated +0.563 with density —
another re-encoding. **But the raw contact counts are the most orthogonal
quantity measured in this project:** within-epitope Spearman against density is
**+0.088** and +0.103, where everything else ever measured sat between +0.26 and
+0.91. They are also too weak to use — standalone AUC 0.588 and 0.559, nested
gain **+0.0008**.

**Conclusion, and a general lesson worth keeping: orthogonality is necessary but
not sufficient.** A component at AUC 0.59 blended into a model at 0.85
contributes almost nothing however independent it is — the achievable
contribution scales with its margin over chance, and 0.09 does not buy anything.

## E34 — Predicted contacts against crystallography

**Goal.** E28 and E31 both rest on predicted interfaces. Check them against
experiment.

**Method.** The PDB yields 14 solved TCR–pMHC complexes for four competition
epitopes, folded through the *same* pipeline used for the competition.
`scripts/28_real_structure_validation.py`, `scripts/29_fold_real.sh`.

**Result: the fingerprints are faithful.** Median Spearman on contact counts per
peptide position 0.992; on minimum distance 0.983; Jaccard of touched positions
1.000; total contacts 102 predicted vs 99 measured.

**Three self-checks, two of which cut the claim down.**

- *Zero inflation.* Of 20 fingerprint slots, 11.5 on average are zero in both
  vectors — buried anchor positions no TCR can reach. Restricting to
  contact-bearing positions drops the median from 0.992 to **0.953**. The
  headline number was partly agreement on what is structurally impossible.
- *Memorisation.* All 14 complexes were released 2003–2022, so all predate
  Boltz-2. **This is a memorisation check, not a prediction test.** The one
  post-training complex available (9x7r, released 2026-09-23) drops to Spearman
  **0.840**, Jaccard **0.700**, with contacts overpredicted by 40%. n=1 is an
  anecdote, but it is the only out-of-distribution point and it points the
  expected way.
- *Between-complex discrimination — the check that matters most.* Over the 36
  pairs of GILGFVFTL complexes, Spearman between predicted and measured pairwise
  fingerprint distance is **+0.948**. Within the memorised regime the
  representation does carry real inter-TCR interface variation.

**Conclusion, and the consequences are asymmetric.** E28's negative result is
**safe and strengthened** — prediction noise can only weaken a real signal.
E31's positive result is **exposed but partly defended by its design**: it
compared discordant pairs against matched same-label controls from the same
pipeline, so a uniform prediction bias affects both arms equally. Producing that
effect from prediction error alone would require the error to correlate with the
mutation itself. A real defence, not immunity.

## E36 — Target against anchor-matched decoy peptides

**Goal.** Test E14's inference directly: is the structural signal peptide-specific
at all?

**Method.** Fold the same TCR against its annotated peptide and against decoys,
then subtract — cancelling everything that is a property of the TCR alone.
**Decoys are anchor-matched by class, not identity**: HLA-A\*02:01 constrains an
aliphatic P2 and C-terminus, A\*01:01 an acidic P3 and aromatic C-terminus. A
first draft asserted identical residues and failed immediately — no reported
A\*02:01 nonamer shares GILGFVFTL's P2 isoleucine *and* its C-terminal leucine —
which also exposed two decoys chosen with a serine at P3.
`scripts/30`, `scripts/31`.

**Result: the margin carries nothing.**

| quantity | target only | target minus decoy |
|---|---:|---:|
| ipTM | **0.7387** | 0.5362 |
| confidence score | 0.7189 | 0.5398 |
| TCRα–peptide ipTM | **0.7429** | 0.5331 |
| TCRβ–peptide ipTM | 0.7329 | 0.5251 |
| MHC–peptide ipTM | 0.6500 | 0.5147 |
| complex pLDDT | 0.6313 | 0.5247 |

Mean target ipTM 0.9152 against decoy 0.9091 — a margin of **+0.0061**. Boltz
assigns a TCR's true target and an unrelated anchor-matched epitope essentially
the same interface confidence.

**Four self-checks, which turn this from "noisy" into "far too small".** The
decoy structures *are* different (spread across decoys 0.036, comparable to 0.062
across TCRs). The margin's mean is a tenth of its between-TCR spread, and
decoy-choice noise alone is **5.9×** the margin. No single decoy carries it —
per-decoy margins run −0.0128 to +0.0247, two negative. And reaching AUC 0.70
would take a margin near 0.042, so the observed value is **14% of what is
needed**: not a power problem more decoys would fix.

**Phase 2 was run after all, and it was worth seven GPU-hours.** I had staged the
private extension behind a positive phase-1 result and then skipped it, arguing
that a feature failing on all five labelled epitopes can never clear the
threshold. The argument stands but answers a different question from the one
phase 2 asks: *is the failure universal, or specific to the epitopes that happen
to carry labels?* Mean absolute margin is 0.0406 on labelled epitopes and 0.0394
on private ones — indistinguishable. **The failure is universal.** Two public
epitopes carry **negative** mean margins: on TTDPSFLGRY and YLQPRTFLL the model
prefers the decoys on average.

**Conclusion.** This closes the structural line with its mechanism identified:
**Boltz interface confidence contains no peptide-specific information in this
setting**, so every structural feature tried (E14, E19, E28, E36) was reading the
same TCR-canonicality signal that density already reads.

## E37 — Does the E31 interface result survive relaxation?

**Goal.** E31 is the project's one positive mechanistic finding and it was
measured on raw Boltz coordinates, which carry clashes and unrelaxed side chains.

**Method.** 69 complexes relaxed under one fixed protocol — PDBFixer for terminal
atoms and hydrogens at pH 7, Amber14 with GBn2 implicit solvent, L-BFGS
minimisation at 200 and 1,000 iterations.
`scripts/32_minimisation_robustness.py`.

**Result — descriptor stability differs sharply**, which decides which
descriptors a conclusion may rest on. *(Numbers corrected by the audit; see
below.)*

| descriptor | median relative change, 0→1,000 steps | Spearman(raw, relaxed) |
|---|---:|---:|
| buried interface area | **1.63%** | +0.960 |
| peptide contacts | 9.21% | +0.793 |
| hydrogen bonds | **33.33%** | +0.788 |

Buried area is robust; hydrogen-bond counts are not, so a conclusion resting on
them would not be safe. **The E31 effect survives**: contacts at the mutated site
read disc 8.61 vs ctl 3.97 (p=0.026) raw and 7.16 vs 3.53 (p=0.045) minimised,
with the contact ratio going 2.17× → 2.03×.

**Conclusion.** E31 is not an artefact of unrelaxed coordinates. But see the two
revisions below — this entry both overstated its evidence and misreported three
of its own numbers.

### E37 revision — the three rows are one finding, and it is positional

Found while preparing E39, which needed the same mutated-site residues. Measured
on E37's own relaxed coordinates, 40 pairs and 126 site measurements,
`scripts/35_md_analysis.py --static`.

Reading a **second descriptor** of the same physical quantity splits the result
apart: the contact count separates the groups (3.581 vs 1.766, p=0.017) but the
minimum site-to-peptide distance does not (6.412 vs 6.354, p=0.899).

**Not an outlier artefact — and the first version of this audit checked that
wrongly.** It rank-transformed the count and re-ran, which returns the identical
p-value by construction, because Mann-Whitney is already a rank test. The real
reason: **77% of mutated sites have no peptide contact at all**, so the count is
three quarters ties at zero and its rank test is in effect a binary comparison,
while the distance uses all 126 sites.

Testing that binary reading directly, and then conditioning on it:

| quantity | disc | ctl | p |
|---|---:|---:|---:|
| **site touches the peptide at all** | **0.323** | **0.141** | **0.020** (Fisher) |
| distance given the site touches (Å) | 3.504 | 3.319 | 0.195 |
| contacts given the site touches | 11.100 | 12.556 | 0.981 |

**The whole effect is the binary one.** Discordant pairs differ at a
peptide-facing position 2.29× as often — numerically E31's headline "2.3× more
peptide contact". Conditional on facing there is no difference at all.

**Two consequences.** **E31 should be restated**: not "discordant negatives
mutate at positions carrying more peptide contact", which implies a graded
energetic difference, but "discordant pairs are about twice as likely to differ
at a peptide-facing position". The mechanism is **positional, not energetic**.
And **E37 overstated its own evidence** — its three rows are not three
independent confirmations; with 77% ties the count row and the facing row are
substantially the same test.

### E37 correction — three numbers in its descriptor table are not reproducible

Found by the data audit, which flagged that `scripts/32` declared
`reports/10_minimisation_robustness.md` and never wrote it, so E37's stability
table existed only in the log.

| descriptor | logged | reproducible | logged ρ | reproducible ρ |
|---|---:|---:|---:|---:|
| buried interface area | 2.1% | **1.63%** | +0.963 | +0.960 |
| peptide contacts | 8.5% | **9.21%** | +0.793 | +0.793 |
| hydrogen bonds | 33.3% | 33.33% | +0.748 | **+0.788** |

The contacts ρ matches exactly and the hydrogen-bond median matches exactly,
which fixes the stage pair and denominator; **no variant** of stage pair,
denominator, mean vs median, or signed vs absolute reproduces 8.5% or 2.1%. The
likely cause is the recurring one: the table was computed in a heredoc that was
never saved, probably against a partially written CSV. The conclusion is
unchanged and if anything strengthened, but **a sixth of E37's reported figures
were wrong in the second decimal because the computation was not in version
control.**

## E39 — Is the interface difference dynamically stable?

**Goal.** The last structural card. E31 and E37 are both static; function depends
on how long a contact persists, not only on whether it exists in one relaxed
snapshot.

**Method.** 69 relaxed complexes from 40 pairs, three velocity seeds each, 2 ns at
310 K with a 4 fs step and hydrogen-mass repartitioning, implicit GBn2, sampled
every 20 ps. **207 trajectories, 414 ns, no failures.**
`scripts/34_md_interface_stability.py`, `scripts/35_md_analysis.py --md`.

**The setup holds the pose.** Interface contact rises a median 7.1% over each
trajectory and 20 of 207 trajectories lose more than a quarter of their contacts.
Side chains settle into contact on heating from the minimised state; nothing is
dissociating.

**Every group comparison at the mutated site is null.**

| quantity | disc | ctl | p |
|---|---:|---:|---:|
| contacts, time-averaged | 2.456 | 1.332 | 0.159 |
| distance, time-averaged (Å) | 6.817 | 6.803 | 0.545 |
| fraction of frames in contact | 0.288 | 0.154 | 0.149 |
| **touches the peptide at any point** | **0.597** | **0.594** | **1.000** |

The last row is the one to read against the static result: on the relaxed
snapshot the binary facing rate was 0.323 vs 0.141 at p=0.020; given 2 ns and
three replicas, both groups reach the peptide at the same rate.

**MD adds nothing beyond the static facing rate, and that is an identity, not an
impression.** Conditioning on whether the site faced the peptide in the relaxed
pose, both groups behave the same: among the 29 facing sites the contact fraction
is 0.758 vs 0.708 (p=0.943), and among the 97 non-facing sites 0.064 vs 0.063
(p=0.579). Decomposing:

```
discordant  0.323 × 0.758  +  0.677 × 0.064  =  0.288   (observed 0.288)
control     0.141 × 0.708  +  0.859 × 0.063  =  0.154   (observed 0.154)
```

The conditional terms are group-independent to three decimals. The dynamic
contact fraction **is** the static facing rate carried through a persistence
factor that does not distinguish the groups.

**The direction is not consistent across epitopes** — the third pre-planned
self-check, and it fails: GILGFVFTL 1.161 vs 0.641 and NLVPMVATV 5.322 vs 0.187
are in the claimed direction, while GLCTLVAML 2.835 vs **9.357** and YLQPRTFLL
**0.000** vs 4.442 are reversed, on 2–8 sites in the small cells.

**Conclusion, against a limit stated before the run.** It was pre-registered that
only 29 of 126 sites face the peptide (20 discordant, 9 control), so a null result
would be **uninformative rather than evidence of absence**. That is what
happened. C11 neither supports nor refutes E31. What it does establish is
narrower and worth having: **the effect is entirely contained in the static
facing rate, so no further structural dynamics will add to it**, and this line
closes on evidence rather than on judgement.

---

# Validation methodology and error anatomy

## E15 — How much does *choosing a configuration* cost?

**Goal.** Price the selection step itself as part of the model.

**Method.** 323 candidate configurations. Compare the oracle (best config scored
on all five epitopes), nested LOEO (config re-chosen inside each fold), and a
config fixed a priori then scored by plain LOEO. `scripts/13_strict_loeo.py`.

**Result.**

| | mean |
|---|---:|
| oracle — tuned on the evaluation set | 0.8508 |
| **nested LOEO — tune by CV, then submit** | **0.8460** |
| best a priori config, plain LOEO | 0.8508 |

**Conclusion.** Cross-validated selection is **0.0048 worse than fixing the
configuration in advance**. Mechanism, visible in the TTDPSFLGRY fold: on the
other four epitopes `pae_a_mhc_min@0.3` looked best and returned 0.8328 on the
held-out epitope, where the fixed config gave 0.8504. **With five folds the
selection step has worse signal-to-noise than the thing being selected. Make
fewer decisions, not better-cross-validated ones.** Every earlier number quoted
for a *chosen* configuration is optimistic by roughly this much.

## E18 — Master cross-epitope comparison

Everything in one protocol: five public epitopes, LOEO, mean per-epitope AUC,
nothing tuned on the epitope it is scored on. `scripts/15_master_comparison.py`.

| method | mean |
|---|---:|
| FUSION: 3-space density + structure @0.2 | **0.8508** |
| FUSION: 2-space density + structure @0.2 | 0.8496 |
| FUSION: 3-space density | 0.8430 |
| FUSION: 2-space density | 0.8415 |
| density: SCEPTR default | 0.8396 |
| density: SCEPTR small | 0.8354 |
| FUSION: 3-space density + specificity @0.2 | 0.8348 |
| density: tcrdist3 kernel | 0.8329 |
| tcrdist3 nndist (organisers' baseline) | 0.8247 |
| density: SCEPTR α-chain only | 0.7986 |
| structure: TCR–peptide min-PAE | 0.7565 |
| density: SCEPTR β-chain only | 0.7505 |
| structure: whole-complex confidence | 0.7467 |
| metadata: vdjdb_score | 0.4849 |
| specificity model | 0.4812 |
| generation probability | 0.4165 |

## E21 — Cohort thinning: does density survive a sparser repertoire?

**Goal.** The private cohorts are sparser than anything validated on, so the LOEO
estimate is an extrapolation. Simulate it.

**Method.** Randomly thin each public cohort to fractions 1.0 down to 0.12,
recompute density within the thinned cohort, score against real labels, 40
repeats. Two variants, because the naive design confounds a worse density
estimate with a noisier AUC on fewer rows. `scripts/19_cohort_thinning.py`.

**Result 1 — robustness to sparsity is wildly epitope-dependent.**

| epitope | dense end | sparse end | change |
|---|---:|---:|---:|
| GILGFVFTL | 0.868 | 0.869 | **0.000** |
| YLQPRTFLL | 0.928 | 0.920 | −0.008 |
| GLCTLVAML | 0.713 | 0.688 | −0.025 |
| TTDPSFLGRY | 0.845 | 0.724 | **−0.121** |
| NLVPMVATV | 0.844 | 0.699 | **−0.145** |

**Result 2 — and this corrects E20.** At essentially identical sparsity
(close-neighbour fraction ≈0.013–0.014) YLQPRTFLL scores **0.920** and NLVPMVATV
**0.699**. **Sparsity does not determine performance; intrinsic separability
does.** So E20's claim that the private epitopes "look like TTDPSFLGRY,
NLVPMVATV and YLQPRTFLL where density scores 0.84–0.93" **does not hold** — at
private-like sparsity those same three span 0.70–0.93. Cohort geometry cannot
forecast a private epitope's AUC.

**Result 3 — thinning is not strongly pessimistic.** The singleton-to-negative
association barely moves (gap 0.486 → 0.454 at 20% coverage), so the simulation
is reasonably faithful rather than a loose bound.

**Result 4 — structure is insurance, and it changed the submission.**

| | full | 45% | 20% |
|---|---:|---:|---:|
| gain on the two fragile epitopes | +0.0059 | +0.0120 | **+0.0202** |
| gain on the three robust epitopes | +0.0112 | +0.0103 | +0.0099 |

**Decision: structure weight 0.2 → 0.3**, on a dominance argument rather than by
optimising a regime. w=0.3 ties w=0.2 at full coverage (−0.0002) and beats it at
both sparse levels; the private cohorts are the sparse case. w=0.5 wins at 20%
coverage but loses 0.013 at full coverage, so its risk is asymmetric. **If the
sparsity argument is simply wrong, the cost of 0.3 is bounded at 0.0002.**

## E22 — Pairwise AUC error decomposition

**Goal.** Mean per-epitope AUC is exactly the probability a validated TCR
outranks a non-validating one from the same cohort, so 1−AUC is a set of specific
inverted pairs. Find out which ones. The inference this guards against: 54% of
positives are singletons, which does **not** imply 54% of the error involves them.

**Method.** Partition every positive–negative pair within an epitope by cluster
state, mutually exclusive and complete, so the class losses sum to 1−AUC. The
script asserts that reconstruction; residual 1.4×10⁻¹⁶.
`scripts/20_auc_error_decomposition.py`.

**Result.**

| class | share of loss | share of pairs | error rate within class |
|---|---:|---:|---:|
| both_singleton | **62.0%** | 50.5% | 0.207 |
| pos_sing_neg_clu | **21.2%** | **4.5%** | **0.673** |
| diff_cluster | 10.0% | 6.2% | 0.231 |
| pos_clu_neg_sing | 6.4% | 38.6% | 0.021 |
| same_cluster | 0.5% | 0.3% | 0.309 |

**Conclusion.** **The adversarial class is systematically inverted**: when the
positive is isolated and the negative is clustered, the model gets the order
wrong 67% of the time — worse than chance. It is 4.5% of pairs and 21.2% of the
error. **This corrects E20's account of the 34 near-positive discordant
negatives**: same-cluster pairs contribute only 0.5% of loss, so relabelling them
did not gain 0.024 by fixing within-cluster discrimination. It gained it because a
relabelled negative stops being a negative in *every* comparison, and most of its
damage was in `pos_sing_neg_clu`. **The mechanism is over-scoring of clustered
negatives, not failure to resolve a neighbour.**

## E23 — Can the adversarial class be fixed? No, and the arithmetic says why

**Goal.** Act on E22's diagnosis.

**Result — per-class feature strength.** Structure is better than chance exactly
where density is catastrophically worse:

| feature | pos_sing_neg_clu | pos_clu_neg_sing |
|---|---:|---:|
| density | **0.697** | **0.013** |
| blend w=0.3 | 0.602 | 0.023 |
| structure min-PAE | **0.422** | 0.156 |

**Attempt — cluster-conditional weighting.** Cluster state is label-free, so
weighting structure differently for clustered and isolated TCRs is legitimate at
inference. Strict nested LOEO over both weights: **0.8418 against 0.8486 for the
flat weight, −0.0068.** And the oracle over the whole 8×8 grid peaks at 0.8486 at
(0.2, 0.3) — *identical to flat w=0.3*. Even with perfect knowledge, conditioning
buys nothing.

**Why, quantitatively.** Scores are per TCR, not per pair, so a TCR cannot know
which comparison it is in. Raising structure weight from 0.3 to 0.6 gives
**+0.0040** in the adversarial class (4.5% of pairs) and **−0.0336** elsewhere —
net −0.0296. At w=1.0 the best available weighting is worth +0.0081
while costing −0.107. *(Corrected in E40's review: the class's actual ceiling is
its mean loss of 0.0318, so eliminating its error entirely would be worth
+0.032. The +0.0081 is what the tested weighting family achieves.)* **The trade is about 13:1 against.**

**Conclusion.** The decomposition succeeded diagnostically and failed
prescriptively. What would be needed is a feature informative *within* the
clustered subset and roughly uninformative *between* subsets. Structure is not
that: it is informative everywhere, merely worse than density in most places.

## E24 — Is probabilistic label modelling available? No

**Goal.** A reviewer suggested modelling
`P(Y_true=1 | Y_observed, assay evidence)` rather than a hard binary.

**Checked.** The released TCRvdb file has 24 columns and the only label
information is `pos`/`neg`/`ambiguous`, strictly one-hot integers. There is no
enrichment value, read or UMI count, replicate, p-value, CD69 readout or any
other continuous quantity; the sole numeric field is `vdjdb_score` (0–3), which
is the *prior annotation* confidence and already known to anti-correlate with
validation (E4).

**Conclusion.** The route is **closed with the data available**, not merely
unpromising. It would need a richer release from the authors.

## E29 — Family holdout and chain discordance

**(a) S06 — leave-one-clonotype-family-out.** Each TCR scored against a cohort
with its own single-linkage family removed, so no relative can supply its density:

| epitope | in a family | standard | family-masked | delta |
|---|---:|---:|---:|---:|
| NLVPMVATV | 0.183 | 0.8436 | 0.7417 | **−0.1019** |
| GILGFVFTL | 0.476 | 0.8678 | 0.8520 | −0.0158 |
| GLCTLVAML | 0.552 | 0.7130 | 0.7000 | −0.0130 |
| TTDPSFLGRY | 0.047 | 0.8449 | 0.8356 | −0.0093 |
| YLQPRTFLL | 0.184 | 0.9285 | 0.9272 | −0.0012 |
| mean | 0.288 | 0.8396 | **0.8113** | **−0.0282** |

Family recurrence carries part of the ranking, and once again the dependence is
epitope-specific rather than structural: NLVPMVATV and YLQPRTFLL have
near-identical family rates (0.183, 0.184) and lose 0.102 against 0.001.

**(b) D03 — the α/β asymmetry is not mispairing.** The mispairing hypothesis
predicts that a mispaired TCR shows large single-chain rank disagreement.
Measured, signed disagreement scores 0.536 and absolute disagreement **0.429 —
below chance**, meaning *larger* α/β disagreement associates with validating, the
opposite of what mispairing would produce. Rejected, and the mispairing
explanation gets no support.

## E30 — Why leave-one-epitope-out overestimates the private score

**This is the most consequential result in the project.**

**Goal.** S06 prompted a comparison of cohort *structure* rather than cohort
size. It revealed something LOEO cannot see.

**Result 1 — the private epitopes have four times less family structure**, 7.2%
against 28.8% of TCRs sitting in a family of two or more:

| epitope | split | in a family |
|---|---|---:|
| CINGVCWTV | PRIVATE | 0.025 |
| LLWNGPMAV | PRIVATE | 0.044 |
| LTDEMIAQY | PRIVATE | 0.146 |
| TTDPSFLGRY | public | 0.047 |
| NLVPMVATV | public | 0.183 |
| YLQPRTFLL | public | 0.184 |
| GILGFVFTL | public | 0.476 |
| GLCTLVAML | public | 0.552 |

That matters because of E22: `pos_clu_neg_sing` is 38.6% of public pairs at a
**2% error rate**. Those easy pairs hold a public AUC up, and CINGVCWTV will have
about 5.6% of them.

**Result 2 — decomposing and reprojecting.** Taking each class's measured
within-class error rate for the submitted model and recombining with each
epitope's own class composition **reproduces every public epitope's actual AUC to
within 0.0003**, so the arithmetic is sound. Applying it to the private epitopes'
label-free compositions:

| epitope | both-singleton share | easy share | projected AUC |
|---|---:|---:|---:|
| CINGVCWTV | 92.6% | 5.6% | 0.795 |
| LLWNGPMAV | 87.2% | 9.8% | 0.797 |
| LTDEMIAQY | 60.7% | 29.2% | 0.806 |
| **macro** | | | **0.799** |

against the LOEO estimate of 0.846.

**Result 3 — the whole projection collapses to one unknown.** Since the private
epitopes are 61–93% both-singleton pairs, the projected macro tracks the
within-both-singleton AUC almost exactly:

| if private singletons separate like | both-singleton AUC | projected private macro |
|---|---:|---:|
| GLCTLVAML | 0.642 | **0.668** |
| NLVPMVATV | 0.748 | 0.758 |
| GILGFVFTL | 0.827 | 0.831 |
| TTDPSFLGRY | 0.830 | 0.832 |
| YLQPRTFLL | 0.915 | **0.907** |

**Conclusion.** **LOEO holds out epitope identity but not cohort structure.** It
is unbiased for "an epitope drawn like the public ones" and biased upward for an
epitope with private-like family sparsity. That bias is about **−0.047**, and no
amount of leave-one-epitope-out can detect it, because all five public epitopes
are family-richer than the private three.

**Revised position at the time: central estimate 0.80, plausible range
0.67–0.91.** **E40 subsequently withdrew the 0.80**: the projection is not
identifiable without label-dependent assumptions, its "reproduces every public
AUC" check is an identity, and it fails an out-of-sample transfer test. The
qualitative finding stands; the number does not. The spread is **not sampling noise** — it is one substantive
unknown, how separable the private epitopes' isolated TCRs are, and E21 already
showed that cannot be forecast from cohort geometry. This is the convergence
point of E20, E21, E22 and E32's D10.

---

# Mechanism of the irreducible error

## E31 — The discordant negatives look like fine specificity, not assay dropout

**Goal.** E20 found negatives that are single-residue CDR3 variants of validated
positives; E22 relocated their damage. Nothing had tested *why* they differ, and
the two candidates imply opposite things: a mutation abolishing a peptide contact
is real fine specificity, so no sequence method can be blamed; a mutation at a
solvent-facing position points at the assay or the annotation.

**Method.** The test needs a null, so each discordant mutation is compared
against **matched same-label pairs** — positive/positive neighbours sharing one
CDR3 exactly and differing by one or two residues on the other, within the same
epitope. Without that control, "mutations tend to be peptide-facing" would just
be a property of CDR3 mutation in general. 71 discordant pairs (103 mutated
sites) against 522 control pairs (767 sites). Permutation is **stratified by
epitope** over 2,000 draws, so an effect produced only by GILGFVFTL dominating
the pool would wash out. `scripts/26_discordant_pair_anatomy.py`.

**Result.**

| mutated-site property | discordant | control | stratified permutation p |
|---|---:|---:|---:|
| peptide contacts at that position | **6.53** | 2.90 | **0.0070** |
| change in peptide contacts between members | **3.37** | 1.27 | **0.0010** |
| position is peptide-facing at all | 0.340 | 0.231 | 0.0790 |
| volume change | 53.7 | 41.9 | (Mann-Whitney 0.023) |
| charge change | 0.27 | 0.23 | (Mann-Whitney 0.314) |

**Conclusion.** This supports genuine fine specificity over assay dropout: if
these were dropouts, the mutated position would be random with respect to the
peptide interface and therefore indistinguishable from the control. It is not. So
near-positive discordance is at least partly **real biology**, which means
sequence-neighbourhood methods have a **biologically grounded ceiling** rather
than merely a noise ceiling — and it cannot be fixed by better sequence
modelling, only by reading the interface, which E28 showed is too weak to use.

Limits stated at the time: these are predicted structures; charge change shows no
effect, so the simple "salt-bridge lost" story is not supported; and this explains
the ceiling without moving the competition score at all.

> **Restated by the E37 revision.** The surviving effect is the **binary** one —
> discordant pairs differ at a peptide-facing position 2.29× as often — and *not*
> a graded energetic difference. Conditional on the site facing the peptide there
> is no difference in distance or contact count. The mechanism is **positional**.
> E39's 414 ns of MD then showed the effect is entirely contained in that static
> binary, with dynamics adding nothing.

---

# Independent cross-check

## E25 — Cross-check against an independent pipeline

**Goal.** A teammate ran a separate pipeline (HeidAI-immrep26) on the same
competition. Recorded for what it confirms, what it exposes as a gap, and where
the two disagree. **Nothing from it was adopted into the submission.**

**Independent replication, which is the most valuable part.** Two pipelines built
without contact reached the same numbers on the scaffolding: macro-AUC SE 0.0116
vs 0.0110; per-epitope SEs matching to the third decimal; baseline recomputation
0.8245 vs 0.8247. And **seven dead ends closed independently**: annotation
metadata, OLGA Pgen, background normalisation, graph topology, chain
decomposition, a learned representation failing to beat classical distances
(ESM-2 650M at 0.8182 for them, SCEPTR for us), and most importantly **inner
hyperparameter search damaging cross-epitope transfer** — their result went
negative with a 400-point grid and positive when only one parameter was tuned,
the same effect E15 measured at 0.0048.

**Where they are ahead.** Candidate-pool expansion: they ran what E12 rejected on
design grounds, quantified the private-realisable gain by downsampling (+0.011,
or +0.013 with IEDB) rather than reporting the pool-inflated public number, and
used **targeted** expansion — merging IEDB only into epitopes whose VDJdb
multiple is below 2 — lifting the LTDEMIAQY pool from 1.04× to 1.72× at exactly
zero public cost. **We never examined IEDB, McPAS or ITRAP at all. That is a real
gap.** E12's error was conflating "this will not help the private side much"
(true) with "this is a trap" (overstated).

**Where we are ahead.** They list a self-run structure branch as pending with
negative expected value, reasoning that the AF3 baseline's 0.787 trails
tcrdist3's 0.824. **That is a standalone comparison, not a fusion test.** E14 and
E21 ran all 2,636 complexes: +0.007 honest nested gain, and the gain **triples on
the epitopes where density collapses** while staying flat on the robust ones.
Their own diagnosis is that the private side is pool-limited — exactly the
condition under which density weakens and structure pays. Also unique to this
pipeline: the pairwise AUC error decomposition (E22) and the near-positive
discordant clonotypes (E20).

**A caution on comparing ledgers.** Their best public figure (0.8552) and ours
(0.8498) are not comparable, because theirs is lifted by pool expansion the
private epitopes mostly cannot receive. On a private-realisable basis theirs is
roughly 0.840 and ours is the 0.846 nested estimate, which needs no discount.
**Mixed denominators inside a single ledger are easy to misread.**

**One framing of theirs that is better than ours.** They map the paper's causal
taxonomy of non-validation onto similarity-space signatures: multimer background
binding gives an isolated outlier that density catches; α-chain mispairing gives
a typical β with a random α; single-base sequencing error or V-gene
misannotation puts the TCR *immediately adjacent to a true positive*; and genuine
low affinity below the Jurkat threshold does the same. The last two are the same
population as our near-positive discordant clonotypes, but their frame explains
*why* that population must exist and why no sequence method can reach it.

## E26 — Where does the whole-V-domain alignment gain come from?

**Goal.** E25 replicated a teammate's alignment distance over the assembled
variable domains and found it the most decorrelated ensemble member, leaving a
candidate submission pending on a +0.0026 blend gain. Attribute that gain to a
region of the domain, with the **pre-declared criterion** that if framework /
V-gene explains all of it, the result is germline confounding rather than
recognition signal.

**Method.** CDR3 is found verbatim inside `full_seq` for all 2,636 rows, at a
median 79% into the domain and occupying about 12% of it, so the split is clean:
full domain, CDR3 masked out (framework + CDR1/2, the other 88%), and CDR3 alone.

**Result.**

| region | mean AUC | within-epitope Spearman with full |
|---|---:|---:|
| full V domain | 0.8379 | — |
| CDR3 masked out | 0.8287 | **+0.990** |
| CDR3 only | 0.7999 | +0.730 |

**The full-domain ranking is the CDR3-masked ranking.** CDR3 contributes almost
nothing to the order; the signal sits in the germline-encoded 88%. The
decorrelation comes from the same place.

**Conclusion.** The pre-declared criterion is met, so the honest reading is
**germline confounding**, not a new comparison of CDR3 chemistry. There is a
second, less comfortable reading: V-gene usage tracks study and donor and
correlates with convergent-recombination propensity, so this may be re-encoding
"how typical is this TCR" for the **third** time, after density itself and after
structure confidence.

**Decision: `submissions/density4_struct03.csv` is not adopted.** Already below
the +0.02 threshold at +0.0026, it now also has a mechanism arguing against it.

**Worth noting what this says about the method, not just the feature:** a result
that looked like useful decorrelation survived a correlation check, a nested LOEO
check and a 5-of-5 fold-selection check, and **only failed when asked which part
of the molecule it came from. Statistical screening did not catch it; biological
attribution did.**

---

# Review response: E40 and E41

An external review of the 2026-10-23 edition made three criticisms that proved
correct, and running them down changed two of the project's headline claims.

## E40 — Does E30's private projection transfer? It does not

**Goal.** Check the reviewer's charge that E30's private composition is not
label-free. `scripts/37_private_projection_validation.py`.

**First finding, before any statistics: E30 was computed in a heredoc that was
never saved.** Its log entry cites no script, so the project's most consequential
number rested on a computation absent from version control.

**The composition is not label-free — and the charge holds.** A class such as
`pos_clu_neg_sing` needs to know which pair member is positive. With
*a* = P(clustered | positive), *b* = P(clustered | negative), π = P(positive):

```
pos_clu_neg_sing = a(1−b)       pos_sing_neg_clu = (1−a)b
both_singleton   = (1−a)(1−b)   clustered-both   = ab
```

and the only label-free observable is q = πa + (1−π)b. One equation, two
unknowns. **π is a third unknown the review did not name**, running 0.287 to
0.663 across the public epitopes; splitting the both-clustered mass needs a
fourth. λ = a/b spans 1.9 to 10.2 where it is defined at all.

**The reconstruction E30 cited as validation is an identity.** Recombining each
epitope's own error rates with its own shares returns its own AUC: residual
**1.1×10⁻¹⁶**, machine epsilon.

**A basis mismatch nobody had noticed.** `scripts/20` clusters only the
*labelled* rows; `scripts/25`, whose family rates E30 reasoned from, clusters
the *whole* cohort. E30 combined error rates from one basis with compositions
from the other.

**The test: leave-one-epitope-out transfer.** Estimate everything transferred on
four epitopes, predict the fifth from its label-free q alone.

| held out | actual | full transfer | oracle composition | oracle errors |
|---|---:|---:|---:|---:|
| GILGFVFTL | 0.869 | 0.954 | 0.876 | 0.938 |
| GLCTLVAML | 0.741 | **0.954** | 0.833 | 0.929 |
| NLVPMVATV | 0.865 | 0.856 | 0.902 | 0.782 |
| TTDPSFLGRY | 0.840 | 0.847 | 0.852 | — |
| YLQPRTFLL | 0.933 | **0.830** | 0.864 | 0.912 |
| | | **MAE 0.083** | MAE 0.043 | MAE 0.090 |

**The method cannot predict a held-out public epitope.** Its out-of-sample MAE
of 0.083 is nearly twice the 0.047 bias it claimed to detect. Splitting the two
transfers shows the composition is the weaker half — exactly the
identifiability problem, now measured.

**And done honestly, it does not give 0.80.** Sweeping the unobservables over
the public range gives a private macro of **0.829–0.857**, median 0.846,
indistinguishable from LOEO. The mean transfer bias is +0.038 with a standard
error of about 0.054, so bias-correcting is not justified.

**A factual error in E30, which the reviewer caught.** It said all five public
epitopes are family-richer than the private three. False: TTDPSFLGRY at 0.047
sits below LTDEMIAQY at 0.146. Two of three private epitopes are below every
public one; the third is inside the public range.

**Conclusion — the revised position.**

| number | what it is |
|---|---|
| **0.846** | the estimate: nested LOEO over five public epitopes |
| 0.799 | **withdrawn** — assumptions never recorded, and it fails its own out-of-sample test |
| 0.67–0.91 | retained as a **scenario range**, not a confidence interval |

E30's qualitative finding survives: the private cohorts differ from the public
ones on a label-free observable, and LOEO cannot see it. What does not survive
is any number attached to that difference. **With five epitopes the size of the
shift is not measurable**, which is itself the result and agrees with E21.

## E41 — Is E31 robust to non-independence? Partly

**Goal.** E31 treated each mutated site as a draw. Re-test with the clonotype
family as the sampling unit. `scripts/38_e31_family_robustness.py`.

**The effective sample size was overstated about tenfold.**

| group | sites | pairs | TCRs | **families** |
|---|---:|---:|---:|---:|
| control | 767 | 522 | 302 | **79** |
| discordant | 103 | 71 | 87 | **33** |

The three largest families contribute 129, 92 and 73 sites on their own, and 15
families appear in both groups.

**On E31's full raw-structure set**, with a cluster bootstrap over families and
a permutation that shuffles the group label between families:

| quantity | effect | family 95% CI | family p | E31's site-level p |
|---|---:|---|---:|---:|
| peptide contacts at site | +3.629 | [+0.43, +7.83] | **0.043** | 0.0070 |
| site faces the peptide | +0.109 | [−0.03, +0.27] | **0.236** | 0.0790 |

**On the 40-pair relaxed subset E37 and E39 used** (126 sites, 30 families):

| quantity | effect | family 95% CI |
|---|---:|---|
| peptide contacts at site | +1.815 | **[−1.16, +5.19]** crosses zero |
| site faces the peptide | +0.182 | [+0.04, +0.38] |

**The two samples disagree about which reading survives.** All four effects are
positive, but the contact count clears family-level resampling only on the full
set and the binary facing rate only on the relaxed subset. **No single reading
survives in both** — which undercuts the E37 revision's "the whole effect is the
binary one", a conclusion read off the relaxed subset alone.

**Leave-one-epitope-out does not reverse the direction**, the one thing that
holds: the contact effect runs 2.35 to 4.60 across the five drops against 3.63
on the full data. Direction is stable; magnitude is not pinned.

**Conclusion, revised wording for E31.** *Near-identical TCR pairs with
discordant functional labels show an enrichment of substitutions at predicted
peptide-facing positions relative to concordant pairs. The direction is
consistent across epitopes and across both structural samples. With the
clonotype family as the sampling unit the enrichment is significant in one
reading per sample and the two samples disagree on which, so the effect size is
not established and causality is not demonstrated.*

## E42 — Framework or CDR1/2? The last low-cost exploration

**Goal.** The reviewer objected that E26's germline-confounding reading is not
forced: CDR1 and CDR2 are germline-encoded **and** genuinely contact
peptide-MHC. Separate the two. `scripts/39_v_region_ablation.py`.

**Method.** Split the CDR3-masked region into the part that cannot touch the
interface (framework) and the part that can (CDR1 + CDR2 + CDR2.5). The loops
come from tcrdist3's IMGT-derived reference keyed on the V gene, so no numbering
scheme is invented: 0 genes absent, CDR1 located verbatim in 99.5% of chains,
CDR2 98.0%, CDR2.5 99.98%. E26's own split was another unsaved heredoc, so all
four regions are recomputed over a common σ grid.

| region | mean α length | σ=0.02 | 0.05 | 0.10 | 0.20 | best |
|---|---:|---:|---:|---:|---:|---:|
| full V domain | 111 | 0.8243 | 0.8304 | **0.8354** | 0.8343 | 0.8354 |
| CDR3 masked | 98 | 0.8199 | 0.8260 | **0.8276** | 0.8264 | 0.8276 |
| **framework only** | 79 | 0.8154 | 0.8246 | **0.8270** | 0.8223 | **0.8270** |
| loops only | 18 | 0.8074 | 0.8092 | 0.8122 | **0.8173** | 0.8173 |

Within-epitope Spearman against the CDR3-masked ranking: **framework +0.994**,
loops +0.970.

**Conclusion — the alternative is tested and rejected.** Framework alone
recovers 0.8270 of 0.8276 at rank agreement +0.994, from a region that makes
**no contact with peptide or MHC**. The loops, which do contact both, recover
less and agree less. E26's reading is now demonstrated rather than inferred:
the ranking reads V-gene identity through the scaffold, not recognition
chemistry. This does not say CDR1/CDR2 are biologically unimportant — it says
*this feature's ranking* is explained by a region that cannot participate in
recognition at all.

---

# Predictor frozen — 2026-10-25

No further feature search, configuration selection or submissions.

| | |
|---|---|
| submitted file | `submissions/density3_struct03.csv` |
| sha256 (first 16) | `523a4692df02467b` |
| built by | `scripts/18_final_submission.py` |
| fitted parameters | none |
| public score | 0.8498 (leaderboard 0.849) |
| private estimate | **0.846**, downside risk real but unquantified |
| deadline | 2026-11-02 15:00 UTC |

**Why freeze here.** Every remaining avenue is closed on evidence rather than
judgement: the signal saturates at 0.83–0.84 in every representation (E6, E16,
E17, E18); supervised routes lose to variance on five folds (E10, E13, E32);
pool expansion has no measurable gain on this configuration (E27, E33, E35,
E38); the structural line is closed with its mechanism identified (E36) and its
one mechanistic lead reduced to an unestablished effect size (E41). E15 priced
further configuration selection at −0.0048, so additional search has negative
expected value.

**What remains open is not a modelling question.** Whether substitutions at
peptide-facing positions *cause* the functional discordance cannot be settled by
more computation on this dataset — E41 showed the statistical support does not
survive consistently at the clonotype-family level, and E34 showed the
structures it rests on are validated only in a memorised regime. That question
is left to independent data or functional experiments.

---

## Other corrections from the same review

- **E20's "zero fitted parameters, so nothing to overfit" is too strong.** There
  are no fitted coefficients, but feature and configuration selection are still
  choices made against the public epitopes, and E15 priced that at 0.0048.
- **E21's "cost of 0.3 is bounded at 0.0002"** applies to the five measured
  epitopes, not in general.
- **E26's germline-confounding reading was one of two live explanations** —
  **settled by E42 above**, in E26's favour: framework alone, which contacts
  nothing, recovers the ranking at Spearman +0.994.
- **E39's "no further structural dynamics will add to it" overreaches.** What
  was shown is that 2 ns of implicit-solvent dynamics around the predicted pose
  adds nothing beyond the static facing rate.

---

# Self-review findings

Mistakes caught in my own work, listed because several of them had already
produced a number that looked like a result. This is the part of the record most
likely to be useful to someone repeating the work.

| where | what was wrong | effect |
|---|---|---|
| E6 | density computed on *labelled* rows only, where inference uses the whole cohort | reverses that table's operator ordering; conclusion and submissions unaffected |
| E16 | two wrong diagnoses before the right one — scale coupling, then loss shape | the real cause is that the surrogate only aligns with AUC as it stops being differentiable |
| E19 | claimed 5 diffusion samples were immune to the canonicality confound because they share sequence input | wrong: they share the computed trunk representation, so consistency is *downstream* of trunk confidence |
| E20 | leaned on a 3.9% agreement between two rates with different denominators | coincidence; carries no weight |
| E20 | inferred private AUC from cohort geometry | **disproved by E21**: sparsity does not determine separability |
| E20 | relabelling arithmetic presented without its selection effect | the rows were chosen *because* they rank high |
| E27 | private-gain point estimate computed on a model not submitted | +0.0013 → **+0.0004**, hiding a 26-fold spread and optimistic subsampling |
| E34 | pooled contact counts and minimum distances into one Spearman | inflated it; the blocks are inversely related and occupy different ranges |
| E35 | claimed in narration, before measuring, that the gain would be smaller on the organisers' snapshot | it is 0.0012 smaller — inside the noise |
| E36 | staged phase 2 behind a positive phase 1, then skipped it | the argument was sound but answered a different question; running it showed the failure is **universal** |
| E37 | three rows reported as independent confirmations | with 77% ties they are substantially one test |
| E37 | descriptor table computed in an unsaved heredoc | three of six numbers not reproducible from the CSV |
| E39 | mutated-site residues taken as the **union** of both pair members | 5 of 40 pairs mis-probed, all discordant; a smoke test showed a clean 6.2 vs 0.0 that was entirely contamination |
| E39 | replica-spread guard compared against the median per-site replica sd | zero for most sites, so every quantity "exceeded the spread" by construction |
| E39 | static-to-MD join not deduplicated on the structural key | 126 rows inflated to 150, manufacturing an apparent reversal and a p=0.047 |
| E39 | groups not disjoint at site level | 3 of 117 sites are discordant in one pair and control in another; applies to E37 equally |
| E30 | private projection quoted as a measurement when its composition is not label-free | **withdrawn by E40**: fails leave-one-epitope-out transfer at MAE 0.083 |
| E30 | "reproduces every public AUC to 0.0003" offered as validation | an identity; residual is machine epsilon |
| E30 | "all five public epitopes are family-richer than the private three" | false — TTDPSFLGRY 0.047 < LTDEMIAQY 0.146 |
| E30 | error rates and compositions taken from two different clustering bases | `scripts/20` clusters labelled rows, `scripts/25` the whole cohort |
| E31 | each mutated site treated as an independent draw | 767 control sites come from 79 families; family-level CIs are far wider |
| E23 | "fixing the adversarial class is worth at most +0.0081" | that is the best tested weighting; the class ceiling is +0.032 |

Three of these had already produced a publishable-looking number. The recurring
cause is computation in an unsaved heredoc rather than in version control.

---

# Final state

**Predictor frozen 2026-10-25. 46 of 49 plan cards executed.** The three not run are wet-lab validation
(excluded by the PI) and two cards requiring artifact exchange with the teammate
(declined).

| | |
|---|---|
| experiments recorded | **42** (E1–E42, plus two E37 revisions) |
| scripts | 39 |
| reports | 17 |
| experiment log | 1,920 lines |
| data-layer audit | 57 checks, **all passing** |
| GPU work | 2,636 single-sample folds (7.4 h), 2,636 five-sample folds (13.2 h), 2,706 decoy folds (7 h), 414 ns MD (14 h) |

**The central scientific finding.** The extractable signal is **local
neighbourhood density inside the epitope cohort**, and it saturates at 0.83–0.84
in every representation tried. Roughly 24 negative results, each with an
identified mechanism rather than merely a number.

**Two positive mechanistic results.** E31/E37, now restated as positional:
discordant pairs differ at a peptide-facing position about twice as often, which
gives sequence methods a biologically grounded ceiling. And E38/R04: pool
expansion carries real antigen-specific information, since all three
size-matched negative controls actively *hurt*.

**The headline revision.** E30 — LOEO holds out epitope identity but not cohort
structure — a real, label-free difference LOEO cannot see. But **E40 withdrew
the 0.80 attached to it**: the projection fails its own leave-one-epitope-out
test. **Private estimate 0.846**, downside risk real but unquantified, with
0.67–0.91 as a scenario range.

**One item outstanding.** Verify in the Kaggle UI that the final selected
submission is `density3_struct03.csv` and not the earlier density-only file.

---

## Data provenance and licensing

TCRvdb functional-validation labels are request-only from `tcrvdb.nkiai.nl`.
**Redistribution or public sharing of the dataset, in whole or in part, is not
permitted** — the file is gitignored and never committed. Work using it must cite
doi:10.1016/j.immuni.2026.09.007.

VDJdb releases 2024-06-13, 2025-02-21 and 2026-06-03; IMGT/HLA allele sequences;
RCSB PDB for the 14 solved complexes in E34.
