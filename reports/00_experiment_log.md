# Experiment log

Every experiment run on IMMREP26, in order, with what it was asking and what it
settled. Numbers are mean per-epitope ROC AUC over the five public epitopes
unless stated. Reference points: organisers' tcrdist3 baseline **0.8247**,
AlphaFold3 min-PAE baseline **0.787**. Standard error of a five-epitope
estimate is about **0.011**; on the three private epitopes it is about
**0.023**.

---

## E1 — Is the competition's metric knowable without the Evaluation page?

**Method.** Score the organisers' own `tcrdist_predictions.csv`, whose public
score is 0.824, against the published TCRvdb labels under every plausible
metric. `scripts/02_confirm_metric.py`.

**Result.** Exactly one candidate matches:

| candidate | value |
|---|---:|
| **mean per-epitope ROC AUC, ambiguous dropped** | **0.824728** |
| mean per-epitope average precision | 0.8338 |
| n-weighted mean per-epitope ROC AUC | 0.8466 |
| pooled ROC AUC | 0.8577 |
| pooled average precision | 0.8602 |
| mean per-epitope AUC, ambiguous as positive | 0.7908 |
| mean per-epitope AUC, ambiguous as negative | 0.7769 |

**Settled.** Only within-epitope ranking is read, so the 28.7%–66.3% spread in
validation rate across epitopes needs no calibration. Each epitope carries
equal weight regardless of size, so the 130-row LTDEMIAQY counts as much as the
243-row CINGVCWTV. Later confirmed end to end: a submission predicted locally
at 0.8396 scored 0.839. **The local harness is the leaderboard**, so no
submission slot ever needs to be spent measuring public-epitope performance —
which matters, as the daily limit is 1.

Also resolves the competition text's internal contradiction: the "~68%" public
split is 1,782/2,636 = 67.6%, i.e. the 2,034 public rows minus 252 ambiguous,
not the 2,034/2,636 = 77.2% the row counts suggest.

---

## E2 — Can the public labels simply be recovered? (asked by the PI)

**Method.** Join TCRvdb to `test.csv` and score the recovered label.

**Result.** `tcr_id` is shared **verbatim**: 2,034/2,034 rows match exactly, no
fuzzy matching needed. Submitting the recovered label scores **1.000 on every
public epitope**. Simulating partial match failure, even 5% unmatched rows
still gives 0.9988.

**Settled.** The public leaderboard measures a `merge`, not skill. The 1.000
and 0.997 standings are reproducible in one line. Public standing was dropped
from consideration entirely after this.

**Follow-up — do the private labels leak too?** No. TCRvdb's changelog lists
exactly which epitopes each release contains: YLQ+GLC (2025-05-01), plus TTD,
GIL and NLV (2026-10-01). CINGVCWTV, LLWNGPMAV and LTDEMIAQY have never been
released, are absent from the downloaded file, and no deposit containing them
was found. **The private leaderboard is sound**, so modelling is what decides
the ranking.

---

## E3 — Can annotation provenance rank the private epitopes at all?

**Method.** For every VDJdb metadata field, count distinct levels *within* each
epitope cohort, and flag a field as unusable if it is constant, or if its
dominant level holds >95% of rows, in any private epitope.
`scripts/01_audit_metadata_viability.py`.

**Result.** **Not one field survives.** `vdjdb_score` is constant 0 for
LLWNGPMAV and LTDEMIAQY; `method.verification` is null for all 602 private
rows; LTDEMIAQY is 97.8% a single `reference.id` and a single
`method.identification`.

**Settled.** A provenance feature can look strong in cross-validation on the
public epitopes — where GILGFVFTL and NLVPMVATV span `vdjdb_score` 0–3 — and
contribute literally zero on the private side. **Every new feature now gets
checked for within-epitope variance on all three private epitopes before its
LOEO score is believed.**

---

## E4 — Metadata-only baseline: how much does evidence quality explain?

**Method.** Score `vdjdb_score` directly, pooled and within epitope.

**Result.** Pooled AUC **0.455** — below chance — and the direction *flips
between epitopes*: 0.581 on GLCTLVAML, 0.399 on NLVPMVATV, 0.443 on
YLQPRTFLL, undefined on TTDPSFLGRY (constant). The 158 score-3 NLVPMVATV
entries validate at 15.8%, against 23.5% for score-0 in the same cohort.

**Settled.** VDJdb's own confidence score **anti-correlates** with functional
validation. A model that learns "higher confidence means more reliable" on one
epitope is actively harmed on the next. Metadata is kept as a documented
ablation, never as a submission component.

---

## E5 — What is the signal, structurally?

**Method.** Within each epitope, compare median distance from each TCR to its
nearest validated neighbour, split by label; cluster at tcrdist ≤ 90 and
measure label purity and the singleton rate.

**Result.**

| epitope | pos → nearest validated | neg → nearest validated | singletons | singleton validation rate vs cohort |
|---|---:|---:|---:|---|
| GILGFVFTL | 24 | 169 | 34% | 0.31 vs 0.66 |
| YLQPRTFLL | 42 | 198 | 58% | 0.10 vs 0.41 |
| NLVPMVATV | 24 | 212 | 70% | 0.13 vs 0.29 |
| TTDPSFLGRY | 108 | 183 | 83% | 0.39 vs 0.50 |
| GLCTLVAML | 24 | 144 | 37% | 0.48 vs 0.66 |

Clusters of ≥3 members have median label purity 0.89–1.00.

**Settled.** Validated TCRs sit in tight, label-pure clusters; non-validating
ones are isolated. Every method in this competition is really answering **"is
this TCR part of a convergent, antigen-driven cluster"**. GLCTLVAML is the
weak epitope everywhere because being a singleton barely predicts anything
there (0.48 vs a 0.66 cohort rate), and density scores its singletons at AUC
0.482 — no information at all.

---

## E6 — Does the representation matter?

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

**Settled.** A learned paired-chain representation beats a hand-built distance
by 0.007, inside one standard error. **The metric space does not matter.** Two
durable facts fell out: the **α chain carries more of this signal than β**
(0.799 vs 0.751), inverting the usual expectation; and α-only is the only
space that helps GLCTLVAML.

---

## E7 — Generation probability

**Method.** OLGA Pgen per chain and paired, as a neighbour-free alternative.
`scripts/04_pgen.py`.

**Result.** `-log10 Pgen(paired)` scores **0.417**, i.e. *higher* Pgen means
more likely to validate — the opposite of the hypothesis that spurious
annotations are typical repertoire draws. Not orthogonal either: Spearman
correlation with density reaches −0.63 on GILGFVFTL and −0.59 on TTDPSFLGRY.
Within singletons, AUC 0.50–0.53.

**Settled.** Dense convergent clusters are built from easy-to-generate
sequences, so Pgen mostly re-encodes density, with an unstable sign. Rejected.

---

## E8 — Background normalisation (TCRNET/ALICE-style enrichment)

**Method.** Divide cohort density by density against the other epitopes' TCRs,
to separate antigen-driven convergence from "this region of TCR space is just
crowded".

**Result.** Own-cohort density 0.8396; difference 0.8357; log ratio 0.8262.
Background density **alone** scores 0.60 *positive*.

**Settled.** Crowded regions of TCR space are where genuinely public clonotypes
live, so normalising subtracts signal rather than noise. Rejected.

---

## E9 — Cluster shape: can a real motif family be told from a co-annotated artifact?

**Method.** Prompted by a diagnostic: NLVPMVATV's 23 high-density
non-validating TCRs cluster tightly with *each other* (median cosine 0.148) but
sit far from that epitope's validated TCRs (0.312), whereas on other epitopes
the high-density mistakes sit right on the positives (GILGFVFTL 0.050,
GLCTLVAML 0.088). So: measure cluster size, internal tightness, V/J purity,
per-position CDR3 entropy, centrality.

**Result.** Restricted to rows inside a cluster, per-position CDR3β entropy
scores **0.773** against density's 0.768, and internal tightness scores
**0.275** — meaning **tighter clusters validate less**. Imputed across the full
cohort, though, entropy collapses to 0.5553.

**Settled.** A real signal with a clean mechanism — convergent recombination
produces loose, diverse motif families, while a knot of near-duplicates that
entered VDJdb together is tight, and plain density scores both as "dense" — but
it is largely a restatement of cluster membership and adds nothing on top.

---

## E10 — Fitted fusion

**Method.** Logistic regression and gradient boosting over all density and
cluster-shape features, LOEO. `scripts/07_loeo_fusion.py`.

**Result.** Best single unsupervised feature 0.8362; logistic regression on
density only 0.8359; on density + cluster shape 0.8299; depth-2 booster 0.8188;
cluster shape alone 0.6846.

**Settled.** Five folds cannot support a fitted combiner. **Hard constraint
adopted: any fusion is a zero-parameter rank average**, and a component enters
only if it adds >0.02 in LOEO (`scripts/12_final_fusion.py` enforces this).

---

## E11 — Graph refinement

**Method.** Replace degree centrality (plain density) with eigenvector
centrality, PageRank, and 1–10 rounds of label propagation on the similarity
graph.

**Result.** Degree 0.8362; PageRank(0.85) 0.8329; propagation ×1 0.8297, ×3
0.8095, ×10 0.7417; eigenvector centrality 0.7486.

**Settled.** Monotonically worse with more propagation. The signal is
**strictly local**; eigenvector centrality collapses onto the largest cluster
and stops ranking the others. Rejected.

---

## E12 — Can the density estimate be sharpened with unscreened VDJdb TCRs?

**Method.** Count same-epitope paired TCRs in VDJdb that are not in the test
cohort and could pad the reference set.

**Result.** Public: 3,005 extra for GILGFVFTL, 9,881 for NLVPMVATV (28×),
459 for GLCTLVAML. Private: **11 for CINGVCWTV, 5 for LTDEMIAQY**, 191 for
LLWNGPMAV.

**Settled.** Enrichment would inflate public scores while leaving two of the
three deciding epitopes untouched, destroying the correspondence between LOEO
and private performance — the same trap as E3. The organisers appear to have
chosen epitopes where VDJdb coverage is spent. Rejected on design grounds,
without being run.

---

## E13 — Epitope-conditioned specificity (requested by the PI)

**Method.** LightGBM over SCEPTR TCR embeddings concatenated with a BLOSUM62
peptide encoding, trained on 21,112 VDJdb paired pairs across 119 epitopes
(capped at 1,500 per epitope), negatives being the same TCRs paired with wrong
epitopes so TCR identity carries no label. **The target epitope is removed from
training entirely** — it must be, since the test TCRs are themselves VDJdb
entries annotated to it. `scripts/11_specificity_model.py`.

**Result.** NLVPMVATV 0.553, GILGFVFTL 0.533, GLCTLVAML 0.531, TTDPSFLGRY
0.507, YLQPRTFLL **0.283**. Mean **0.481 — below chance**. Rank-averaged with
density it drags 0.8415 down to 0.7293. Cost: 78 min wall, 1,153 CPU-min.

**Settled.** Predicting binding to a peptide never seen in training remains
unsolved. What the model does learn is *how generic a TCR looks against VDJdb*,
and on YLQPRTFLL the generic-looking TCRs are precisely the spurious
annotations — the same failure shape as `vdjdb_score` in E4. Two independent
failures with one mechanism: **a quantity whose sign is not stable across
epitopes is worse than no feature.**

---

## E14 — Structure

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

`iptm_mhc_pep` is the diagnostic: that interface contains no TCR and is
identical for every complex in a cohort, yet scores 0.647. Within-epitope
Spearman correlation with density: `confidence_score` +0.541,
`pae_tcr_pep_min` +0.528.

**Settled.** Structure confidence here largely reports **how canonical the TCR
is** — typical sequences fold confidently everywhere through the shared pair
representation — so it is not the orthogonal signal it was hoped to be. It
still adds a little: blended at weight 0.2–0.3, nested LOEO gives **+0.0070**
(0.8392 → 0.8461), with the weight stably selected at 0.2–0.3 in every fold.

Four of five epitopes improve; the sole regression is TTDPSFLGRY (−0.011),
**the only 10-mer in the competition**. Peptide length changes MHC-I groove
geometry, and all three private epitopes are 9-mers, so restricted to the four
9-mers the gain is +0.0115. Caution the other way: LLWNGPMAV's structure
confidence has sd 0.40 against 1.0–2.0 elsewhere, so structure has little to
rank with there.

---

## E15 — How much does *choosing a configuration* cost?

**Method.** 323 candidate configurations (representation, bandwidth, structure
feature, blend weight). Compare three quantities: the oracle (best config
scored on all five epitopes), nested LOEO (config re-chosen inside each fold),
and a config fixed a priori then scored by plain LOEO.
`scripts/13_strict_loeo.py`.

**Result.**

| | mean |
|---|---:|
| oracle — tuned on the evaluation set | 0.8508 |
| **nested LOEO — tune by CV, then submit** | **0.8460** |
| best a priori config, plain LOEO | 0.8508 |

**Settled.** Cross-validated selection is **0.0048 worse than fixing the
configuration in advance**. Mechanism, visible in the TTDPSFLGRY fold: on the
other four epitopes `pae_a_mhc_min@0.3` looked best and returned 0.8328 on the
held-out epitope, where the fixed config gave 0.8504. With five folds the
selection step has worse signal-to-noise than the thing being selected.

**Make fewer decisions, not better-cross-validated ones.** And: every earlier
number quoted for a *chosen* configuration — including 0.8415 — is optimistic
by roughly this much. The defensible private estimate is **0.846**.

---

## E16 — Learning the metric that density lives in

**Method.** Learn the kernel bandwidth alone (1 parameter) and a diagonal metric
over the 64 SCEPTR dimensions (64 parameters, shrunk toward identity) by gradient
descent on a within-epitope pairwise ranking surrogate. Strictly nested: fitted
on four epitopes with an inner LOEO for shrinkage, applied to the fifth.
`scripts/14_metric_learning.py`.

**Result.** Plain cosine sigma=0.1 **0.8396**; scale-only **0.8199**; diagonal
**0.8144**. Learned bandwidths came out at 0.20-0.26 in every fold, against a
true optimum of 0.07.

**Diagnosis, in three steps, because the first two were wrong.**

*First guess — scale coupling.* The original loss was
`softplus(-margin x 50)`, whose fixed scaling interacts with sigma, since sigma
also sets the scale of the density values. Standardising density within each
epitope before the loss removes that coupling. It did not help: sigma still
converged to 0.20-0.26.

*Second guess — loss shape.* softplus penalises negative margins roughly
linearly and is unbounded, so the loss is dominated by its worst-ordered pairs,
whereas AUC counts every pair equally. A bounded sigmoid surrogate should be
closer. Measuring both against AUC over a sigma grid killed this too: softplus
minimises at sigma=0.20 and sigmoid at sigma=0.30, and **both are negatively
correlated with AUC** across the grid (-0.39 and -0.51). The sigmoid was worse,
not better.

*What is actually going on.* Sweeping the sigmoid temperature shows the
surrogate only aligns with AUC as it stops being smooth:

| temperature | loss argmin (sigma) | corr(-loss, AUC) |
|---|---|---:|
| 1.0 | 0.30 | -0.51 |
| 0.2 | 0.20 | -0.37 |
| 0.05 | 0.15 | -0.12 |
| 0.02 | 0.10 | **+0.03** |

As the temperature goes to zero the sigmoid becomes exactly the 1-AUC step
function, and only there does its optimum approach the real one at 0.07 -- by
which point gradients have vanished almost everywhere.

The reason this problem is unusually hostile to margin surrogates is that the
score distribution is **heavy-tailed**: a minority of TCRs sit in dense clusters
and take extreme density values. AUC gives full credit for ranking those at the
very top; any smooth surrogate saturates there and is instead dominated by the
ambiguous bulk. So the surrogate is least sensitive exactly where the metric
gets most of its signal.

**Settled.** Gradient descent on a margin-based ranking surrogate is
structurally mismatched to this metric and this score distribution -- not an
implementation flaw. Gradient-free direct search on AUC is viable for one
parameter (that is the sigma grid, which gains nothing: nested CV 0.8353 against
an oracle 0.8396) and would almost certainly overfit four epitopes at 64.

**What the sweep did settle.** Per-epitope optimal bandwidth varies five-fold:
GILGFVFTL 0.15, GLCTLVAML 0.10, NLVPMVATV 0.03, TTDPSFLGRY 0.05, YLQPRTFLL 0.15.

---

## E17 — Label-free adaptive bandwidth

**Method.** E16 showed real per-epitope bandwidth heterogeneity. Set σ from each
cohort's own distance geometry instead — no labels, no selection, so it
transfers to private epitopes unchanged. Ten rules tested.

**Result.**

| rule | mean |
|---|---:|
| fixed 0.07 | **0.8397** |
| fixed 0.10 | 0.8396 |
| median pairwise / 12 | 0.8385 |
| median nn1 | 0.8054 |
| q10 pairwise | 0.7927 |

**Settled.** The rules that *work* barely adapt — `median pairwise/8` picks
0.091–0.115 across all epitopes, i.e. it is fixed σ in disguise. The rules that
genuinely vary adapt in the **wrong direction**: NLVPMVATV's optimum is 0.03
and `median nn1` hands it 0.48. Mechanism: singleton-dominated cohorts (NLV
70%, TTD 83%) have large nearest-neighbour distances and need *small* σ so the
few genuinely close pairs dominate, while any bandwidth-scaling heuristic gives
sparse cohorts *large* σ and washes the signal out. Bandwidth heterogeneity is
real but **not predictable from cohort geometry**. Fixed σ≈0.07–0.10 stands.

---

## E18 — Master cross-epitope comparison (PI's experiment 03)

Everything in one protocol: five public epitopes, LOEO, mean per-epitope AUC,
nothing tuned on the epitope it is scored on. `scripts/15_master_comparison.py`,
table in `04_master_comparison.md`.

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

---

## Where this leaves things

Eighteen experiments, one consistent picture: **the extractable signal is
local neighbourhood density inside the epitope cohort, and almost nothing
survives being layered on top of it.** Structure adds a real but small amount
because it partially re-encodes the same "is this TCR canonical" quantity.
Everything that reads annotation provenance, generation probability, or
peptide identity is at or below chance, and in two cases its sign flips
between epitopes, which is worse than having no feature.

Best honest estimate of private performance: **0.846** (nested LOEO), against
the organisers' 0.8247 baseline. The gap to the leaderboard's real competitors
is within the ±0.023 noise of a three-epitope private score, so final standing
will be substantially luck.

**Methodology that earned its keep, in order of how much it changed decisions:**

1. Check every feature for within-epitope variance on the *private* epitopes
   before trusting any cross-validation estimate (E3, E12).
2. Treat configuration selection as part of the model and measure its cost
   (E15). With five folds, fix choices a priori.
3. Distrust any feature whose sign flips between epitopes, however good its
   pooled score (E4, E13).
4. Diagnose before concluding — E16 looked like "learning fails" and was
   actually a broken surrogate.

5. Verify that a surrogate objective actually tracks the metric before
   optimising it (E16). Two plausible fixes were both wrong, and the real
   answer was that no smooth surrogate works here.
