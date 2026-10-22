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

## E19 — Docking consistency across diffusion samples

**Method.** Re-fold all 2,636 complexes with 5 diffusion samples (Boltz runs its
trunk once and samples only the diffusion head, so 5 samples cost 1.36x one:
18 s each, 13.2 h total). Superpose each pair of samples on the MHC and peptide
C-alphas, then take TCR C-alpha RMSD. The question is whether a pair docks to
the *same pose* when resampled, i.e. whether one binding mode exists.
`scripts/16`, `scripts/17`.

**The hypothesis, and why it was wrong.** E14 showed Boltz *confidence* largely
reports how canonical a TCR sequence is. Consistency was expected to be immune
to that, on the argument that all five samples share identical sequence input,
so any difference between them must come from the docking solution itself.

**That argument was wrong.** The five samples share more than the input: they
share the single already-computed trunk representation. The diffusion head
samples a distribution conditioned on it, so a confident trunk yields a narrow
conditional and therefore consistent samples. Consistency is *downstream* of
trunk confidence. The shared input is the reason it is not orthogonal, not a
reason it would be.

**Result.** Docking RMSD spans 0.32-42.7 A, so the dynamic range is real, but:

| feature | mean | GLCTLVAML |
|---|---:|---:|
| density (reference) | 0.8396 | 0.713 |
| `ms_iptm_mean` (5-sample mean confidence) | 0.7665 | 0.735 |
| single-sample min-PAE | 0.7565 | 0.737 |
| `consist_tcr_rmsd_mean` | 0.7405 | **0.775** |

Within-epitope Spearman against density is **+0.530** — indistinguishable from
the +0.528 of the confidence features it was supposed to improve on.

Blended into the current best (3-space density + structure @0.2 = 0.8508), the
best gain is **+0.0010**, against a pre-registered inclusion threshold of
+0.02. **Rejected.**

One thing worth keeping: consistency scores 0.775 on GLCTLVAML, the highest of
any single feature on the project's weakest epitope. But it buys nothing free --
weighting it heavily enough to lift GLCTLVAML to 0.765 drops the overall mean to
0.8128, because the gain is taken from the other four epitopes.

**Settled.** 13.2 GPU-hours for +0.001. The pre-registered rule did its job:
a mechanistic hypothesis was wrong, the data said so, and the threshold stopped
a near-useless component from entering the submission.

---

## E20 — Four questions from the PI about what density is actually doing

### Does the method lean on density more than on the labels?

It uses **no labels at all.** Density and structure are both unsupervised, and
all 1,782 labels serve only to validate. That is a strength -- there is nothing
to overfit, which is why the local estimate reproduced the leaderboard exactly
-- and an unresolved gap: 1,782 labels are available and the predictor extracts
nothing from them. Every supervised route tried (E10 fitted fusion, E13
specificity, E16 metric learning) lost to variance or scored below chance.

### Are clustered TCRs merely label-pure, or genuinely positive-leaning?

Genuinely positive-leaning. Of 45 clusters with >=3 members (SCEPTR cosine
<= 0.12):

| cluster validation rate | clusters | TCRs |
|---|---:|---:|
| 100% (pure positive) | 26 | 205 |
| 80-100% | 8 | 136 |
| 60-80% | 8 | 31 |
| 20-60% | **0** | **0** |
| 0-20% | 0 | 0 |
| 0% (pure negative) | 3 | 13 |

The 20-60% band is **empty**, and only 13 of 385 clustered TCRs (3.4%) sit in a
pure-negative cluster. So cluster membership is positive evidence in its own
right, not a statement about homogeneity. Positives are 41-62% clustered
against negatives at 0-30%.

The ceiling is equally clear: only 487 of 1,782 TCRs are in any cluster, and
**about 54% of positives are singletons**.

### Is anything available for singletons?

Density remains the best thing there -- mean AUC 0.758 within singletons,
against structure 0.722, multi-sample confidence 0.729, consistency 0.692, Pgen
0.497. "Singleton" at threshold 0.12 is not a binary state; the continuous
nearest-neighbour distance still varies and still ranks (`nn1_dist` alone
scores 0.757).

One exception, and it is the weak epitope. On GLCTLVAML's 81 singletons density
collapses to **0.551** while every structure-derived feature beats it:
multi-sample confidence 0.676, docking consistency 0.668, min-PAE 0.659.

**Near-neighbour label discordance (not, on reflection, a proven ceiling).** Some negatives are single-residue
variants of validated positives -- identical CDR3 on one chain with one or two
substitutions on the other, sitting inside otherwise-pure positive clusters:

| epitope | negatives | near-identical to a positive | % of that epitope's positives | AUC now | AUC if relabelled |
|---|---:|---:|---:|---:|---:|
| GILGFVFTL | 194 | 14 | 3.7% | 0.8708 | 0.9081 |
| GLCTLVAML | 56 | 8 | 7.3% | 0.7018 | 0.7495 |
| NLVPMVATV | 226 | 11 | 12.1% | 0.8384 | 0.8704 |
| TTDPSFLGRY | 172 | 0 | 0.0% | 0.8420 | 0.8420 |
| YLQPRTFLL | 224 | 1 | 0.6% | 0.9279 | 0.9310 |
| total | 872 | 34 (**3.9%**) | | 0.8362 | **0.8602** |

3.9% of negatives. The pooled screen's own reported false-negative rate is also
3.9%, but that is expressed against positives, not negatives, so the agreement
is a coincidence of denominators and carries no weight -- it should not have
been leaned on.

Two further cautions, both raised in review and both correct. The relabelling
arithmetic is **selection-affected**: these rows were picked *because* they sit
close to positives, which is the same property that makes the model rank them
high, so of course relabelling them raises AUC. And the discordance has two
possible causes that this analysis cannot separate -- assay dropout, or genuine
fine specificity, since one CDR3 residue can abolish recognition. So 0.024 is
not a demonstrated unreachable ceiling, and E22 shows the figure also does not
arise where it appeared to: same-cluster pairs contribute only 0.5% of the total
loss, and the damage from these rows is instead concentrated in
`pos_sing_neg_clu`.

The honest description is **near-positive discordant clonotypes**: negatives
indistinguishable from validated positives by sequence, of unresolved cause.

### Why is GLCTLVAML the weak epitope, and does that threaten the private side?

GLCTLVAML is the most clonally converged cohort in the competition. Density
fails there because **everything is dense**:

| epitope | split | median nn1 | fraction with nn1 <= 0.12 |
|---|---|---:|---:|
| GLCTLVAML | public | 0.096 | **0.552** |
| GILGFVFTL | public | 0.131 | 0.476 |
| YLQPRTFLL | public | 0.398 | 0.184 |
| NLVPMVATV | public | 0.481 | 0.183 |
| TTDPSFLGRY | public | 0.497 | 0.047 |
| LTDEMIAQY | PRIVATE | 0.553 | 0.146 |
| LLWNGPMAV | PRIVATE | 0.361 | 0.044 |
| CINGVCWTV | PRIVATE | 0.404 | 0.025 |

Half of GLCTLVAML's cohort has a neighbour inside the clustering radius, so the
feature has almost no dynamic range left; its mixed clusters turn out to be 26:1
and 15:2 positive, with the negatives being the single-residue variants above.

**None of the three private epitopes share this pathology.** A first pass using
top-V-gene share flagged LLWNGPMAV (0.61, higher than GLCTLVAML's 0.50) and
CINGVCWTV (0.44), but that is a red herring: those cohorts share a V gene while
their CDR3s stay diverse, and only 2.5-4.4% of their TCRs have a close
neighbour. The right diagnostic is the close-neighbour fraction, on which
GLCTLVAML and GILGFVFTL are the converged pair and all three private epitopes
look instead like TTDPSFLGRY, NLVPMVATV and YLQPRTFLL -- where density scores
0.84-0.93.

**Two consequences.** Do not build a GLCTLVAML-specific adaptive weighting: the
pathology does not transfer, so tuning for it is fitting on a single epitope.
And the 0.846 private estimate is probably slightly **pessimistic**, because
GLCTLVAML's 0.73 drags the five-epitope mean and nothing on the private side
resembles it.

---

## E21 — Cohort thinning: does density survive a sparser repertoire? (PI's proposal)

**Why it mattered.** Density is computed within the epitope cohort, so its
quality depends on how well that cohort covers the relevant TCR space. The
private cohorts are sparser than anything validated on: the fraction of TCRs
with a neighbour inside the clustering radius is 2.5% for CINGVCWTV and 4.4%
for LLWNGPMAV, against 4.7-55% across the public five. CINGVCWTV is below every
public epitope, so the LOEO estimate is an extrapolation.

**Method.** Randomly thin each public cohort to fractions 1.0 down to 0.12,
recompute density within the thinned cohort, score against the real labels,
40 repeats. Two variants, because the naive design confounds a worse density
estimate with a noisier AUC on fewer rows: `thin_both` drops TCRs and scores
the survivors; `thin_reference` drops from the reference set only and still
scores every labelled row. `scripts/19_cohort_thinning.py`.

**Result 1 — density's robustness to sparsity is wildly epitope-dependent.**
Within-epitope, AUC against the close-neighbour fraction thinning produced:

| epitope | dense end | sparse end | change |
|---|---:|---:|---:|
| GILGFVFTL | 0.868 (cf 0.449) | 0.869 (cf 0.130) | **0.000** |
| YLQPRTFLL | 0.928 (cf 0.194) | 0.920 (cf 0.013) | -0.008 |
| GLCTLVAML | 0.713 (cf 0.542) | 0.688 (cf 0.081) | -0.025 |
| TTDPSFLGRY | 0.845 (cf 0.041) | 0.724 (cf 0.002) | **-0.121** |
| NLVPMVATV | 0.844 (cf 0.205) | 0.699 (cf 0.014) | **-0.145** |

**Result 2 — and this corrects E20.** At essentially identical sparsity, cf
about 0.013-0.014, YLQPRTFLL scores **0.920** and NLVPMVATV **0.699**. Sparsity
does not determine performance; intrinsic separability does. So the claim in
E20 that the private epitopes "look like TTDPSFLGRY, NLVPMVATV and YLQPRTFLL
where density scores 0.84-0.93" **does not hold**: at private-like sparsity
those same three epitopes span 0.70-0.93. Cohort geometry cannot forecast a
private epitope's AUC, and the real uncertainty on the private score is wider
than the +/-0.023 sampling error alone implies.

**Result 3 — thinning is not strongly pessimistic.** The worry was that thinning
manufactures false singletons by deleting a positive's true cluster-mates. The
singleton-to-negative association barely moves: the gap between clustered and
singleton validation rates goes 0.486 at full coverage to 0.454 at 20%. So the
simulation is reasonably faithful rather than a loose bound.

**Result 4 — structure is insurance, and it changed the submission.** The
structure blend's gain concentrates exactly where density fails:

| | full | 45% | 20% |
|---|---:|---:|---:|
| gain on the two fragile epitopes (NLVPMVATV, TTDPSFLGRY) | +0.0059 | +0.0120 | **+0.0202** |
| gain on the three robust epitopes | +0.0112 | +0.0103 | +0.0099 |

Sweeping the weight at each coverage level:

| coverage | w=0.2 | w=0.3 | w=0.5 | best |
|---|---:|---:|---:|---|
| full | **0.8487** | 0.8485 | 0.8360 | 0.2 |
| 45% | 0.8263 | **0.8305** | 0.8252 | 0.3 |
| 20% | 0.8117 | 0.8125 | **0.8144** | 0.5 |

**Decision: structure weight 0.2 -> 0.3**, on a dominance argument rather than
by optimising a regime. w=0.3 ties w=0.2 at full coverage (-0.0002) and beats it
at both sparse levels (+0.0042, +0.0008), and the private cohorts are the sparse
case (cf 0.025-0.146, between the 20% and 45% thinning levels). w=0.5 wins at
20% coverage but loses 0.013 at full coverage, so its risk is asymmetric. If the
sparsity argument is simply wrong, the cost of 0.3 is bounded at 0.0002.

Public-epitope score moves 0.8508 -> 0.8498 under the new weight, which is the
expected small price for the insurance.

---

## E22 — Pairwise AUC error decomposition (reviewer's P0)

**Why.** Mean per-epitope AUC is exactly the probability a validated TCR
outranks a non-validating one from the same cohort, so 1 - AUC is a set of
specific inverted pairs. The inference this guards against: 54% of positives
are singletons, which does **not** imply 54% of the error involves them.

**Method.** Partition every positive-negative pair within an epitope by cluster
state, mutually exclusive and complete, so the class losses sum to 1 - AUC.
The script asserts that reconstruction; residual 1.4e-16.
`scripts/20_auc_error_decomposition.py`.

**Result.**

| class | share of loss | share of pairs | error rate within class |
|---|---:|---:|---:|
| both_singleton | **62.0%** | 50.5% | 0.207 |
| pos_sing_neg_clu | **21.2%** | **4.5%** | **0.673** |
| diff_cluster | 10.0% | 6.2% | 0.231 |
| pos_clu_neg_sing | 6.4% | 38.6% | 0.021 |
| same_cluster | 0.5% | 0.3% | 0.309 |

**The adversarial class is systematically inverted.** When the positive is
isolated and the negative is clustered, the submitted model gets the order
wrong 67% of the time -- worse than chance -- on GLCTLVAML 0.819 and NLVPMVATV
0.825. It is 4.5% of pairs and 21.2% of the error.

**This corrects E20's account of the 34 near-positive discordant negatives.**
Same-cluster pairs -- the fine-specificity comparisons -- contribute only 0.5%
of loss. So relabelling those negatives did not gain 0.024 by fixing
within-cluster discrimination. It gained it because a relabelled negative stops
being a negative in *every* comparison, and most of its damage was in
`pos_sing_neg_clu`: a clustered negative collects a high density score and
outranks isolated positives elsewhere in the cohort. The mechanism is
over-scoring of clustered negatives, not failure to resolve a neighbour.

---

## E23 — Can the adversarial class be fixed? No, and the arithmetic says why

**Per-class feature strength.** Structure is better than chance exactly where
density is catastrophically worse:

| feature | pos_sing_neg_clu | pos_clu_neg_sing |
|---|---:|---:|
| density | **0.697** | **0.013** |
| blend w=0.3 | 0.602 | 0.023 |
| structure min-PAE | **0.422** | 0.156 |
| multi-sample confidence | 0.447 | 0.125 |

Per epitope in the adversarial class, density against structure: YLQPRTFLL
0.638/0.219, GLCTLVAML 0.931/0.405, NLVPMVATV 0.968/0.468.

**Attempt 1 — cluster-conditional weighting.** Cluster state is label-free, so
weighting structure differently for clustered and isolated TCRs is legitimate at
inference. Strict nested LOEO over both weights: **0.8418 against 0.8486 for the
flat weight, -0.0068.** And the oracle over the whole 8x8 grid peaks at 0.8486
at (0.2, 0.3) -- *identical to flat w=0.3*. Even with perfect knowledge,
conditioning buys nothing.

**Why, quantitatively.** Scores are per TCR, not per pair, so a TCR cannot know
which comparison it is in. Raising structure weight to help the adversarial
class necessarily moves the same TCRs in the easy class:

| structure weight 0.3 -> 0.6 | contribution |
|---|---:|
| pos_sing_neg_clu (4.5% of pairs) | **+0.0040** |
| pos_clu_neg_sing (38.6%) | -0.0146 |
| both_singleton (50.5%) | -0.0161 |
| diff_cluster (6.2%) | -0.0029 |
| **net** | **-0.0296** |

Fixing the adversarial class completely, at w=1.0, is worth at most +0.0081
while costing -0.107 elsewhere. The trade is about 13:1 against.

**Settled.** The decomposition succeeded diagnostically and failed
prescriptively. What would be needed is a feature informative *within* the
clustered subset and roughly uninformative *between* subsets, so it can demote
clustered negatives without disturbing the clustered-versus-singleton ordering
that supplies 38.6% of correct pairs at a 2% error rate. Structure is not that:
it is informative everywhere, merely worse than density in most places.

---

## E24 — Is probabilistic label modelling available? No

The reviewer's GuanRank analogy suggests treating functional labels as evidence
of varying strength, modelling `P(Y_true = 1 | Y_observed, assay evidence)`
rather than a hard binary, and correctly flags that this depends on TCRvdb
exposing per-TCR detail.

**Checked.** The released TCRvdb file has 24 columns and the only label
information is `pos` / `neg` / `ambiguous`, strictly one-hot integers. There is
no enrichment value, read or UMI count, replicate, p-value, CD69 readout or any
other continuous quantity; the sole numeric field is `vdjdb_score`, taking
values 0-3, which is the *prior annotation* confidence and is already known to
anti-correlate with validation (E4). So the route is closed with the data
available, not merely unpromising. It would need a richer release from the
authors.

---

## E25 — Cross-check against an independent pipeline (teammate's HeidAI-immrep26)

A teammate ran a separate pipeline on the same competition. Recorded here for
what it confirms, what it exposes as a gap, and where the two disagree. Nothing
from it was adopted into the submission.

### Independent replication, which is the most valuable part

Two pipelines built without contact reached the same numbers on the scaffolding:

| | theirs | ours |
|---|---|---|
| macro-AUC standard error | 0.0116 | 0.0110 |
| per-epitope SE, GIL/YLQ/TTD/NLV/GLC | .016/.015/.022/.027/.041 | .014/.015/.021/.026/.039 |
| baseline recomputation | 0.8245 | 0.8247 |

And seven dead ends closed independently: annotation metadata (their S2, our
E3-E4), OLGA Pgen (M40, E7), background normalisation (S1b, E8), graph topology
(M70, E11), chain decomposition (M20), a learned representation failing to beat
classical distances (their M60 with ESM-2 650M at 0.8182, ours with SCEPTR), and
most importantly **inner hyperparameter search damaging cross-epitope transfer**
-- their M10 went negative with a 400-point grid and positive when only one
parameter was tuned, which is the same effect E15 measured at 0.0048.

### Where they are ahead, and where E12's reasoning was wrong

**Candidate-pool expansion.** E12 rejected padding cohorts with unscreened
same-epitope TCRs on design grounds, *without running it*. They ran it and did
two things we did not:

- Quantified the private-realizable gain by downsampling: **+0.011, or +0.013
  with IEDB**, rather than reporting the pool-inflated public number.
- **Targeted** expansion -- merging IEDB only into epitopes whose VDJdb multiple
  is below 2 -- which lifts the LTDEMIAQY pool from 1.04x to 1.72x at exactly
  zero public cost.

We never examined IEDB, McPAS or ITRAP at all. That is a real gap.

E12's error was conflating "this will not help the private side much" (true,
about +0.011) with "this is a trap" (overstated). It is only a trap for someone
reading the public score, which they explicitly refuse to do.

Their compact statistic for the same phenomenon: corr(cohort size, baseline
AUC) = +0.706 across the five public epitopes, with clustering measures nearly
uncorrelated. On five points that is weak evidence, but it points to the right
action, and it agrees with E21's within-epitope thinning curves.

### Where we are ahead

**Structure.** They list a self-run structure branch as pending with negative
expected value, reasoning that the AF3 baseline's 0.787 trails tcrdist3's 0.824.
That is a *standalone* comparison, not a fusion test. E14 and E21 ran all 2,636
complexes: +0.007 honest nested gain, and the gain **triples on the epitopes
where density collapses** (+0.006 to +0.020) while staying flat on the robust
ones. Their own diagnosis is that the private side is pool-limited, with
LTDEMIAQY smallest and sparsest -- and pool-limited is exactly the condition
under which density weakens and structure pays. So their diagnosis calls for the
insurance their prior rules out.

Also unique to this pipeline: the pairwise AUC error decomposition (E22) and the
near-positive discordant clonotypes (E20).

### Replicated and measured: whole-V-domain alignment distance (their M50)

Worth replicating because it carries no private discount -- it changes how two
TCRs are compared, not how many there are to compare against. Result in
`scripts/22_fullseq_alignment_distance.py`: standalone 0.8379, the same band as
SCEPTR's 0.8396 and tcrdist3's 0.8329, but the **most decorrelated** member
available (within-epitope Spearman +0.867 against SCEPTR and +0.862 against
tcrdist3, where those two correlate +0.909), and the best of the three on
GLCTLVAML at 0.734.

As a fourth density member the blend reads 0.8525 against the submitted 0.8498,
and strict nested LOEO selects it in all five folds. **Not adopted**: +0.0026 is
inside noise, so it is not a demonstrated improvement, and the candidate file
`submissions/density4_struct03.csv` is kept unsubmitted.

### A caution on comparing ledgers

Their best public figure (0.8552) and ours (0.8498) are not comparable, because
theirs is lifted by pool expansion that the private epitopes mostly cannot
receive. On a private-realizable basis theirs is roughly 0.8245 + 0.011 + 0.0044
= 0.840 and ours is the 0.846 nested estimate, which needs no discount. Mixed
denominators inside a single ledger are easy to misread.

### One framing of theirs that is better than ours

They map the paper's causal taxonomy of non-validation onto similarity-space
signatures: multimer background binding gives an isolated outlier that density
catches; alpha-chain mispairing gives a typical beta with a random alpha;
single-base sequencing error or V-gene misannotation puts the TCR *immediately
adjacent to a true positive*; and genuine low affinity below the Jurkat
threshold does the same. The last two are the same population as our
near-positive discordant clonotypes, but their frame explains *why* that
population must exist and why no sequence method can reach it, which is a better
account than our "false cluster" language. Worth carrying forward.

---

## E26 — Where does the whole-V-domain alignment gain come from? (plan card D02)

**Why.** E25 replicated a teammate's alignment distance over the assembled
variable domains, found it the most decorrelated ensemble member, and left a
candidate submission pending on a +0.0026 blend gain. A joint experiment plan
asked the right follow-up: attribute that gain to a region of the domain, with
the pre-declared criterion that if framework / V-gene explains all of it, the
result is germline confounding rather than recognition signal.

**Method.** CDR3 is found verbatim inside `full_seq` for all 2,636 rows, at a
median 79% into the domain and occupying about 12% of it, so the split is clean:
full domain, CDR3 masked out (framework + CDR1/2, the other 88%), and CDR3
alone. Same alignment, same normalisation, same bandwidth.

**Result.**

| region | mean AUC | within-epitope Spearman with full |
|---|---:|---:|
| full V domain | 0.8379 | -- |
| CDR3 masked out | 0.8287 | **+0.990** |
| CDR3 only | 0.7999 | +0.730 |

The full-domain ranking is the CDR3-masked ranking. CDR3 contributes almost
nothing to the order, and the signal sits in the germline-encoded 88%.

The decorrelation comes from the same place: `CDR3-masked` against tcrdist3 is
+0.812 while `CDR3-only` against tcrdist3 is +0.897. What made this distance the
most independent ensemble member is framework and CDR1/2 similarity -- which is
largely V-gene and allele identity, measured continuously.

**Settled.** The pre-declared criterion is met, so the honest reading is
germline confounding, not a new comparison of CDR3 chemistry. tcrdist3 already
carries V-gene identity as a discrete term, and the teammate's own ledger closed
a within-cohort V/J frequency feature at +0.0012 with the note that nndist
already covers it. There is a second, less comfortable reading too: V-gene usage
tracks study and donor, and correlates with convergent-recombination propensity,
so this may be re-encoding "how typical is this TCR" for the third time, after
density itself and after structure confidence (E14).

**Decision: `submissions/density4_struct03.csv` is not adopted.** It was already
below the +0.02 inclusion threshold at +0.0026; it now also has a mechanism that
argues against it. The submitted file stays `density3_struct03.csv`.

Worth noting what this says about the method, not just the feature: a result
that looked like useful decorrelation survived a correlation check, a nested
LOEO check and a 5-of-5 fold-selection check, and only failed when asked *which
part of the molecule it came from*. Statistical screening did not catch it;
biological attribution did.

---

## E27 — Pool expansion, audited and quantified (plan cards Q02, R01, R04-adjacent)

E12 rejected candidate-pool expansion without running it; E25 conceded that as a
gap. This runs it, and the gap turns out to be procedural rather than
substantive: the conclusion holds, but it is now measured instead of asserted,
and the mechanism is different from what E12 assumed.

### Q02 -- the query finds itself

Every test TCR came from VDJdb, so expanding with VDJdb risks counting the query
as its own neighbour. `scripts/23_pool_selfhit_audit.py`:

| epitope | queries finding an exact self-hit | mean self-hits | records per distinct clonotype | top-study share of pool |
|---|---:|---:|---:|---:|
| LTDEMIAQY | **1.000** | 1.000 | 1.000 | **0.978** |
| TTDPSFLGRY | 1.000 | 1.021 | 1.029 | 0.599 |
| YLQPRTFLL | 0.998 | 1.038 | 1.042 | 0.413 |
| GILGFVFTL | 0.990 | **2.193** | **1.650** | 0.349 |
| GLCTLVAML | 0.958 | 1.151 | 1.585 | 0.508 |
| NLVPMVATV | 0.945 | 1.383 | 1.249 | **0.948** |
| CINGVCWTV | 0.942 | 0.942 | 1.000 | 0.625 |
| LLWNGPMAV | 0.878 | 0.882 | 1.002 | 0.468 |

96.4% of queries find themselves, and the count is **not constant** -- GILGFVFTL
averages 2.19, so repeatedly curated clonotypes collect extra zero-distance
"neighbours".

**Curation multiplicity does track validation**, which is why that matters:

| VDJdb records for the clonotype | n | validation rate |
|---|---:|---:|
| 1 | 1470 | 49.0% |
| 2 | 195 | 43.6% |
| 3 | 39 | **82.1%** |
| 4+ | 77 | **93.5%** |

As a feature it is the familiar failure shape -- mean AUC 0.534 with the sign
flipping (GILGFVFTL 0.665, NLVPMVATV 0.433) -- and it is **constant on the
private epitopes**: `n_records` takes one distinct value on CINGVCWTV and
LTDEMIAQY. So it is both unusable and, left undeduplicated, a label proxy
masquerading as density.

Confirmed directly: raw expansion **loses** 0.0024 of mean AUC, while
query-removed expansion gains 0.0048.

### The aggregation operator decides whether expansion helps at all

Our kernel-mean density barely benefits; top-k does. Mean over the five public
epitopes, pool deduplicated:

| operator | no pool | pooled | delta |
|---|---:|---:|---:|
| top-8 | 0.8414 | **0.8526** | **+0.0112** |
| top-3 | 0.8334 | 0.8476 | +0.0142 |
| top-16 | 0.8377 | 0.8488 | +0.0111 |
| kernel sigma=0.1 | 0.8396 | 0.8432 | +0.0036 |

The mechanism is visible per epitope. NLVPMVATV, whose pool is 94.8% a single
dextramer study of ~10k records, goes **-0.0354 under kernel-mean and +0.0131
under top-8**: a mean over the whole pool is drowned by one study's records,
while top-k never looks at them. This explains the teammate's +0.0194 against
our +0.0048 -- they used top-k, which is the right operator once a pool is
involved, and our kernel is the right operator for a fixed cohort.

### What the private side can actually realise: +0.0013

Downsampling each public pool to the private epitopes' own *distinct-clonotype*
multipliers, 30 repeats, top-8:

| private epitope | multiplier | top-8 density | **submitted blend** |
|---|---:|---:|---:|
| LTDEMIAQY | 1.038 | +0.0002 | **-0.0003** |
| CINGVCWTV | 1.045 | +0.0002 | **-0.0002** |
| LLWNGPMAV | 1.834 | +0.0036 | **+0.0018** |
| **macro** | | **+0.0013** | **+0.0004** |

against +0.0112 on the full public pools. **88% of the public gain cannot
transfer**, because two of three private epitopes can add only about 4% more
distinct clonotypes.

**How that estimate was made, and three reasons not to quote the point value
alone.** For each private multiplier m, every public pool was topped up by
(m-1) x its own cohort size, 30 random subsamples, and the five resulting gains
averaged; the three private results were then averaged.

1. *It was first computed on a model we do not submit.* Top-8 is the operator
   where expansion works; our submission is kernel density blended 0.7/0.3 with
   structure, where the dilution is compounded. Recomputed on the submitted
   configuration the macro estimate falls to **+0.0004**, naive 95% interval
   [-0.0022, +0.0030] -- indistinguishable from zero, and negative on two of the
   three private epitopes.
2. *The mean hides a 26-fold spread.* At LLWNGPMAV's multiplier the five public
   estimates run GLCTLVAML +0.0158, GILGFVFTL +0.0038, YLQPRTFLL -0.0000,
   TTDPSFLGRY -0.0030, NLVPMVATV -0.0076 (sd 0.0089, five times the mean). And
   that interval is itself optimistic: averaging across public epitopes assumes
   the gain depends only on the multiplier, which is exactly the exchangeability
   E21 disproved. The defensible statement is that **five epitopes cannot
   estimate this quantity**, only that the point estimate sits near zero with
   the spread straddling both signs.
3. *Random subsampling is optimistic.* A random draw from GILGFVFTL's 3,005
   diverse clonotypes is not what LTDEMIAQY would get -- its 5 available extras
   come from a pool that is 97.8% one reference. Repeating the exercise drawing
   only from each pool's dominant study gives LTDEMIAQY **-0.0002**, also zero.

So the verdict stands and strengthens: not "the gain is small" but "on our
configuration there is no measurable gain".

This is nearly an order of magnitude below the teammate's +0.011 estimate. The
likely cause is the matching unit: downsampling by *records* rather than
*distinct clonotypes* lets a public subsample keep duplicate clonotypes, which
supply near-zero distances that the private pools -- at a record:clonotype ratio
of 1.000 -- have none of. Even allowing their targeted IEDB merge, which lifts
LTDEMIAQY to 1.72x, the estimate reaches only about +0.0025.

**Not adopted.** +0.0013 is a fifth of the structure gain and far below the
inclusion threshold. E12's verdict stands, though its reasoning was incomplete:
the problem is not only that few extra TCRs exist, but that the ones that do are
duplicates from single studies.

### A methodology slip this exposed

E6's representation comparison computed density using only the *labelled* rows
as the reference set, where inference uses the whole cohort. On the full cohort
top-8 reads 0.8414 against kernel's 0.8396, rather than the reverse ordering E6
reported. The conclusion is unaffected -- everything still lands in the same
0.83-0.85 band, so operator and space still do not matter much -- and the
submitted scripts always used the full cohort, so no submission was affected.
But the comparison table in E6 should be read with that caveat.

---

## E28 — Contact fingerprints and docking geometry (plan cards D07, C01, C02, C06)

All four extracted in one pass over the 2,636 structures,
`scripts/24_structure_contacts.py`. Peptide position is used as the common
coordinate -- the peptide is identical within an epitope -- so no IMGT numbering
is needed for the one thing that has to be exactly right.

**C01 geometry QC passes.** Peptide-MHC contacts: median 334, minimum 260, so
every peptide sits in the groove. CDR3 supplies 74-76% of all TCR-peptide
contacts, which is the biologically expected share and independently validates
the CDR3 localisation. 1,929 of 2,636 complexes carry at least one inter-chain
heavy-atom pair under 2.2 A, but the median count is 1, which is ordinary for
predicted structures. One complex, `1_6_B11_22_YLQPRTFLL`, has **no TCR-peptide
contact at all** -- a genuine docking failure, in a public epitope.

**D07 contact-fingerprint density: 0.7476**, below the confidence features and
far below sequence density, and correlated +0.563 with density. Another
re-encoding.

**But the raw contact counts are the most orthogonal quantity measured in this
project.** Within-epitope Spearman against density: `contacts_a_pep` **+0.088**,
`contacts_b_pep` +0.103. Everything else ever measured sat between +0.26 and
+0.91. That makes sense: how many peptide atoms a TCR actually touches is
constructionally unrelated to how crowded its neighbourhood is.

They are also too weak to use. Standalone AUC 0.588 and 0.559; best geometry
feature `tcr_height_over_peptide` (negated) 0.657; crossing angle 0.530. Under
strict nested LOEO with component and weight both chosen on the training
epitopes, the gain over the submitted blend is **+0.0008**. Rejected.

**General lesson worth keeping: orthogonality is necessary but not sufficient.**
A component at AUC 0.59 blended into a model at 0.85 contributes almost nothing
however independent it is -- the achievable contribution scales with its margin
over chance, and 0.09 does not buy anything. This closes the structural line
completely: interface confidence (+0.007, correlated), docking consistency
(+0.001, correlated), fingerprint density (worse, correlated), raw contact
counts (+0.0008, orthogonal but weak), docking geometry (weak).

---

## E29 — Family holdout and chain discordance (plan cards S06, D03)

**S06 -- leave-one-clonotype-family-out.** Each TCR scored against a cohort with
its own single-linkage family removed, so no relative can supply its density:

| epitope | in a family | standard | family-masked | delta |
|---|---:|---:|---:|---:|
| NLVPMVATV | 0.183 | 0.8436 | 0.7417 | **-0.1019** |
| GILGFVFTL | 0.476 | 0.8678 | 0.8520 | -0.0158 |
| GLCTLVAML | 0.552 | 0.7130 | 0.7000 | -0.0130 |
| TTDPSFLGRY | 0.047 | 0.8449 | 0.8356 | -0.0093 |
| YLQPRTFLL | 0.184 | 0.9285 | 0.9272 | -0.0012 |
| mean | 0.288 | 0.8396 | **0.8113** | **-0.0282** |

Family recurrence does carry part of the ranking, and once again the dependence
is epitope-specific rather than structural: NLVPMVATV and YLQPRTFLL have
near-identical family rates (0.183, 0.184) and lose 0.102 against 0.001.

**D03 -- the alpha/beta asymmetry is not mispairing.** Alpha-only density 0.799
against beta-only 0.751 (E6), and the mispairing hypothesis predicts that a
mispaired TCR shows large single-chain rank disagreement. Measured, the signed
disagreement scores 0.536 and the absolute disagreement **0.429 -- below
chance**, meaning *larger* alpha/beta disagreement associates with validating,
the opposite of what mispairing would produce. It has ample within-epitope
variance on all three private epitopes (sd 0.33-0.36, so it passes E3's rule)
and still adds only +0.0026 at weight 0.1, negative at 0.2. Rejected, and the
mispairing explanation for the chain asymmetry gets no support.

---

## E30 — Why leave-one-epitope-out overestimates the private score

S06 prompted a comparison of cohort structure rather than cohort size, and it
revealed something LOEO cannot see.

**The private epitopes have four times less family structure than the public
ones**, 7.2% against 28.8% of TCRs sitting in a family of two or more:

| epitope | split | in a family | implied cross-stratum pair share |
|---|---|---:|---:|
| CINGVCWTV | PRIVATE | 0.025 | ~0.048 |
| LLWNGPMAV | PRIVATE | 0.044 | ~0.084 |
| LTDEMIAQY | PRIVATE | 0.146 | ~0.250 |
| TTDPSFLGRY | public | 0.047 | ~0.089 |
| NLVPMVATV | public | 0.183 | ~0.299 |
| YLQPRTFLL | public | 0.184 | ~0.300 |
| GILGFVFTL | public | 0.476 | ~0.499 |
| GLCTLVAML | public | 0.552 | ~0.495 |

That matters because of E22: `pos_clu_neg_sing` -- clustered positive against
isolated negative -- is 38.6% of public pairs at a **2% error rate**. Those easy
pairs are what holds a public AUC up, and CINGVCWTV will have about 5.6% of them.

**Decomposing and reprojecting.** Taking each class's measured within-class
error rate for the submitted model and recombining with each epitope's own class
composition reproduces every public epitope's actual AUC to within 0.0003, so
the arithmetic is sound. Applying it to the private epitopes' label-free
compositions:

| epitope | both-singleton share | easy share | projected AUC |
|---|---:|---:|---:|
| CINGVCWTV | 92.6% | 5.6% | 0.795 |
| LLWNGPMAV | 87.2% | 9.8% | 0.797 |
| LTDEMIAQY | 60.7% | 29.2% | 0.806 |
| **macro** | | | **0.799** |

against the LOEO estimate of 0.846.

**So LOEO holds out epitope identity but not cohort structure.** It is unbiased
for "an epitope drawn like the public ones" and biased upward for an epitope
with private-like family sparsity. That bias is about **-0.047**, and no amount
of leave-one-epitope-out can detect it, because all five public epitopes are
family-richer than the private three.

**And the whole projection collapses to one unknown.** Since the private
epitopes are 61-93% both-singleton pairs, the projected macro tracks the
within-both-singleton AUC almost exactly. Substituting each public epitope's
value for that quantity:

| if private singletons separate like | both-singleton AUC | projected private macro |
|---|---:|---:|
| GLCTLVAML | 0.642 | **0.668** |
| NLVPMVATV | 0.748 | 0.758 |
| GILGFVFTL | 0.827 | 0.831 |
| TTDPSFLGRY | 0.830 | 0.832 |
| YLQPRTFLL | 0.915 | **0.907** |

**Revised position: central estimate 0.80, plausible range 0.67-0.91**, in place
of the single 0.846. The spread is not sampling noise -- it is one substantive
unknown, how separable the private epitopes' isolated TCRs are, and E21 already
showed that cannot be forecast from cohort geometry. This is the convergence
point of E20 (54% of positives are singletons), E21 (sparsity does not predict
AUC), E22 (62% of error is both-singleton) and this experiment.

---

## E31 — The discordant negatives look like fine specificity, not assay dropout

Plan cards S03, C09, C03 and C04 on one pair set,
`scripts/26_discordant_pair_anatomy.py`. E20 found negatives that are
single-residue CDR3 variants of validated positives and E27 relocated their
damage, but nothing had tested *why* they differ. The two candidates imply
opposite things: a mutation abolishing a peptide contact is real fine
specificity, so no sequence method can be blamed for missing it; a mutation at a
solvent-facing position points at the assay or the annotation.

**The test needs a null**, so each discordant mutation is compared against
**matched same-label pairs** -- positive/positive neighbours sharing one CDR3
exactly and differing by one or two residues on the other, within the same
epitope. Without that control, "mutations tend to be peptide-facing" would just
be a property of CDR3 mutation in general. 71 discordant pairs (103 mutated
sites) against 522 control pairs (767 sites).

| mutated-site property | discordant | control | stratified permutation p |
|---|---:|---:|---:|
| peptide contacts at that position | **6.53** | 2.90 | **0.0070** |
| change in peptide contacts between members | **3.37** | 1.27 | **0.0010** |
| position is peptide-facing at all | 0.340 | 0.231 | 0.0790 |
| volume change | 53.7 | 41.9 | (Mann-Whitney 0.023) |
| charge change | 0.27 | 0.23 | (Mann-Whitney 0.314) |

Permutation is stratified by epitope over 2,000 draws, so an effect produced
only by GILGFVFTL dominating the pool would wash out. Two things follow.

**The continuous effect survives and the binary one does not.** It is not that
discordant mutations are more often at a peptide-contacting position -- that
drops to p=0.079 once epitope composition is held fixed. It is that they sit at
positions carrying **2.3x more peptide contact** and produce **2.7x larger
changes** in it. The more precise claim is the quantitative one.

**Direction is consistent where there is data:** GILGFVFTL (65 sites),
GLCTLVAML (10) and NLVPMVATV (26) all agree; YLQPRTFLL has 2 sites and goes the
other way, TTDPSFLGRY has none. Non-independence is mild -- 71 discordant pairs
come from 87 unique TCRs, so most appear once, though one TCR appears in 14.

**Settled, with its limits.** This supports genuine fine specificity over assay
dropout: if these were dropouts, the mutated position would be random with
respect to the peptide interface and therefore indistinguishable from the
control. It is not. So the near-positive discordance is at least partly real
biology, which means sequence-neighbourhood methods have a **biologically
grounded ceiling** rather than merely a noise ceiling -- and that it cannot be
fixed by better sequence modelling, only by reading the interface, which E28
showed is too weak to use here.

Limits worth stating: these are predicted structures, so a predicted contact is
not a measured one; charge change shows no effect, so the simple
"salt-bridge lost" story is not supported; and this explains the ceiling without
moving the competition score at all.

---

## E32 — Six remaining sequence cards, all rejected (D04, D05, D10, S01, S02, S07)

`scripts/27_remaining_sequence_cards.py`, all judged by strict nested LOEO
against the submitted blend at 0.8498.

| card | result | verdict |
|---|---|---|
| D04 substitution cost on CDR3 | identity 0.8007, BLOSUM62 0.7999, BLOSUM45 0.8014 | the cost matrix does not matter |
| D05 multi-scale neighbourhood | three-scale average 0.8378 against single sigma=0.1 at 0.8396; blended +0.0007 | worse than one scale |
| S01 residual learnable information | best extra feature +0.0002, `iptm` **-0.0193** | third independent confirmation that the labels add nothing |
| S02 singleton rescue model | global **-0.0186**; within-singleton 0.7704 against the frozen score's 0.793 | a fitted singleton model is worse than the frozen score even on singletons |
| S07 abstention layer | close-neighbour group AUC 0.720, isolated group 0.778 | distance-based confidence does not map to ranking reliability |

**D10 is the one with lasting value: it says plainly what the model is.**
Binning by nearest-neighbour distance and scoring within each bin:

| nn1 band | n | validation rate | **within-band AUC** |
|---|---:|---:|---:|
| <0.06 | 310 | 0.823 | 0.791 |
| 0.06-0.12 | 189 | 0.921 | **0.553** |
| 0.12-0.25 | 221 | 0.805 | **0.523** |
| 0.25-0.5 | 526 | 0.432 | 0.690 |
| >0.5 | 534 | 0.139 | 0.659 |

Within-band discrimination is poor everywhere, 0.52 to 0.79, while the
validation rate falls 0.92 to 0.14 across bands. **The model is essentially a
monotone function of nearest-neighbour distance and barely discriminates inside
any narrow band of it.**

And the private epitopes sit in the far bands: 86.8% of CINGVCWTV, 75.1% of
LLWNGPMAV and 80.0% of LTDEMIAQY are beyond 0.25, where within-band AUC is
0.690 and 0.659. That is an independent and simpler confirmation of E30's
reprojection, reached without the pair-class decomposition.

Also useful for S07's actual purpose: the group with close neighbours has a
*lower* AUC (0.720) than the isolated group (0.778), because its validation rate
is already 0.86 and there is little left to order. So a distance-based
confidence layer would flag where the prior is extreme, not where the ranking is
unreliable -- worth knowing before carrying such a layer into IRIS.

---

## E33 — Pool and provenance cards, with two results that correct each other (R02, I02, R08, S05)

**I02 -- expansion does rescue positives preferentially, so it is not pure
volume.** When the pool gives an isolated TCR a neighbour inside the clustering
radius, that TCR is about three times more likely to be a validated positive:

| epitope | P(rescued \| positive) | P(rescued \| negative) | ratio |
|---|---:|---:|---:|
| GILGFVFTL | 0.134 | 0.041 | **3.2x** |
| GLCTLVAML | 0.055 | 0.036 | 1.5x |
| pooled | 0.047 | 0.015 | ~3x |

So expansion carries genuine evidence. But outside GILGFVFTL the absolute counts
are negligible -- NLVPMVATV, TTDPSFLGRY and YLQPRTFLL rescue two TCRs each -- so
this is entirely consistent with E27: the mechanism is real and there is nothing
to rescue on the private side.

**R02 -- and yet throttling kills the gain, so there is a volume component too.**

| pool policy | mean |
|---|---:|
| no pool | 0.8414 |
| all records | **0.8526** |
| one share per reference | 0.8504 |
| one per subject | 0.8428 |

Deduplicating by subject nearly erases the benefit. The gain therefore has both
an evidential component (I02) and a volume component (R02), which is why the
raw-record version of expansion looked better than it deserved in E27.

**R08 -- how much of "clustered" just means "easy to generate" is
epitope-specific.** Matching on terciles of paired Pgen crossed with CDR3-beta
length and re-measuring the clustering-to-validation gap:

| epitope | raw gap | gap within matched strata |
|---|---:|---:|
| **GILGFVFTL** | 0.400 | **0.067** |
| GLCTLVAML | 0.281 | 0.375 |
| NLVPMVATV | 0.658 | **0.711** |
| TTDPSFLGRY | 0.518 | 0.416 |
| YLQPRTFLL | 0.576 | 0.412 |

On GILGFVFTL the association nearly vanishes once generation probability and
length are matched -- its clusters are largely high-Pgen, easy-to-generate
motifs. On NLVPMVATV it is untouched. This is the fine-grained version of E7's
finding that Pgen is not orthogonal: **the degree of Pgen confounding is itself
epitope-specific**, which is one more reason a single global correction (E8)
could only hurt. Caveat: few TCRs per stratum, so these are noisy.

**S05 -- provenance weak supervision, dead for the third time.** Records per
clonotype 0.534, distinct references 0.557, distinct assay classes 0.514, and
every one of them takes only one or two distinct values within each private
epitope.

---

## E34 — Predicted contacts against crystallography (plan card D09)

E28 built contact fingerprints from Boltz predictions and E31 used them to argue
that discordant negatives differ at peptide-contacting positions. Both rest on
predicted interfaces, which E31 flagged as its main limitation. This checks them
against experiment: `scripts/28_real_structure_validation.py`, with
`scripts/29_fold_real.sh` folding the crystal sequences through the *same*
pipeline used for the competition.

The PDB yields 14 solved TCR-pMHC complexes for four competition epitopes
(GILGFVFTL 9, YLQPRTFLL 3, NLVPMVATV 2, GLCTLVAML 1), plus one out-of-training
anchor described below. CINGVCWTV and LLWNGPMAV do have crystal structures
(3mrg, 5n6b, 9sl0) but all are peptide-MHC only, with no TCR.

**Result on the solved complexes: the fingerprints are faithful.**

| statistic | median |
|---|---:|
| Spearman, contact counts per peptide position | 0.992 |
| Spearman, restricted to contact-bearing positions | **0.953** |
| Spearman, minimum distance per position | 0.983 |
| Jaccard of touched positions | 1.000 |
| total contacts, predicted against measured | 102 vs 99 |

**Three self-checks, two of which cut the claim down.**

*Zero inflation.* Of 20 fingerprint slots, 11.5 on average are zero in both
vectors -- the buried anchor positions no TCR can reach. Restricting to
contact-bearing positions drops the median from 0.992 to **0.953**. Still high,
but the headline number was partly agreement on what is structurally
impossible.

*Memorisation.* Every one of the 14 complexes was released between 2003 and
2022, so all of them predate Boltz-2. **This is a memorisation check, not a
prediction test.** One post-training complex exists -- 9x7r, released
2026-09-23, a TCR-pMHC complex on an unrelated peptide -- and on it accuracy
drops materially: contact-count Spearman **0.840**, distance Spearman 0.848,
Jaccard **0.700**, and total contacts overpredicted by 40% (146 against 104).
n=1 is an anecdote, not evidence, but it is the only out-of-distribution point
available and it points the expected way.

*Between-complex discrimination -- the check that matters most.* Reproducing
each structure individually would be useless if the model could not reproduce
how they *differ*, which is the only thing a fingerprint feature uses. Over the
36 pairs of GILGFVFTL complexes, Spearman between predicted and measured
pairwise fingerprint distance is **+0.948**, with comparable magnitudes (real L1
65.3 +/- 36.2, predicted 69.9 +/- 39.1). Within the memorised regime, the
representation does carry real inter-TCR interface variation.

**What this does and does not license.** The fingerprint machinery is validated
as a representation of real interfaces for TCRs the model has seen. Accuracy on
the competition's 2,636 TCRs, almost none of which are in the PDB, remains
unestablished, and the single OOD point suggests degradation.

The consequences are asymmetric:

- **E28's negative result is safe, and strengthened.** Prediction noise can only
  weaken a real signal, so "contact fingerprints add nothing" holds more firmly
  if the predictions are noisier out of distribution than here.
- **E31's positive result is exposed, but partly defended by its design.** It
  compared discordant pairs against matched same-label controls drawn from the
  *same* prediction pipeline, so a uniform prediction bias affects both arms
  equally. Producing that effect from prediction error alone would require the
  error to correlate with the mutation itself. That is a real defence, not
  immunity.

A code review of the comparison also caught one flaw of my own: the first
version pooled contact counts and minimum distances into a single Spearman. The
two blocks are inversely related and occupy different value ranges, so they
separate in the rank order of both vectors regardless of whether the
per-position pattern agrees, inflating the correlation. Fixed to report the
blocks separately. Two further caveats recorded rather than hidden: chain "A" in
these fingerprints is whichever TCR chain contacts the peptide more, not
necessarily the alpha chain, so chain-specific claims are out of scope; and the
RCSB search for recent entries returned computed-structure models despite
requesting experimental content, so the out-of-distribution set was filtered
against the sequence database rather than trusted from the query.

---

## E35 — Reference-database version, source and leakage (plan card R06)

Three VDJdb releases downloaded -- 2024-06-13 (32,183 paired human rows),
2025-02-21 (35,672) and 2026-06-03 (86,221) -- and the pool measurement from
E27/E33 rerun on each with the prediction function fixed.

**Leakage audit, which had not been done and should have been.** If VDJdb had
absorbed records from the Schumacher re-screening, the expansion pool would
contain the answer key. It has not: searching every human reference for the
Immunity 2026 DOI, the 2025 bioRxiv preprint, "TCRvdb", the Schumacher lab repo
and "T-RAP" returns **zero records each**. The pool is independent of the
labels.

**The competition's screening set most resembles a ~2025-02 snapshot.** Over
the five epitopes the organisers appear to have taken whole rather than
subsampled, the median absolute difference between VDJdb clonotype count and
test-set size is 4 for 2025-02-21, against 10 for 2024-06-13 and 25 for
2026-06-03 -- and CINGVCWTV matches exactly, 243 against 243. Suggestive, not
proof: VDJdb removes records between releases too, and the organisers applied
their own QC for synthesis failure and ambiguous calls. It does explain
GILGFVFTL (3,282 available, 670 screened) and GLCTLVAML (487, 192) as
deliberate subsampling.

**The pool gain is not resolvably version-dependent.**

| reference pool | mean | gain |
|---|---:|---:|
| none | 0.8414 | -- |
| VDJdb 2024-06-13 | 0.8466 | +0.0052 |
| VDJdb 2025-02-21 | 0.8513 | +0.0099 |
| VDJdb 2026-06-03 | 0.8526 | +0.0112 |

The trend looks monotone, and the self-check says it is not interpretable. Paired
across the same five epitopes: 2026 minus 2024 is +0.0060 (se 0.0048), 2026
minus 2025 is +0.0012 (se 0.0023), 2025 minus 2024 is +0.0048 (se 0.0064). **No
pairwise difference reaches two standard errors.** The per-epitope detail shows
why: the spread is carried almost entirely by GLCTLVAML (+0.005, +0.035, +0.030)
while GILGFVFTL, NLVPMVATV, TTDPSFLGRY and YLQPRTFLL are flat to within 0.003.
A 2.7-fold increase in database size produces no resolvable change in the gain,
which is itself the result.

It also falsifies a claim I made in narration before measuring -- that the gain
would be smaller with the snapshot the organisers actually used. It is 0.0012
smaller, inside the noise.

**The durable finding is a count, not an AUC, so noise does not apply.**
Available pool multiplier per private epitope across the three releases:

| epitope | 2024-06 | 2025-02 | 2026-06 |
|---|---:|---:|---:|
| CINGVCWTV | 1.01 | 1.00 | 1.05 |
| LTDEMIAQY | 1.04 | 1.02 | 1.04 |
| LLWNGPMAV | 1.08 | 1.03 | **1.83** |

CINGVCWTV and LTDEMIAQY have gained essentially no new paired TCRs in eighteen
months, so **their exhaustion is a durable property of the epitope rather than
an artefact of which snapshot is compared against** -- which makes E27's and
E33's conclusion version-robust for two of the three deciding epitopes.
LLWNGPMAV is the exception and its surplus is recent, added after the
organisers' snapshot.

---

## E36 — Target against anchor-matched decoy peptides (plan cards D08, C07)

E14 inferred that Boltz confidence mostly reports how canonical a TCR is, from
the fact that `iptm_mhc_pep` -- an interface containing no TCR at all -- scored
0.647. This tests it directly. Folding the same TCR against its annotated
peptide and against decoys and subtracting cancels everything that is a property
of the TCR alone, leaving only what is peptide-specific.

**Decoys are anchor-matched by class, not by identity.** HLA-A*02:01 constrains
an aliphatic P2 and C-terminus; HLA-A*01:01 an acidic P3 and aromatic
C-terminus. A first draft of the script asserted identical residues and failed
immediately -- no reported A*02:01 nonamer shares GILGFVFTL's P2 isoleucine
*and* its C-terminal leucine -- which also exposed two decoys I had chosen with
a serine at P3, where A*01:01 needs acid. 300 TCRs (30 positive and 30 negative
per public epitope) x 3 decoys = 900 folds; targets were already in the main run.
`scripts/30_control_peptide_inputs.py`, `scripts/31_fold_decoys.sh`.

**Result: the margin carries nothing.**

| quantity | target only | target minus decoy |
|---|---:|---:|
| ipTM | **0.7387** | 0.5362 |
| confidence score | 0.7189 | 0.5398 |
| TCR alpha-peptide ipTM | **0.7429** | 0.5331 |
| TCR beta-peptide ipTM | 0.7329 | 0.5251 |
| MHC-peptide ipTM | 0.6500 | 0.5147 |
| complex pLDDT | 0.6313 | 0.5247 |

Target-only confidence predicts validation at 0.63-0.74. Remove the
TCR-generic component by differencing against a matched decoy and **every
margin collapses to chance**. The structural signal is a property of the TCR and
its fold, not of the peptide it is paired with -- E14's inference, now
demonstrated rather than deduced.

**The raw magnitudes say the same thing more bluntly.** Mean target ipTM 0.9152
against decoy 0.9091: a margin of **+0.0061**. Boltz assigns a TCR's true target
and an unrelated anchor-matched epitope essentially the same interface
confidence. The MHC-peptide margin is +0.0001, which is at least a good internal
check -- both peptides occupy the groove equally, as anchor matching intended.

**Self-checks, which turn this from "noisy" into "far too small".**

*Are the decoy structures actually different?* Yes. The spread of decoy ipTM
across the three decoys of one TCR is 0.036, comparable to the 0.062 spread
across different TCRs. Boltz is responding to the peptide -- just not in a way
that tracks the real target.

*Signal to noise.* The margin's mean (+0.0061) is a tenth of its between-TCR
spread (0.0602), and decoy-choice noise alone (0.036) is **5.9 times** the
margin.

*Is any one decoy carrying it?* No. Per-decoy margins run from **-0.0128 to
+0.0247** and two decoys are negative, both A*01:01 ten-mers -- so for
TTDPSFLGRY, two of three decoys score *higher* than the true target. Per-decoy
AUCs scatter from 0.488 to 0.676 around chance.

*How big would it need to be?* Reaching AUC 0.70 would take a margin near
0.042; the observed 0.0061 is **14% of that**. So this is not a power problem
that more decoys would fix -- the effect is roughly seven times too small.

**Phase 2 was run after all, and it was worth the seven GPU-hours.** I had
staged the private extension behind a positive phase-1 result and then skipped
it, on the argument that a feature failing on all five labelled epitopes can
never clear the pre-registered threshold -- since that threshold is only
measurable where labels exist. The argument stands, but it answers a different
question from the one phase 2 asks: *is the failure universal, or specific to
the epitopes that happen to carry labels?* Without the data that could only be
asserted. All 902 TCRs x 3 decoys = 2,706 folds completed.

| | target ipTM | decoy ipTM | margin | sd | fraction positive |
|---|---:|---:|---:|---:|---:|
| GILGFVFTL | 0.9317 | 0.9174 | +0.0143 | 0.0512 | 0.650 |
| GLCTLVAML | 0.9284 | 0.9183 | +0.0101 | 0.0543 | 0.617 |
| NLVPMVATV | 0.9154 | 0.8988 | +0.0166 | 0.0499 | 0.750 |
| TTDPSFLGRY | 0.9024 | 0.9062 | **-0.0038** | 0.0693 | 0.533 |
| YLQPRTFLL | 0.8979 | 0.9048 | **-0.0069** | 0.0711 | 0.617 |
| CINGVCWTV (private) | 0.9022 | 0.8973 | +0.0049 | 0.0609 | 0.617 |
| LLWNGPMAV (private) | 0.9525 | 0.9238 | +0.0287 | 0.0242 | 0.939 |
| LTDEMIAQY (private) | 0.9233 | 0.9035 | +0.0199 | 0.0700 | 0.731 |

**The failure is universal.** Mean absolute margin is 0.0406 on the labelled
epitopes and 0.0394 on the private ones -- indistinguishable -- so there is no
reason to think the quantity behaves differently where it cannot be scored. It
does have within-epitope variance everywhere, so it passes E3's formal
requirement; the variance is simply noise, as the labelled-side AUC of 0.536
establishes.

Two public epitopes carry **negative** mean margins: on TTDPSFLGRY and
YLQPRTFLL the model prefers the decoys on average. That is about as direct a
statement of absent peptide specificity as this setup can produce.

**Settled.** This closes the structural line with its mechanism identified:
Boltz interface confidence contains no peptide-specific information in this
setting, so every structural feature tried (E14, E19, E28, E36) was reading the
same TCR-canonicality signal that density already reads.

---

## E37 — Does the E31 interface result survive relaxation? (plan cards C10, C04)

E31 is the project's one positive mechanistic finding and it was measured on raw
Boltz coordinates, which carry clashes and unrelaxed side chains. The plan's own
standard is that an effect unstable to pose or parameters warrants no
conclusion. 69 complexes relaxed under one fixed protocol -- PDBFixer for
terminal atoms and hydrogens at pH 7, Amber14 with GBn2 implicit solvent, L-BFGS
minimisation at 200 and 1,000 iterations. `scripts/32_minimisation_robustness.py`.

**Descriptor stability differs sharply, which matters for which ones to trust.**

| descriptor | median relative change, 0 to 1,000 steps | Spearman(raw, relaxed) |
|---|---:|---:|
| buried interface area | **2.1%** | **+0.963** |
| peptide contacts | 8.5% | +0.793 |
| hydrogen bonds | **33.3%** | +0.748 |

Buried area is robust; hydrogen-bond counts are not, and 16% of complexes move
their contact count by more than 20%. So a conclusion resting on hydrogen-bond
counts would not be safe, while one resting on contacts needed checking -- which
is what E31's does.

**The E31 effect survives.** Measured on the same 20 discordant and 20 control
pairs (31 and 32 mutated sites), before and after relaxation:

| quantity | raw | minimised |
|---|---|---|
| peptide contacts at the mutated site | disc **8.61** vs ctl 3.97, p=0.026 | disc **7.16** vs ctl 3.53, p=0.045 |
| change in contacts between members | 4.61 vs 0.97, p=0.028 | 4.39 vs 1.34, p=0.041 |
| site is peptide-facing | 0.419 vs 0.156, p=0.022 | 0.419 vs 0.188, p=0.048 |

Effect sizes shrink modestly -- the contact ratio goes 2.17x to 2.03x -- and all
three remain significant with direction unchanged. **E31 is not an artefact of
unrelaxed coordinates.** The p-values here are weaker than E31's only because
this subset has 31 sites against E31's 103; the comparison that matters is raw
against minimised on identical pairs.

What this deliberately does not produce is an interface binding energy. Doing
that properly needs separated-state topologies whose termini and atom counts
differ from the complex, so the subtraction would not be comparable, and a
single-point interaction energy on an unrelaxed separated state is not a binding
free energy whatever it is called.

**Two code-review findings from my own work.** The first version of this script
stored only whole-interface descriptors, which cannot test E31 at all: E31
measures contact at the *mutated CDR3 position*, and both members of a pair are
near-identical TCRs whose whole interfaces are bound to look alike -- duly
giving p=0.971. It had to be rerun saving relaxed coordinates. And OpenMM
renames chains on write, M and P becoming C and D, so the position-specific
analysis found no peptide chain until it was changed to identify the peptide by
residue count rather than by chain letter. The hydrogen filter was verified
separately: OpenMM does write the element column, so the 3,192 added hydrogens
per complex are correctly excluded.

---

## E38 — Negative controls and the pool x structure factorial (plan cards R04, R05)

Both cards I had skipped on my own judgement. Running them corrected me twice.

**R04: the pool gain is real information, not volume.** E27 audited whether the
*query* finds itself in the pool, which is a different worry from whether adding
any sizeable batch of sequences raises density regardless of content. Three
size-matched controls, all 2,684 added sequences, the same count as the real
pool:

| reference added | mean | gain |
|---|---:|---:|
| none | 0.8414 | -- |
| **real same-epitope** | **0.8526** | **+0.0112** |
| other-epitope, shape-matched on CDR3-beta length | 0.8333 | **-0.0081** |
| other-epitope, unmatched | 0.8350 | -0.0064 |
| OLGA synthetic repertoire | 0.8381 | -0.0063 |

**All three controls actively hurt.** Adding the same number of irrelevant
sequences does not merely fail to help, it degrades the density estimate. So the
expansion gain is antigen-specific information, which together with I02
(positives rescued about three times as often) and R02 (subject-level throttling
nearly erases it) completes the picture: a real mechanism, carried by multiple
records from the same subjects, that the private epitopes cannot access because
only ~4% more clonotypes exist for CINGVCWTV and LTDEMIAQY (E27).

**R05: pool and structure are additive, and my reason for skipping it was
wrong.** I had argued the factorial would be degenerate because the pool arm
measures near zero -- but that was the *private-realizable* figure; on the
public epitopes the pool gives +0.0112, so the design is informative.

| | original pool | expanded pool |
|---|---:|---:|
| no structure | 0.8414 | 0.8526 |
| structure @0.3 | 0.8487 | **0.8583** |

| effect | value |
|---|---:|
| structure without pool (B-A) | +0.0073 |
| structure with pool (D-C) | +0.0058 |
| pool without structure (C-A) | +0.0112 |
| pool with structure (D-B) | +0.0096 |
| **interaction** | **-0.0015** |

The interaction is essentially zero: the two contribute independently, neither
redundant nor synergistic. One limit on reading that as good news -- the
additivity is a public-side fact. The pool arm does not transfer to the private
epitopes (+0.0004, E27) while the structure arm does (+0.0070 nested, E21), so
on the deciding epitopes only one arm of this factorial survives.

---

## Where this leaves things

Thirty-eight experiments, one consistent picture: **the extractable signal is
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
