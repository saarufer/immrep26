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

### E37 revision — the three rows are one finding, and it is positional

Found while preparing C11, which needed the same mutated-site residues. Measured
on E37's own relaxed coordinates, 40 pairs and 126 mutated-site measurements,
`scripts/35_md_analysis.py --static`.

E37 read one descriptor, a 4.5 A contact count. Reading a second descriptor of
the same physical quantity -- the minimum heavy-atom distance from the mutated
site to the peptide -- splits the result apart:

| descriptor | disc | ctl | p |
|---|---:|---:|---:|
| peptide contacts at the site (count) | **3.581** | 1.766 | **0.017** |
| site-to-peptide minimum distance (A) | 6.412 | 6.354 | 0.899 |

**This is not an outlier artefact, and the first version of this audit "checked"
that wrongly.** It rank-transformed the count and re-ran, which returns the
identical p-value by construction, because Mann-Whitney is already a rank test
and is invariant to monotone transforms. The check verified nothing. The real
reason the two disagree is that **77% of mutated sites have no peptide contact
at all**, so the count is three quarters ties at zero and its rank test is in
effect a binary comparison, while the distance uses all 126 sites.

Testing that binary reading directly, and then conditioning on it:

| quantity | disc | ctl | p |
|---|---:|---:|---:|
| **site touches the peptide at all** | **0.323** | **0.141** | **0.020** (Fisher) |
| distance given the site touches (A) | 3.504 | 3.319 | 0.195 |
| contacts given the site touches | 11.100 | 12.556 | 0.981 |

**The whole effect is the binary one.** Discordant pairs differ at a
peptide-facing position 2.29x as often -- which is numerically E31's headline
"2.3x more peptide contact". Conditional on the site facing the peptide there is
no difference at all: the discordant site is marginally *farther* and carries
marginally *fewer* contacts, both far from significance.

Two consequences.

**E31 should be restated.** Not "discordant negatives mutate at positions
carrying more peptide contact", which implies a graded energetic difference, but
"discordant pairs are about twice as likely to differ at a peptide-facing
position". The mechanism is positional -- which residue the pair differs at --
not energetic. That is a sharper claim and the data support it better.

**E37 overstated its own evidence.** Its three rows -- contact count, change in
contacts between members, site is peptide-facing -- are not three independent
confirmations. With 77% ties at zero the count row and the facing row are
substantially the same test, and E37's "all three remain significant with
direction unchanged" treated redundant tests as corroborating. The conclusion
that E31 survives relaxation stands; the weight of evidence behind it is one
test, not three.

This also sets what C11 can add. The quantity that survived is binary facing, so
the dynamic question is whether a site that faces the peptide keeps facing it --
the fraction of frames in contact. Power is the limit: only 29 of 126 sites face
the peptide at all, 20 discordant against 9 control, so a null MD result will be
uninformative rather than evidence of absence, and that is stated before the run
rather than after it.

---

### E37 correction — three numbers in its descriptor table are not reproducible

Found by the data audit (`scripts/36_data_audit.py`), which flagged that
`scripts/32` declared `reports/10_minimisation_robustness.md` and never wrote
it, so E37's stability table existed only here. Rebuilding it from
`data/interim/minimisation_robustness.csv` (69 complexes x 3 stages, no missing
cells) with `scripts/32_minimisation_robustness.py --report`:

| descriptor | logged | reproducible | logged rho | reproducible rho |
|---|---:|---:|---:|---:|
| buried interface area | 2.1% | **1.63%** | +0.963 | +0.960 |
| peptide contacts | 8.5% | **9.21%** | +0.793 | +0.793 |
| hydrogen bonds | 33.3% | 33.33% | +0.748 | **+0.788** |

The contacts rho matches exactly and the hydrogen-bond median matches exactly,
which fixes the stage pair as 0 to 1,000 and the denominator as the raw value;
no variant of stage pair, denominator, mean against median, or signed against
absolute reproduces 8.5% or 2.1%. The likely cause is the same one that has
bitten twice already: the table was computed in a heredoc that was never saved,
probably while the background minimisation was still filling the CSV.

**The conclusion is unchanged and if anything strengthened** -- buried area is
the robust descriptor at 1.63%, hydrogen bonds are unusable at 33%, contacts sit
between at 9.21% with 17% of complexes moving more than a fifth. Only the three
numbers above should be cited, and a sixth of E37's reported figures were wrong
in the second decimal because the computation was not in version control.

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

## E39 — Is the interface difference dynamically stable? (plan card C11)

The last structural card, and the one I twice called infeasible for the wrong
reason. 69 relaxed complexes from 40 pairs, three velocity seeds each, 2 ns at
310 K with a 4 fs step and hydrogen-mass repartitioning, implicit GBn2, sampled
every 20 ps. 207 trajectories, 414 ns, no failures.
`scripts/34_md_interface_stability.py`, `scripts/35_md_analysis.py --md`.

**The setup holds the pose.** Interface contact rises a median 7.1% over each
trajectory and 20 of 207 trajectories lose more than a quarter of their
contacts. Side chains settle into contact on heating from the minimised state;
nothing is dissociating, so the trajectories are sampling the pose rather than
losing it.

**Every group comparison at the mutated site is null.**

| quantity | disc | ctl | p |
|---|---:|---:|---:|
| contacts, time-averaged | 2.456 | 1.332 | 0.159 |
| distance, time-averaged (A) | 6.817 | 6.803 | 0.545 |
| fraction of frames in contact | 0.288 | 0.154 | 0.149 |
| **touches the peptide at any point** | **0.597** | **0.594** | **1.000** |

The last row is the one to read against the static result. On the relaxed
snapshot the binary facing rate was 0.323 against 0.141 at p=0.020; given 2 ns
and three replicas, both groups reach the peptide at the same rate. The loose
criterion washes the difference out entirely.

**MD adds nothing beyond the static facing rate, and that is an identity, not an
impression.** Conditioning on whether the site faced the peptide in the relaxed
pose, both groups behave the same under dynamics: among the 29 facing sites the
contact fraction is 0.758 against 0.708 (p=0.943), and among the 97 non-facing
sites it is 0.064 against 0.063 (p=0.579). Decomposing the overall fraction:

```
discordant  0.323 x 0.758  +  0.677 x 0.064  =  0.288   (observed 0.288)
control     0.141 x 0.708  +  0.859 x 0.063  =  0.154   (observed 0.154)
```

The conditional terms are group-independent to three decimals. So the dynamic
contact fraction *is* the static facing rate carried through a persistence
factor that does not distinguish the groups. The trajectories neither amplify
nor erode E31; they reprice it at exactly what the snapshot already said.

**The direction is not consistent across epitopes**, which was the third
pre-planned self-check and it fails:

| epitope | disc | ctl | sites (disc/ctl) |
|---|---:|---:|---:|
| GILGFVFTL | 1.161 | 0.641 | 34/46 |
| NLVPMVATV | 5.322 | 0.187 | 18/8 |
| GLCTLVAML | 2.835 | **9.357** | 6/2 |
| YLQPRTFLL | **0.000** | 4.442 | 4/8 |

Two epitopes in the claimed direction, two reversed, on 2 to 8 sites in the
small cells. No epitope-level conclusion is available from this.

**Verdict, against the limit stated before the run.** I pre-registered that only
29 of 126 sites face the peptide at all, 20 discordant against 9 control, so a
null result here would be uninformative rather than evidence of absence. That is
what happened. C11 does not refute E31 and does not support it. What it does
establish is narrower and worth having: the effect is entirely contained in the
static facing rate, so no further structural dynamics will add to it, and this
line is closed on evidence rather than on judgement.

**Three code-review findings in my own work, two of which had manufactured a
result.**

The replica-spread guard the plan requires passed everything by construction. It
compared the group-mean difference against the median per-site replica sd, but
most sites never touch the peptide in any replica, so their sd is trivially zero
and the median with it. Rebuilt as the spread of the group difference when each
replica is used alone, it correctly reports the distance effect as no result --
the sign flips, -0.025, +0.200, -0.134 -- while the contact fraction agrees
across all three. Found by running the analysis against partial trajectories at
83 of 207 rather than waiting, which is the only reason it was fixed before the
write-up rather than after.

Joining the static measurements to the MD table inflated 126 site rows to 150,
because the same TCR appears in several pairs and the same residue can be the
mutated site in more than one. Before the fix this showed the binary facing rate
*reversing* under dynamics, 0.481 against 0.606, and gave formation among
non-facing sites p=0.047. Both were artefacts of the duplicated join; corrected,
they are p=1.000 and p=0.579.

The groups are not disjoint at the level of sites. Three of 117 distinct sites
are a discordant site in one pair and a control site in another, since a mutated
position can be discordant with one partner and same-label with a different one.
The effect is small here but the point applies to E37's tests equally, and
neither entry had noted it.


## E40 — Does E30's private projection transfer? It does not (reviewer's P0)

A reviewer identified a problem with E30 and the problem is real. Checking it
turned up two more, and the conclusion is that **E30's 0.799 is withdrawn**.
`scripts/37_private_projection_validation.py`.

**First: E30 was computed in a heredoc that was never saved.** The log entry
cites no script. So the project's single most consequential number -- the one
written into long-term memory -- rested on a computation absent from version
control. This is the fourth time that pattern has caused a problem.

**The composition is not label-free, which is the reviewer's point.** A class
such as `pos_clu_neg_sing` needs to know which pair member is positive. With
a = P(clustered | positive), b = P(clustered | negative), pi = P(positive):

    pos_clu_neg_sing = a(1-b)    pos_sing_neg_clu = (1-a)b
    both_singleton   = (1-a)(1-b)    clustered-both = ab

and the only label-free observable is q = pi*a + (1-pi)*b. One equation, two
unknowns. **And pi is a third unknown the review did not name**: it runs 0.287
to 0.663 across the public epitopes. Splitting the both-clustered mass needs a
fourth parameter. So the private composition is assumed, not measured.

| epitope | q | pi | a | b | lambda=a/b |
|---|---:|---:|---:|---:|---:|
| GILGFVFTL | 0.476 | 0.663 | 0.601 | 0.149 | 4.02 |
| GLCTLVAML | 0.552 | 0.663 | 0.645 | 0.339 | 1.90 |
| NLVPMVATV | 0.183 | 0.287 | 0.571 | 0.058 | 9.93 |
| TTDPSFLGRY | 0.047 | 0.497 | 0.082 | 0.000 | undefined |
| YLQPRTFLL | 0.184 | 0.414 | 0.411 | 0.040 | 10.24 |

lambda spans 1.9 to 10.2 where it is defined at all.

**The reconstruction E30 cited as validation is an identity.** Recombining each
epitope's own class error rates with its own class shares returns its own AUC;
measured residual **1.1e-16**, machine epsilon. E30 reported "within 0.0003" and
read it as evidence the arithmetic transfers. It is evidence the arithmetic is
arithmetic.

**A basis mismatch nobody had noticed.** `scripts/20` clusters only the
*labelled* rows; `scripts/25`, whose family rates E30 reasoned from, clusters
the *whole* cohort. E30 combined error rates from one basis with compositions
from the other. Everything here is recomputed whole-cohort, the only basis a
private epitope can supply.

**P0: leave-one-epitope-out transfer.** Estimate every transferred quantity on
four epitopes, predict the fifth from its label-free q alone.

| held out | actual | full transfer | oracle composition | oracle errors |
|---|---:|---:|---:|---:|
| GILGFVFTL | 0.869 | 0.954 | 0.876 | 0.938 |
| GLCTLVAML | 0.741 | **0.954** | 0.833 | 0.929 |
| NLVPMVATV | 0.865 | 0.856 | 0.902 | 0.782 |
| TTDPSFLGRY | 0.840 | 0.847 | 0.852 | -- |
| YLQPRTFLL | 0.933 | **0.830** | 0.864 | 0.912 |
| | | **MAE 0.083** | MAE 0.043 | MAE 0.090 |

**The method cannot predict a held-out public epitope.** Its out-of-sample MAE
of 0.083 is nearly twice the 0.047 bias it claimed to have detected. GLCTLVAML
is overpredicted by 0.213. Splitting the two transfers shows the composition is
the weaker half (0.090 against 0.043), which is exactly the identifiability
problem the reviewer identified, now measured.

**And done honestly, the projection does not give 0.80.** Sweeping the
unobservables over the public range gives a private macro of **0.829 to 0.857**,
median 0.846 -- indistinguishable from LOEO. The mean transfer bias is +0.038
with a standard error of about 0.054 across five folds, so it is not
distinguishable from zero and bias-correcting is not justified.

**A factual error in E30, which the reviewer caught.** It said all five public
epitopes are family-richer than the private three. False: TTDPSFLGRY at 0.047
sits below LTDEMIAQY at 0.146. Correct statement: two of three private epitopes
are below every public epitope, and the third is inside the public range.

**Revised position, replacing E30's.**

| number | what it is |
|---|---|
| **0.846** | the estimate: nested LOEO over five public epitopes |
| 0.799 | **withdrawn** -- a projection whose assumptions were never recorded and which fails its own out-of-sample test |
| 0.67-0.91 | retained, but as a **scenario range** over within-both-singleton separability, not a confidence interval |

The real finding of E30 survives and is worth keeping: the private cohorts
differ from the public ones on a label-free observable, and LOEO cannot see it.
What does not survive is any number attached to that difference. **With five
epitopes the size of the shift is not measurable**, which is itself the result
and agrees with E21.

---

## E41 — Is E31 robust to non-independence? Partly (reviewer's P1)

E31 treated each mutated site as a draw. `scripts/38_e31_family_robustness.py`
re-tests with the clonotype family as the sampling unit -- single-linkage at
SCEPTR cosine 0.12 over the whole cohort -- using a cluster bootstrap over
families and a permutation that shuffles the group label between families.

**The effective sample size was overstated about tenfold.**

| group | sites | pairs | TCRs | **families** |
|---|---:|---:|---:|---:|
| control | 767 | 522 | 302 | **79** |
| discordant | 103 | 71 | 87 | **33** |

The three largest families contribute 129, 92 and 73 sites on their own, and 15
families appear in both groups.

**On E31's full raw-structure set:**

| quantity | disc | ctl | effect | family 95% CI | family p | E31's site-level p |
|---|---:|---:|---:|---|---:|---:|
| peptide contacts at site | 6.534 | 2.905 | +3.629 | [+0.43, +7.83] | **0.043** | 0.0070 |
| site faces the peptide | 0.340 | 0.231 | +0.109 | [-0.03, +0.27] | **0.236** | 0.0790 |

**On the 40-pair relaxed subset E37 and E39 used** (126 sites, 30 families):

| quantity | effect | family 95% CI |
|---|---:|---|
| peptide contacts at site | +1.815 | **[-1.16, +5.19]** crosses zero |
| site faces the peptide | +0.182 | [+0.04, +0.38] |

**The two samples disagree about which reading survives.** All four effects are
positive, but the contact count clears family-level resampling only on the full
set and the binary facing rate only on the relaxed subset. **No single reading
survives in both.** That directly undercuts the E37 revision's claim that "the
whole effect is the binary one" -- that was read off the relaxed subset alone,
where the binary happens to be the surviving reading, and it reverses on the
larger sample.

**Leave-one-epitope-out does not reverse the direction**, which is the one thing
that holds: the contact effect runs 2.35 to 4.60 across the five drops against
3.63 on the full data. Direction is stable; magnitude is not pinned.

**Revised wording for E31**, close to the reviewer's: *near-identical TCR pairs
with discordant functional labels show an enrichment of substitutions at
predicted peptide-facing positions relative to concordant pairs. The direction
is consistent across epitopes and across both structural samples. With the
clonotype family as the sampling unit the enrichment is significant in one
reading per sample and the two samples disagree on which, so the effect size is
not established and causality is not demonstrated.*

---

## Corrections to earlier entries, from the same review

**E23's ceiling was mislabelled.** It said fixing the adversarial class
completely is "worth at most +0.0081". That is what the best weighting in the
tested family achieves; it is not the class's ceiling. `pos_sing_neg_clu`
contributes a mean loss of **0.0318** (reports/06), so eliminating its error
entirely is worth **+0.032**, four times the quoted figure. The conclusion --
that no per-TCR weighting can capture it, because the trade against the easy
class is about 13:1 -- is unaffected, but the two numbers answer different
questions and the entry conflated them.

**E20's "zero fitted parameters, so nothing to overfit" is too strong.** There
are no fitted coefficients, but feature and configuration selection are still
choices made against the public epitopes, and E15 priced exactly that at 0.0048.

**E21's "cost of 0.3 is bounded at 0.0002" applies to the five measured
epitopes**, not in general.

**E26's germline-confounding reading is one of two live explanations.** V-gene
usage does track study and donor, but germline-encoded CDR1 and CDR2 residues
genuinely contact pMHC, so a V-region signal is not automatically confounding.
The present data cannot separate them. A cheap test exists and has not been run:
split the CDR3-masked distance into framework-only and CDR1/2-only. Signal
concentrated in framework, which contacts nothing, would be confounding; signal
in CDR1/2 would be compatible with real recognition.

**E39's "no further structural dynamics will add to it" overreaches.** What was
shown is that 2 ns of implicit-solvent dynamics around the predicted pose adds
nothing beyond the static facing rate. Longer timescales, explicit solvent and
larger rearrangements were not sampled.


## E42 — Framework or CDR1/2? The last low-cost exploration (reviewer's point on E26)

E26 read the whole-V-domain alignment signal as germline confounding. The
reviewer correctly objected that the reading is not forced: CDR1 and CDR2 are
germline-encoded **and** genuinely contact peptide-MHC, so a V-region signal can
be real recognition. The two were never separated. Separating them is cheap, and
this is the experiment the review implied without proposing.
`scripts/39_v_region_ablation.py`.

**Design.** Split the CDR3-masked region into the part that cannot touch the
interface and the part that can. CDR1, CDR2 and CDR2.5 come from tcrdist3's
IMGT-derived reference keyed on the V gene, so no numbering scheme is invented:
0 genes absent, CDR1 located verbatim in 99.5% of chains, CDR2 98.0%, CDR2.5
99.98%. E26's own split was another unsaved heredoc -- the fifth -- so all four
regions are recomputed, and since its bandwidth is lost the comparison runs over
a common sigma grid and rests on ordering rather than on an absolute value.

| region | mean alpha length | sigma=0.02 | 0.05 | 0.10 | 0.20 | best |
|---|---:|---:|---:|---:|---:|---:|
| full V domain | 111 | 0.8243 | 0.8304 | **0.8354** | 0.8343 | 0.8354 |
| CDR3 masked | 98 | 0.8199 | 0.8260 | **0.8276** | 0.8264 | 0.8276 |
| **framework only** | 79 | 0.8154 | 0.8246 | **0.8270** | 0.8223 | **0.8270** |
| loops only (CDR1+2+2.5) | 18 | 0.8074 | 0.8092 | 0.8122 | **0.8173** | 0.8173 |

Within-epitope Spearman against the CDR3-masked ranking: **framework +0.994**,
loops +0.970.

**Framework carries essentially all of it.** 0.8270 of 0.8276, at a rank
agreement of +0.994, from a region that makes **no contact with peptide or
MHC**. The loops, which do contact both, recover less on their own and agree
less with the combined ranking.

**So the reviewer's alternative is tested and rejected, and E26's reading is now
demonstrated rather than inferred.** The ranking is not reading recognition
chemistry; it is reading V-gene identity through the scaffold. This does not say
CDR1/CDR2 are biologically unimportant -- they obviously are important -- it says
*this feature's ranking* is explained by a region that cannot participate in
recognition at all.

Recomputation also re-checks E26: full 0.8354 here against its 0.8379, CDR3
masked 0.8276 against 0.8287. The ordering reproduces; the small absolute
differences are the lost bandwidth.

---

## Predictor frozen — 2026-10-25

The competition predictor is frozen as of this entry. No further feature search,
no further configuration selection, no further submissions.

| | |
|---|---|
| submitted file | `submissions/density3_struct03.csv` |
| sha256 (first 16) | `523a4692df02467b` |
| built by | `scripts/18_final_submission.py` |
| model | rank average of SCEPTR `default` + `small` kernel density (cosine sigma=0.1) and tcrdist3 alpha+beta (sigma=20), blended 0.7/0.3 with Boltz-2 `pae_tcr_pep_min` |
| fitted parameters | none |
| public score | 0.8498 (leaderboard 0.849) |
| private estimate | **0.846**, downside risk real but unquantified; 0.67-0.91 is a scenario range, not an interval |
| deadline | 2026-11-02 15:00 UTC |

**Why freeze here rather than keep searching.** Every remaining avenue is closed
on evidence rather than on judgement: the signal saturates at 0.83-0.84 in every
representation (E6, E16, E17, E18), supervised routes lose to variance on five
folds (E10, E13, E32), pool expansion has no measurable gain on this
configuration (E27, E33, E35, E38), and the structural line is closed with its
mechanism identified (E36) and its one mechanistic lead reduced to an
unestablished effect size (E41). E15 priced further configuration selection at
-0.0048, so additional search has negative expected value.

**What remains open is not a modelling question.** The molecular mechanism --
whether substitutions at peptide-facing positions cause the functional
discordance -- cannot be settled by more computation on this dataset. E41 showed
the statistical evidence does not survive consistently once the clonotype family
is the sampling unit, and E34 showed the structures it rests on are validated
only in a memorised regime. That question needs independent data or functional
experiments, which is where it is left.

**One item outstanding, and it is not a computation.** Verify in the Kaggle UI
that the selected final submission is `density3_struct03.csv` and not the
earlier density-only file.


## Pipeline delivered — 2026-10-26

Freezing the predictor closed one half of the brief; the other half -- a
pipeline someone else can run -- was not actually delivered. Before today the
repository had **no environment specification of any kind**: no pyproject, no
requirements file, no lockfile, 41 scripts with no entry point and no tests.
Nothing here was reproducible by anyone but this machine.

**Dependency manifest.** `pyproject.toml`, pinned to what produced the recorded
results. The 14 direct dependencies were derived by AST-scanning every import
in `scripts/` and `src/`, not by `pip freeze`, so the list is the project's
actual surface rather than its transitive closure. Boltz-2 and the Kaggle CLI
are optional extras; OpenMM is deliberately **absent**, because the PyPI wheel
ships no CUDA plugin and the working build has to come from conda-forge.

**Reproducibility, stated as a checkable claim rather than a hope.**
`make submission` rebuilds the entered file and verifies it against its
sha256. It does: the regenerated file is **byte-identical** to the frozen
`523a4692...`. That is the strongest available statement that the freeze is
real, and it had never been tested.

**Tests.** `tests/test_frozen.py`, 16 guards that run in seconds: the frozen
file's hash and score, the metric reproducing the organisers' 0.8247, the
mean-per-epitope-versus-pooled distinction E1 turned on, single-class folds
returning NaN rather than raising, the data invariants including that no
private row carries a label, and that `submission.validate` actually rejects
each failure it claims to. Tests needing TCRvdb skip rather than fail, so a
clone without the licensed data still checks what it can.

**`make`** encodes the stage graph with honest costs: `check` is five CPU
minutes, while `fold-main` is 7.4 GPU-hours and `md` is ten. None of the
expensive tiers is needed to rebuild the submission.

**One committed cache, after checking what it contains.**
`data/interim/boltz_features.csv` is 1.8 MB of `tcr_id` plus 34 numeric
descriptors -- no sequences, no labels -- and costs 7.4 GPU-hours to
regenerate, so it is committed to remove that barrier. The 29 MB of SCEPTR and
tcrdist caches are not, because they rebuild in under ten CPU minutes. Writing
this up caught a false claim in my own first draft, which asserted the cache
was already committed when `git ls-files data/interim/` was empty.

**`docs/REPRODUCING.md`** carries the two-environment setup, the stage costs,
and a section of traps that each cost real time here: the gene-name format, the
single ORF SCEPTR rejects, OpenMM renaming chains on write, density reference
sets, the two incompatible clustering bases behind E40, the five unsaved
heredocs, and `pgrep -f` matching its own command line.

**Also closed today, and it should not have been left open.** The last entry
asked the user to verify the selected submission in the Kaggle UI. That was
answerable from the API: only two submissions exist, `density3_struct03.csv`
at public 0.849 and the earlier density-only file at 0.839. The higher public
score is auto-selected absent a manual override, so the risk I flagged does
not exist.


## Amendment to the freeze — the deadline is the stopping rule, not my judgement

The entry above froze the predictor on 2026-10-25 and treated that as the end
of the project. That was my decision presented as an evidence-driven one, and
it was the wrong frame. The brief is to keep trying to do better until
2026-11-02 15:00 UTC. Freezing is only correct when no hypothesis is left worth
testing, and four were supplied immediately afterwards.

What the freeze entry got right stands: the submitted file is unchanged, it
rebuilds byte-identically, and nothing has cleared the inclusion threshold. What
it got wrong is the inference from "nothing has cleared the bar yet" to "stop
looking". The file is the **current entry**, not a closed result, and remains so
until the deadline.

---

## E43 — Four new hypotheses from the PI (H1, H2, H4 run; H3 closed on its own rule)

`scripts/40_new_hypotheses.py`, under the unchanged discipline: E3's viability
gate first, zero-parameter rank-average fusion, strict nested LOEO, +0.02 to
enter. Submitted blend baseline 0.8498.

**H3 closed in ten minutes, by the stop condition it came with.** It asked
whether convergent recombination can be distinguished from repeat curation of
one clone. VDJdb carries **no nucleotide column at all**, and on the epitopes
that decide the ranking the donor information is absent:

| epitope | split | subject.id non-null | distinct subjects | clone.id non-null |
|---|---|---:|---:|---:|
| CINGVCWTV | PRIVATE | 78.3% | 10 | 3.5% |
| LLWNGPMAV | PRIVATE | **0.0%** | **0** | 41.3% |
| LTDEMIAQY | PRIVATE | **0.0%** | **0** | **0.0%** |

Its own rule forbids substituting record count for independent rearrangement,
and E27 already showed record count is a label proxy. Closed without spending
compute on it.

**H1 -- biophysical substitution distance. Best of the matrices, still no gain.**

| CDR3 distance | mean AUC |
|---|---:|
| **Grantham (composition, polarity, volume)** | **0.8022** |
| identity (D04's reference) | 0.7989 |
| BLOSUM62 (D04's reference) | 0.7984 |

Grantham is the best of the three, which sharpens D04 rather than contradicting
it: the physicochemical axis is worth a little more than the evolutionary one,
and the whole spread is 0.004. Within-epitope rho against the blend is +0.773,
more decorrelated than SCEPTR against tcrdist3 at +0.909 -- but nested-LOEO
blend gain is **-0.0022**. Rejected.

**H2 -- alpha-chain partner diversity. Beaten by its own trivial confound.**

| feature | mean AUC |
|---|---:|
| alpha-led partner diversity, r=0.2 | 0.7038 |
| **alpha neighbour count alone, r=0.2** | **0.7493** |
| alpha-led partner diversity, r=0.1 | 0.6465 |
| alpha neighbour count alone, r=0.1 | 0.6952 |
| beta-led partner diversity, r=0.2 (control) | 0.6637 |
| beta-led partner diversity, r=0.1 (control) | 0.5981 |

Diversity scores **below** simply counting the neighbours, at both radii, so it
adds nothing over density. The mechanism is visible in the cohort structure: the
median TCR has **zero** alpha-chain neighbours at r=0.2, and only 36.3% have the
two neighbours the diversity measure needs to be defined at all. The
alpha-metaclone structure the hypothesis requires is largely absent here.

One thing it does confirm, through a measure entirely unlike E6's: alpha-led
beats beta-led at both radii (0.6465/0.7038 against 0.5981/0.6637). The
alpha-over-beta asymmetry replicates on a bipartite-graph statistic, not just on
distance. Nested gain -0.0003. Rejected as a feature, kept as a replication.

**H4 -- expression and assembly liabilities. The liabilities are not there.**

| indicator | prevalence | AUC |
|---|---:|---:|
| non-canonical CDR3 termini | **0.0%** | undefined |
| N-glycosylation sequon in CDR3 | 0.4% | 0.5013 |
| extra cysteine beyond the anchor | 1.3% | 0.5106 |
| length deviation from the modal CDR3 | 93.2% | 0.5128 |
| sequon anywhere in the V domain | 100% | 0.5121 |
| hydrophobic fraction of CDR3 | -- | 0.4958 |
| combined liability count | 25.0% | 0.5155 |

**Every one of 5,272 chains has a canonical CDR3**: starts with C, ends with F
or W, verified directly (beta is 100% F-ending, alpha has 83 W-endings, also
canonical). Extra cysteines appear in 34 TCRs and CDR3 sequons in about ten.
Two indicators were discarded by the viability gate for being constant inside a
private epitope.

So the hypothesis is not merely unsupported, its premise is absent: **this
dataset contains essentially no gross, sequence-visible assembly liabilities.**
Either curation removed them upstream or they do not occur at this rate. If
expression failure produces negatives here, it is not through anything these
indicators can see. That is a clean first test of a mechanism nobody had
examined, and it is worth more as a closed door than as a feature.

**Verdict.** Best nested gain across all three is **+0.0000** against a +0.02
threshold. The entry is unchanged.

**One thing to state so E40 is not misread.** E40 withdrew E30's *projected
private AUC*, not the sparsity observation underneath it. That two of three
private epitopes sit below every public epitope on clustered fraction is
directly observable and unaffected. So E21's dominance argument for structure
weight 0.3 over 0.2 still holds, and the weight is not revisited.


## E44 — Grantham as a fourth density member: rejected, and completely

E43 tested Grantham blended on top of the finished model. The PI's follow-up
asks a different question: make it the fourth member *inside* the density rank
average, before the 0.7/0.3 structure blend -- the slot E25/E26 tried the
whole-V-domain alignment in. `scripts/41_grantham_fourth_member.py`.

**The mechanistic objection that killed the V-domain version does not transfer.**
E26 rejected that one as germline confounding, confirmed by E42 showing
framework alone reproduces its ranking at Spearman +0.994. Grantham scores
CDR3, the junction, so a gain would be recognition-relevant sequence. If it
cleared the bar there would be no standing argument against it.

**Members, standalone, with their correlation to the strongest member:**

| member | mean AUC | rho vs SCEPTR default |
|---|---:|---:|
| SCEPTR default | 0.8396 | -- |
| SCEPTR small | 0.8354 | 0.907 |
| V-domain align | 0.8354 | 0.831 |
| tcrdist3 | 0.8329 | 0.878 |
| **Grantham CDR3** | **0.8022** | **0.779** |

**Configurations, blended with structure at 0.3:**

| configuration | density | blended | vs submitted |
|---|---:|---:|---:|
| 3 members (submitted) | 0.8430 | **0.8498** | -- |
| 4 members + Grantham | 0.8405 | 0.8475 | **-0.0024** |
| 4 members + V-domain (E25/E26) | 0.8460 | 0.8524 | +0.0026 |
| 5 members | 0.8449 | 0.8505 | +0.0007 |
| Grantham replaces tcrdist3 | 0.8387 | 0.8469 | -0.0029 |

Grantham is the most decorrelated member available and the weakest by 0.03. An
unweighted rank average is dragged by its weakest term, and decorrelation does
not rescue it -- **E28's lesson restated**: orthogonality is necessary and not
sufficient, because the achievable contribution scales with the margin over
chance. It is not a substitute for tcrdist3 either.

**The oracle sweep closes it completely.** Unweighted averaging is unfair to a
weak member, so the weight on Grantham inside the density term was swept with
hindsight over all five epitopes -- not adoptable under E10's zero-parameter
rule, but decisive as a bound:

| weight on Grantham | 0 | 0.05 | 0.10 | 0.20 | 0.30 | 0.50 |
|---|---:|---:|---:|---:|---:|---:|
| blended AUC | **0.8498** | 0.8494 | 0.8492 | 0.8483 | 0.8467 | 0.8421 |

**The optimum is w=0 and the decline is monotone.** There is no weight at which
Grantham helps, even tuned on the evaluation set. It carries nothing the three
existing members do not already have.

**Viability was not the problem.** Grantham takes a distinct value for every
private-epitope TCR (243/243, 229/229, 130/130), so it passes E3's gate
comfortably. The rejection is on performance alone.

**A cross-check that fell out.** Nested LOEO over all five configurations
selects the V-domain 4-member set in **5 of 5 folds** at 0.8524, independently
reproducing E25's 0.8525 and its 5-of-5 selection on a separately written
harness. That configuration remains unadopted for the reasons E26 and E42 give,
but the agreement is a useful check that this harness and the earlier one say
the same thing.

**Verdict.** -0.0024 against a +0.02 threshold, with an oracle ceiling of
+0.0000. The entry is unchanged.


## E45 — H5: antigen-conditioned alpha-beta pairing compatibility. Closed.

`scripts/42_ab_pairing_compatibility.py`. The reviewer's design, run as three
gates so an early failure closes the branch cheaply. H2 asked how many diverse
beta partners an alpha has; this asks whether the specific alpha-beta
combination a TCR actually carries is the kind its epitope's TCRs make.

**Gate A -- cross-chain dependence exists.** Scoring real pairs against pairs
shuffled inside the epitope, leave-one-out PMI over V/J combinations:

| epitope | n | real | free shuffle | stratified shuffle | strata mobility |
|---|---:|---:|---:|---:|---:|
| GILGFVFTL | 670 | 2.127 | 1.013 | 2.097 | 59.6% |
| GLCTLVAML | 192 | 5.076 | 2.458 | 5.017 | 55.7% |
| YLQPRTFLL | 419 | 2.644 | 1.774 | 2.622 | 32.0% |
| NLVPMVATV | 366 | 3.783 | 2.599 | 3.776 | 21.9% |
| TTDPSFLGRY | 387 | 2.395 | 1.971 | 2.377 | 22.0% |
| CINGVCWTV | 243 | 3.410 | 2.808 | 3.411 | 23.9% |
| LLWNGPMAV | 229 | 2.911 | 2.397 | 2.875 | 34.5% |
| LTDEMIAQY | 130 | 5.331 | 4.557 | 5.331 | 12.3% |

Mean excess over a free shuffle **+1.0126**, so the dependence the 2026 papers
report is clearly present here.

**Gate A' -- and it is composition.** Mean excess over a shuffle restricted to
V/J and length strata is **+0.0213**, about 2% of the free-shuffle excess. The
"mobility" column is the control's own strength and is reported because it
varies: it is the fraction of TCRs a stratified draw actually moves. On
LTDEMIAQY only 12.3% move, which is why its stratified score equals its real
score to three decimals, so that row carries no information. The claim rests on
GILGFVFTL and GLCTLVAML, where a majority move and the residual is still 0.030
and 0.059 against free-shuffle excesses of 1.114 and 2.618.

**Gate B -- and it does not separate labels.** PMI scores mean AUC **0.4400**,
below chance, with the sign flipping across epitopes (GILGFVFTL 0.358,
NLVPMVATV 0.648) -- the same failure shape as `vdjdb_score` in E4 and the
specificity model in E13.

The continuous score, which asks how close a TCR sits to its cohort's
cross-chain regression line, does better: **0.6336**, and above chance on all
five epitopes, which is rare here. But:

| | mean AUC | GIL | GLC | NLV | TTD | YLQ |
|---|---:|---:|---:|---:|---:|---:|
| pairing compatibility | 0.6336 | 0.784 | 0.585 | 0.560 | 0.637 | 0.601 |
| **plain typicality, no pairing** | **0.6378** | 0.779 | 0.641 | 0.518 | 0.667 | 0.584 |

A score built only from how far a TCR's *own* chain features sit from the
cohort mean -- involving no pairing whatever -- scores **slightly higher**, and
the two correlate at Spearman +0.715. So the part of this statistic that
predicts validation is not measuring pairing.

**That is typicality re-encoded for the fourth time**, after density itself,
structure confidence (E14) and the whole-V-domain alignment (E26, E42). Both
scores sit at only +0.26 and +0.27 Spearman against the submitted blend --
genuinely decorrelated -- and still add nothing, which is E28's lesson again:
decorrelation without margin over chance buys nothing.

**Gate C.** Both pass E3's viability gate on all three private epitopes. Nested
LOEO blend gain **-0.0007** and **-0.0013** against a +0.02 threshold. The entry
is unchanged.

**What the reviewer predicted, and what actually happened.** The predicted
failure mode was that every TCR here is a historically observed pair, so a
non-validating one may have ordinary compatibility and the score would not
separate labels. The PMI score failed exactly that way. The continuous score
failed differently and more informatively: it did separate labels, but not
because of pairing.

---

## Two corrections from the same review

**H3's closure was right but its search was not.** E43 closed H3 after grepping
VDJdb's columns for "nucleotide" and finding none. The reviewer pointed out
that VDJdb names the field `cdr3nt`, which that grep would miss, and warned
that where it exists it can be back-inferred from the amino acid sequence by a
rearrangement model rather than observed. Re-checked against the complete column
lists of all three files in this release -- `vdjdb.txt` (22 columns),
`vdjdb.slim.txt` (17) and `vdjdb_full.txt` (35) -- **no `cdr3nt` column exists
in any of them**. The closure stands and now rests on a full enumeration rather
than a keyword match. The caveat is recorded for any future release: an
inferred nucleotide sequence cannot evidence independent rearrangement.

**E43's `liability_sum` was reported in a way that cannot be audited.** Its
prevalence is 25.0% while the three indicators shown beside it are 0%, 0.4% and
1.3%, whose union cannot exceed 1.7%. The reviewer caught the arithmetic. The
cause is that `liability_sum` has a **fourth** term, `len_deviation > 4`, at
23.9% -- and `len_deviation`'s own prevalence was listed at the *different*
threshold `> 0`, which is 93.2%. Two thresholds under one name.

| term in `liability_sum` | prevalence |
|---|---:|
| extra cysteine beyond the anchor | 1.29% (34) |
| N-glycosylation sequon in CDR3 | 0.38% (10) |
| non-canonical CDR3 termini | 0.00% (0) |
| **length deviation > 4 from modal** | **23.86% (629)** |
| union | 25.04% |

E43's conclusion that none of these ranks is unaffected -- the combined score
was 0.5155 -- but the composition should have been stated, and the two
thresholds should not have shared a label.


## E46 — Code review of E40-E45, and a correction that strengthens E41's negative

A review pass over the six scripts written since the review response. Four
suspicious constructs were checked by measurement rather than by reading. Two
were fine, one is immaterial, one was a real defect that invalidated a reported
p-value and, once fixed, makes E41's conclusion stronger.

**Defect 1, material: E41's permutation null was not its observed statistic.**
`family_permutation` shuffled the group label between families and then applied
it to every site in each family, so the null was computed on *family-collapsed*
labels while the observed effect was computed on *site-level* labels. **Fifteen
of 97 families carry both groups, covering 34.1% of sites**, and collapsing
them alone moves the contact effect from **+3.629 to -0.149**. The null and the
observed were therefore different statistics and the p-values of 0.043 and
0.236 were not interpretable. Both are withdrawn.

The cluster bootstrap was unaffected -- it resamples whole families while each
site keeps its own label -- so the confidence intervals stand.

**The corrected family-level analysis, and it goes further than E41 did.**
Replacing the permutation with two tests that have no collapsing step: each
family contributing one mean per group, and a paired comparison inside the 15
families that carry both.

| quantity | site-weighted effect | bootstrap 95% CI | family as unit | paired within mixed families |
|---|---:|---|---|---|
| peptide contacts | +3.629 | [+0.426, +7.832] | 9.043 vs 6.192, **p=0.464** | mean +0.087, **median 0.000**, p=0.674, 5/15 positive |
| site faces peptide | +0.109 | [-0.030, +0.270] | 0.457 vs 0.344, **p=0.249** | mean +0.123, **median 0.000**, p=0.063, 6/15 positive |

**Weighting families equally, neither quantity separates the groups**, and the
paired comparison -- which cannot be confounded by which families fall in which
arm -- has a median difference of exactly zero for both, with fewer than half
the families pointing the right way.

The site-weighted bootstrap CI that excluded zero is explained by the weighting:
three families contribute 129, 92 and 73 sites on their own, so the
site-weighted mean is largely theirs. E41 concluded that the contact count
survives family-level resampling and the binary does not. **The corrected
reading is that neither survives once families are weighted equally.** E31's
effect is weaker than E41 said, which is weaker again than E31 said.

**Defect 2, immaterial: a ddof mismatch in E45.** `continuous_score` computed a
regression slope as `np.cov(xs, ys)[0, 1] / xs.var()`, where `np.cov` defaults
to ddof=1 and `.var()` to ddof=0, inflating the slope by n/(n-1) -- 0.25% at
n=400. Fixed. It cannot affect E45's conclusion, which rests on plain typicality
outscoring the pairing statistic (0.6378 against 0.6336), not on the slope.

**Checked and clean.** `scripts/39` locates CDR1/2/2.5 by substring and removes
the first occurrence, which would be wrong if a loop appeared twice in its own
domain: exactly **1** chain of 5,272 has that, and **0** chains fall back to the
one-residue placeholder, so the V-region ablation is unaffected. `scripts/41`'s
nested selection and rank-average combination were re-read and are correct.

**What this run of self-review says about the method.** Three of the last four
material errors in this project -- E39's union of site residues, E39's
replica-spread guard, and now E41's permutation null -- share one shape: **a
null or a baseline computed on a different data construction from the thing it
is compared against.** That is worth naming, because it is not caught by any of
the checks already in place; the data audit verifies files, the tests verify the
frozen artifact, and neither looks at whether two numbers in the same sentence
were produced the same way.


## E47 — The private shift is published. E40's retreat to 0.846 was too far.

`scripts/43_published_baseline_anchor.py`. Found while reading material the PI
supplied; the finding itself comes from the Immunity paper this competition is
built on, which we have cited since E1.

**E30 projected the private macro at 0.799 and E40 withdrew it**, concluding
that with five epitopes the shift is not measurable. That conclusion was wrong
in an avoidable way: the paper reports tcrdist3's per-epitope AUROC for **all
eight epitopes, private ones included**, in its Figure 4 results text:

> "tcrdist3 showed high performance for the GILGFVFTL (GIL) and YLQ epitopes
> (AUROC 0.87-0.93), good performance for the CIN, NLV, LLWNGPMAV (LLW), and
> TTDPSFLGRY (TTD) epitopes (AUROC 0.80-0.83), and modest performance for the
> GLC and LTDEMIAQY (LTD) epitopes (AUROC 0.66-0.68)."

**The anchor validates against our own measurement**, which is what makes it
usable -- the bands have to describe the estimator we actually hold:

| epitope | published | our tcrdist3 | |
|---|---|---:|---|
| GILGFVFTL | 0.87-0.93 | 0.871 | inside |
| YLQPRTFLL | 0.87-0.93 | 0.928 | inside |
| NLVPMVATV | 0.80-0.83 | 0.809 | inside |
| TTDPSFLGRY | 0.80-0.83 | 0.826 | inside |
| GLCTLVAML | 0.66-0.68 | 0.690 | outside by 0.010 |

Four of five inside, the fifth over by 0.010 -- consistent, given the paper
quotes grouped ranges.

**The shift, measured rather than projected.**

| | |
|---|---:|
| tcrdist3 public macro, our measurement | 0.8247 |
| tcrdist3 private macro, published | **0.753-0.780** (mid 0.767) |
| **public-to-private shift** | **-0.071 to -0.045** (mid -0.058) |
| our public | 0.8498, margin +0.0251 over tcrdist3 |
| **our private estimate** | **0.778-0.805** (mid **0.792**) |

**So E30 was closer than E40 left it.** E30's 0.799 was right in magnitude and
wrong in method; E40 correctly demolished the method and then retreated to
0.846, which the published baseline says is too high by about 0.05. The honest
position is **0.79**, from a measurement of the baseline plus one assumption,
rather than 0.846 from an estimator that cannot see the shift at all.

**The one remaining assumption** is that our +0.0251 margin over tcrdist3
carries to the private epitopes. Untestable here, and it cuts both ways: the
margin comes mostly from density, which E29 and E30 show leans on family
structure the private cohorts have less of, and LTD is where tcrdist3 is
weakest while our model reads the same neighbourhood signal. Against that, the
paper puts AlphaFold3 at **0.84-0.88 on LLW** against tcrdist3's 0.80-0.83, and
our blend carries structure at weight 0.3 -- so on one of the three private
epitopes the structure term is the stronger of the two. That is independent
published support for E21's choice of 0.3 over 0.2, which until now rested only
on our own thinning curves.

**The methodological lesson is the uncomfortable one.** Four experiments
(E21, E30, E32-D10, E40) and a large amount of compute went into estimating a
quantity whose baseline was stated in the paper we had been citing since E1.
Nobody read the results section for it. **Before modelling an unknown, check
whether the literature has already measured it** -- especially the literature
the dataset comes from.


## E48 — Rebuilding the rival configuration: the public lead is pool, and pool does not transfer

`scripts/44_rival_configuration.py`. A parallel effort reports public macro
0.8552 against our 0.8498, from pool-expanded top-8 density rank-averaged with
whole-V-domain alignment top-8. Both components are ones we measured and
rejected, so the question is not whether they lift the public score -- E27 and
E25 say they do -- but whether the lead survives transfer.

Rebuilt entirely from our own implementations: the pool machinery from
`scripts/33`, the V-domain distances cached by `scripts/22`, our metric and our
downsampling protocol.

**The harnesses agree, which is what makes the comparison mean anything.**

| configuration | public |
|---|---:|
| m30 pool-expanded top-8 | 0.8526 |
| m50 V-domain top-8 | 0.8319 |
| **m30 + m50** | **0.8549** (reported 0.8552, gap -0.0003) |
| top-8 with no pool | 0.8414 |
| our frozen entry | 0.8498 |
| tcrdist3 baseline | 0.8247 |

Their public advantage is **+0.0051**.

**After the private-realisable downsampling E27 applied to our own pool gain** --
each public pool cut to the private epitopes' own distinct-clonotype
multipliers, 20 repeats:

| | |
|---|---:|
| their configuration, full public pool | 0.8549 |
| **their configuration, private pool density** | **0.8476** |
| pool advantage that cannot transfer | **+0.0073** |
| our frozen entry, which uses no pool at all | 0.8498 |
| **advantage that survives** | **-0.0022** |

**Their entire public lead, and more, is pool.** At the pool density the private
epitopes can actually supply, the two configurations swap places.

**Under E47's published anchor:**

| configuration | margin carried | private |
|---|---:|---:|
| our frozen entry | +0.0251 | **0.792** |
| theirs, public margin | +0.0302 | 0.797 |
| theirs, realisable margin | +0.0229 | **0.790** |

0.790 against 0.792 is far inside the +/-0.023 standard error of a
three-epitope estimate. **The two configurations are indistinguishable on the
private split**, and the public gap that looks decisive is an artefact of a
resource the private epitopes do not have.

Their own private estimate is 0.775, more conservative than the anchor gives,
the difference being that they take the baseline at 0.76 where E47's published
midpoint is 0.767. The direction agrees: every team's private score will sit
well below its public one.

**One thing this leaves open rather than settles.** Their second component, the
V-domain alignment, carries **no pool discount** -- it changes how two TCRs are
compared, not how many there are to compare against, so it transfers in full.
E44 measured it as a fourth density member at +0.0026 with nested LOEO selecting
it in 5 of 5 folds. That is below the +0.02 threshold and E26 and E42 showed its
ranking comes from framework, a region that contacts nothing. But "below
threshold" and "transfers without discount" are different properties, and the
frozen entry was not changed on the strength of this.


## Standing gates — 2026-10-26

Two checks promoted out of individual scripts into `src/immrep/validation.py`,
so a new feature meets them before its LOEO score is believed rather than
after. Both have tests.

**`private_variance_gate`** is E3's rule, which has been applied by hand since
E3 and caught two H4 indicators in E43: a feature constant inside a private
epitope contributes exactly zero where the ranking is decided, however well it
scores publicly.

**`typicality_gate`** is new, and exists because typicality is this project's
gravity well. Four features that looked like different ideas all turned out to
re-encode *how ordinary a TCR's own chain features are within its cohort*:

| feature | what it looked like | what it was |
|---|---|---|
| neighbourhood density | convergent antigen-driven clustering | partly typicality |
| Boltz interface confidence (E14) | structural plausibility of the complex | how canonical the TCR is |
| whole-V-domain alignment (E26, E42) | a decorrelated sequence distance | framework, which contacts nothing, at rho +0.994 |
| alpha-beta pairing compatibility (E45) | cross-chain selection by antigen | typicality, which outscored it |

E45 is the sharpest: a crude score -- the summed negative z-score of CDR3
length, charge and mean residue mass on both chains, within epitope, **with no
pairing in it at all** -- reaches **0.6378** against the pairing statistic's
0.6336, correlating at rho +0.715.

The gate reports the feature's AUC, typicality's, their within-epitope rho and
a verdict. Validated on the two known cases: the submitted blend beats it
(0.8498 against 0.6378, rho +0.275) and H5's pairing score does not (0.6336,
rho +0.715). Three tests pin that behaviour, including the two-case
discrimination, so the gate cannot quietly stop working.

**Memory corrected at the same time.** It still carried E40's retreat to 0.846,
which E47 superseded. Now records 0.79 anchored on the published bands, with
the history of how E30 and E40 each got half of it right, and a new entry on
the typicality well.


---

## E49 — Is typicality a motif detector? Response diversity, conditioned and oracled

**Why.** E45 left typicality scoring 0.6378 with a very uneven per-epitope
profile (GIL 0.779 down to NLV 0.518). The proposed mechanism was that "typical
shape" proxies "belongs to the epitope's dominant motif", which would make the
weight typicality deserves a function of how concentrated the response is — and
that concentration is observable without labels, hence observable on the private
epitopes. `scripts/46_response_diversity.py`.

Two things pre-registered before the run, because both could otherwise have been
chosen after seeing results: the primary diversity index is
`largest_family_share`, and the primary stratification is TRBV. All four indices
are rarefied to the smallest cohort (130 against 670) since every one of them is
biased by n.

**Part 1 — label-free diversity, all eight epitopes.** Rarefied, 200 draws.

| epitope | n | split | largest_family_share | family_simpson | trbv_simpson | len_simpson |
|---|---:|---|---:|---:|---:|---:|
| GILGFVFTL | 670 | public | 0.1064 | 0.9707 | **0.5699** | 0.5136 |
| GLCTLVAML | 192 | public | **0.1654** | 0.9537 | 0.8945 | 0.7437 |
| NLVPMVATV | 366 | public | 0.0464 | 0.9896 | 0.9558 | 0.8413 |
| TTDPSFLGRY | 387 | public | 0.0177 | 0.9920 | 0.8841 | 0.8031 |
| YLQPRTFLL | 419 | public | 0.0364 | 0.9899 | 0.9238 | 0.7830 |
| CINGVCWTV | 243 | private | 0.0124 | 0.9922 | 0.9441 | 0.8405 |
| LLWNGPMAV | 229 | private | 0.0185 | 0.9920 | 0.9460 | 0.7489 |
| LTDEMIAQY | 130 | private | 0.0538 | 0.9891 | 0.9464 | 0.8232 |

**Part 2 — the between-epitope association.** Exact permutation p throughout:
`scipy.stats.spearmanr` uses a t-approximation that reports p=0.000 for rho=-1
at n=5, where the exact value is 2/120 = 0.0167.

| index | rho vs typicality AUC | exact p | family-wise (Sidak, 4 indices) |
|---|---:|---:|---:|
| largest_family_share (pre-registered) | +0.100 | 0.950 | 1.000 |
| family_simpson | -0.100 | 0.950 | 1.000 |
| **trbv_simpson** | **-1.000** | **0.017** | **0.065** |
| len_simpson | -0.700 | 0.233 | 0.655 |

**The pre-registered index fails and a secondary one orders perfectly.** That is
the multiplicity trap the pre-registration existed to expose, so it is reported
as such. Three things make it more than a lucky draw among four:

- GILGFVFTL is an outlier on `trbv_simpson` (0.57 against 0.88–0.96), but
  dropping it leaves rho = **-1.000** on the remaining four (exact p 0.083).
- rho is -1.000 across all six rarefaction seeds tried, and -1.000 (p 0.017)
  with no rarefaction at all.
- the two indices that fail are both *family*-level (near-duplicate clusters at
  cosine 0.12) while the one that works is *germline*-level. That is the same
  split E42 found when whole-V-domain alignment tracked typicality at rho
  +0.994 through framework that contacts nothing: typicality measures the
  germline half, not the CDR3 half.

All three private cohorts sit at 0.946–0.949, NLV's diverse end, so the
association predicts typicality near chance on every private epitope.

**Part 3 — the within-epitope conditional test, and why it is void.** Typicality's
pair-weighted AUC counting only positive-negative pairs inside the same TRBV
stratum, with a cluster bootstrap over families:

| epitope | AUC all pairs | pairs | AUC within TRBV | pairs | retained | drop | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---|
| GILGFVFTL | 0.779 | 73914 | 0.664 | 19778 | 26.8% | **+0.115** | [0.056, 0.183] |
| GLCTLVAML | 0.641 | 6160 | 0.676 | 519 | 8.4% | -0.035 | [-0.161, 0.157] |
| NLVPMVATV | 0.518 | 20566 | 0.620 | 643 | 3.1% | -0.102 | [-0.200, 0.045] |
| TTDPSFLGRY | 0.667 | 29240 | 0.684 | 2642 | 9.0% | -0.016 | [-0.083, 0.059] |
| YLQPRTFLL | 0.584 | 35392 | 0.585 | 2240 | 6.3% | -0.002 | [-0.070, 0.072] |

Only GILGFVFTL's CI excludes zero, and only GILGFVFTL retains enough pairs to
have power. Family-level conditioning is unusable at any cohort size here: 257,
68, 9, 0 and 8 within-family pairs across the five.

**Part 3b — the control, which voids it.** Added after a first run that lacked
one. The same TRBV conditioning applied to other scores, as the drop in AUC:

| epitope | density alone | submitted blend | typicality |
|---|---:|---:|---:|
| GILGFVFTL | **+0.119** | **+0.141** | +0.115 |
| GLCTLVAML | -0.141 | -0.130 | -0.035 |
| NLVPMVATV | -0.024 | -0.017 | -0.102 |
| TTDPSFLGRY | -0.007 | +0.024 | -0.016 |
| YLQPRTFLL | +0.030 | +0.039 | -0.002 |

On GILGFVFTL, conditioning on TRBV costs density alone **more** than it costs
typicality. The one drop with a CI excluding zero is therefore a property of
TRBV stratification on that cohort, not evidence that typicality acts through
motif membership. **The mediation test fails its own control and establishes
nothing.** Without Part 3b the GIL row would have read as confirmation.

**Part 4 — and it is worth nothing on the leaderboard anyway.** Oracle first,
since an oracle that cannot clear the bar closes the direction whatever the
mechanism turns out to be.

| | score | vs frozen |
|---|---:|---:|
| frozen entry | 0.8498 | — |
| oracle *per-epitope* typicality weight | 0.8523 | **+0.0025** |
| oracle single global weight (w=0.10) | 0.8503 | +0.0004 |

Per-epitope oracle weights: GIL 0.20, GLC 0.30, TTD 0.10, NLV 0.05, YLQ 0.00 —
ordered exactly as `trbv_simpson` predicts, which is the association again and
not independent of it. Conditioning buys an oracle **+0.0020** over a single
global weight, against a +0.02 inclusion bar and E15's -0.0048 selection cost.
No label-free estimate of a per-epitope weight can beat its own oracle, so no
nested rule was fitted. Typicality passes both standing gates
(`private_variance_gate` and E3 viability, 243/229/130 distinct values) and
still has nowhere to go.

**Settled.** The mechanism hypothesis is **not supported**: the only test with
power fails its control. What survives is an association — germline TRBV
concentration rank-orders typicality's per-epitope AUC perfectly across five
epitopes, robustly to the outlier and to the rarefaction, exact p 0.017 and
family-wise 0.065 — which five epitopes cannot promote to a mechanism, and which
an oracle prices at +0.0025. The epitope-conditioned weighting direction named
in memory as this project's gap is closed for typicality specifically. The entry
is unchanged.

**Two method notes worth keeping.** A perfect rank correlation at n=5 needs an
exact p; the asymptotic one is wrong by a factor of 60 here. And a mediation
test needs the same conditioning applied to a score that is not the hypothesis,
or a stratification that is simply hard will read as confirmation.

---

## E50 — Pair-conditioned weights: the third pairwise formulation

**Why.** Three things get called pairwise learning and two were already settled.
E16 put a pairwise ranking surrogate under gradient descent and it failed
structurally — the surrogate is *negatively* correlated with AUC across the
sigma grid (−0.39, −0.51) and only aligns as its temperature goes to zero, where
it becomes the 1−AUC step function and the gradients vanish. E10 and E13 fitted
pointwise learners and lost to the zero-parameter rank average. Never built: a
model whose training example is a *pair*.

The motivation was E22. Error rate varies enormously by a pair's cluster state,
and most of the loss sits in one class, so a pair-conditioned weight could
express "when neither member has relatives, trust density less" — which no
pointwise score can. `scripts/47_pairwise_learning.py`.

**Only four of E22's five classes are usable.** Separating `pos_clu_neg_sing`
from `pos_clu_neg_sing` requires knowing which member is the positive, i.e. the
label. Their union — exactly one member clustered — is label-free, as are the
other three. Conditioning on the five would be leakage, and the rule has to work
on the private split where nothing is labelled.

**The model.** `decision(i,j) = Σ_k w[k, class(i,j)] · (f_k(i) − f_k(j))` over
the frozen entry's own two components plus typicality. Every term is an
antisymmetric difference scaled by a symmetric pair property, so
`P(i beats j) = 1 − P(j beats i)` holds by construction, and `w_struct = 0.3,
w_typ = 0` reproduces the frozen entry — **the anchor, which matched to
0.000001**, so everything below describes the model actually entered.

Because class partitions the pairs, a class's weights move only its own pairs,
so the oracle separates into four independent searches and the ceiling of the
whole idea is cheap to price.

**Where the loss sits, under the label-free partition.**

| class | share of pairs | within-class AUC | share of all loss |
|---|---:|---:|---:|
| both_singleton | 49.2% | 0.818 | **59.8%** |
| one_clustered | 43.6% | 0.903 | 28.1% |
| diff_cluster | 6.9% | 0.745 | 11.7% |
| same_cluster | 0.3% | 0.734 | 0.5% |

E22's 62% for `both_singleton` survives the merge of the two adversarial classes
into one observable one, at 59.8%.

**The oracle, one weight pair per class.**

| class | density | structure | typicality | gain |
|---|---:|---:|---:|---:|
| both_singleton | 0.80 | 0.20 | 0.00 | +0.0010 |
| one_clustered | 0.65 | 0.20 | 0.15 | +0.0016 |
| same_cluster | 0.85 | 0.05 | 0.10 | +0.0002 |
| diff_cluster | **0.25** | **0.40** | **0.35** | +0.0012 |

| | score | vs frozen |
|---|---:|---:|
| frozen entry | 0.8498 | — |
| oracle class-conditioned | 0.8537 | **+0.0039** |
| oracle single global weight (struct 0.20, typ 0.10) | 0.8515 | +0.0017 |

Pair-conditioning buys an oracle **+0.0022** over one global weight. The weight
heterogeneity is real — `diff_cluster` wants density down to 0.25 and
typicality up to 0.35, nothing like the frozen weights — but that class holds
6.9% of pairs, so the whole effect is worth +0.0012.

**Settled, and on a sharper criterion than the bar.** The oracle ceiling
+0.0039 is **smaller than E15's measured selection cost of −0.0048** over a
search of comparable size (4 classes × 231 weight pairs here against 323
configurations there). A fitted version is therefore expected to come out net
*negative*, not merely short of +0.02. Nothing was fitted.

And the motivation is answered directly: `both_singleton` holds 59.8% of all
loss and the oracle moves it by **+0.0010**. **Reweighting cannot reach the place
the loss lives** — when neither member of a pair has relatives, there is no
combination of these three components that orders them better. The finer
conditioning confirms it: mutual distance is the one genuinely pairwise
continuous quantity with no pointwise counterpart, and within `both_singleton`
its terciles want structure 0.20/0.30/0.30 and typicality 0.00 throughout, i.e.
essentially the frozen weights everywhere.

**All three pairwise formulations are now closed, each with a mechanism**: the
ranking objective by surrogate–metric mismatch on a heavy-tailed score
distribution (E16), the pointwise learners by five folds (E10, E13), and
pair-conditioned weights by the loss sitting in the one class that reweighting
cannot improve (here). The entry is unchanged.

---

## E51 — A public-maximising configuration, and where its gain actually sits

**Why.** The rules permit several teams and several submissions, which changes
what is worth building. The frozen entry maximises *expected private* score and
deliberately carries no pool expansion, because E27 priced the public gain
(+0.0112) as almost entirely unavailable to the private epitopes (+0.0013
realisable). A second entry is how to take the other side of that measurement.
`scripts/48_public_maximum.py`.

The opportunity was visible in E48 and nobody had measured it: **HeidAI's
configuration is pool-expanded density + whole-V-domain alignment; the frozen
entry is three-space density + Boltz structure. Neither holds the other's term.**

**Public macro, equal-weight rank averages throughout.**

| configuration | public | realisable | pool-dependent | vs frozen (public / realisable) |
|---|---:|---:|---:|---|
| **D  density3 + pooled m30 + m50 + structure** | **0.8601** | **0.8561** | 0.0039 | **+0.0102 / +0.0063** |
| C  HeidAI rebuild + structure 0.3 | 0.8588 | 0.8534 | 0.0055 | +0.0090 / +0.0035 |
| A  frozen, density3 → pooled m30 | 0.8583 | 0.8499 | 0.0084 | +0.0085 / +0.0001 |
| E  (density3 + pooled m30 + m50)/3 + str 0.3 | 0.8582 | 0.8548 | 0.0034 | +0.0083 / +0.0050 |
| B  pooled m30 + m50 + structure | 0.8572 | 0.8521 | 0.0052 | +0.0074 / +0.0023 |
| F  **D with no pool at all** | 0.8557 | — | — | +0.0059 / — |
| HeidAI rebuild (m30 pooled + m50) | 0.8549 | 0.8476 | 0.0073 | +0.0051 / **−0.0022** |
| frozen entry | 0.8498 | 0.8498 | 0.0000 | 0 / 0 |

`realisable` thins every public pool to the private multipliers (1.038, 1.045,
1.834), 20 repeats, exactly as E27 and E48 priced the discount.

**The column that decides transfer is not the public macro.** The private split
is 61–93% pairs where neither member has relatives, and E30 showed the projected
private macro tracks the within-both-singleton AUC almost exactly.

| configuration | both-singleton AUC | vs frozen | worst per-epitope delta |
|---|---:|---:|---:|
| **D** | **0.8084** | **+0.0143** | **+0.0021** |
| C | 0.8083 | +0.0143 | −0.0004 |
| B | 0.8075 | +0.0134 | −0.0056 |
| E | 0.8056 | +0.0115 | +0.0026 |
| **F (no pool)** | 0.8027 | **+0.0087** | +0.0007 |
| A | 0.8005 | +0.0064 | −0.0038 |
| HeidAI rebuild | 0.7984 | +0.0043 | −0.0076 |
| frozen entry | 0.7941 | 0 | 0 |

Three findings, in order of importance.

1. **The gain is concentrated where the private score is decided.** D improves
   both-singleton AUC by +0.0143 against its +0.0102 public macro — a *larger*
   improvement in the harder class than overall, which is the opposite of the
   pool-expansion pattern E27 rejected.
2. **It holds on every epitope with no sign flip** (worst delta +0.0021). D, E
   and F are the only candidates with that property; C, B, A and HeidAI's own
   configuration each make at least one epitope worse. Sign instability across
   epitopes is the exact failure that killed `vdjdb_score` (E4), Pgen (E7) and
   the specificity model (E13), so its absence here matters more than the macro.
3. **Most of it is not the pool.** F — the same four members with the pool
   deleted entirely — still carries +0.0059 public and **+0.0087 of the
   both-singleton gain, 61% of D's**. So the advantage is mainly the V-domain
   term and the top-8 aggregation, neither of which needs a pool, and the
   transfer argument that rejected pool expansion does not apply to most of it.

On E26/E42's objection that whole-V-domain alignment tracks typicality at rho
+0.994 through framework that contacts nothing: that argued m50 carries no new
*biological* information, and it still stands. What m50 contributes here is a
structurally different *metric* — global alignment rather than embedding cosine
— at a standalone 0.8319, well above chance. That distinguishes it from Grantham
(E44), which had decorrelation but no margin, and is consistent with E28's rule
rather than contrary to it.

**Selection cost, measured.** Nested LOEO over the candidate *choice*: D is
picked in 4 of 5 folds, nested 0.8577 against the fixed candidate's 0.8601, so
selection costs **−0.0024** here rather than E15's −0.0048, because the
candidates are so similar.

**Under E47's published anchor.** D projects to private **0.802** on its public
margin and **0.798** on its realisable margin, against the frozen entry's 0.792.

**Settled, and deliberately not a replacement.** +0.0102 does not clear the
+0.02 bar adopted in E10, so under this project's own standing rule the frozen
entry stays as entered. But the bar was set for a component justified by its
LOEO macro, and this evidence is of a different kind — right class, every
epitope, mostly pool-free. Two additional submission files were written and
validated, with the frozen entry untouched and its byte-identical guard still
passing:

| file | configuration | public | sha256 (first 16) |
|---|---|---:|---|
| `d_density3_m30pool_m50_struct.csv` | D, takes the pool bet | 0.8601 | `64b3dd4076a51fbd` |
| `f_density3_m30_m50_struct_nopool.csv` | F, no pool anywhere | 0.8557 | `3314b7d6319d17b9` |

D passes both standing gates (private variance 236/220/126; typicality
"carries something typicality does not", 0.8601 against 0.6378 at rho +0.265).


**Amended: sensitivity to the pool multiplier (2026-10-27).** A collaborator
reports that adding IEDB targeted at epitopes whose VDJdb multiplier is below 2
raises LTDEMIAQY from 1.04 to **1.72** at no public cost, and separately
estimates M30's private pool gain at about +0.011 against the +0.0013 E27
measured. This repository holds no IEDB download, so the claim cannot be
verified here — but the realisable column can be recomputed under it, since
that multiplier is the single input the pool rejection rests on.

| configuration | realisable, VDJdb-only | realisable, LTD at 1.72 | vs frozen |
|---|---:|---:|---|
| D | 0.8561 | 0.8565 | +0.0063 → +0.0067 |
| C | 0.8534 | 0.8540 | +0.0035 → +0.0042 |
| A (pure pool bet) | 0.8499 | 0.8509 | +0.0001 → +0.0011 |
| HeidAI rebuild | 0.8476 | 0.8484 | −0.0022 → **−0.0015** |

**The conclusion does not depend on the multiplier.** Even under the more
optimistic figure the rival configuration stays negative against the frozen
entry on the private split, and D's advantage moves by +0.0004. LTDEMIAQY is
one of three private epitopes and 1.038 → 1.72 adds a modest number of
neighbours to a top-8 density on a 130-row cohort, so the input has little
leverage. The disagreement with the +0.011 estimate therefore does not need to
be resolved to act: our thinning measurement is direct, theirs appears to apply
a discount factor, and the decision is the same either way.

---

## E52 — Classifying singletons only

**Why, and why it is not E10 again.** The idea is the user's. Every supervised
attempt so far (E10, E13, E50) optimised the *overall* AUC, which the easy class
dominates: `one_clustered` pairs have within-class AUC 0.903 and 43.6% of the
mass (E50), so almost all of a fitted model's pressure goes where density
already works. Restricting the training set to singletons replaces the objective
with the hard class alone, and that class is 59.8% of all loss and 61-93% of the
private split's pairs. `scripts/49_singleton_classifier.py`.

**Two pieces of machinery, both validated before use.**

*The splice.* The metric is a per-epitope AUC over all rows, so a separate
singleton ranking has to be reinserted without disturbing anything else. Within
an epitope the singletons keep the rank slots they already occupy in the base
score and only their ordering among themselves is replaced. Splicing a score
into itself reproduces it to **0.849831 against 0.849831**, so the operation is
lossless and the only thing a result can move is the both-singleton part.

*The features.* Within-epitope percentile ranks computed **among singletons**,
so every epitope's feature distribution is uniform on [0,1] and nothing
epitope-specific can survive into a coefficient. 1,281 labelled singletons:
GIL 317, TTD 328, YLQ 308, NLV 252, GLC 76.

**Part 1 — single features within labelled singletons.**

| feature | mean AUC | worst epitope | neighbour-free |
|---|---:|---:|---|
| the frozen score itself | **0.7941** | 0.6549 | no |
| nn3_dist | 0.7793 | 0.5994 | no |
| nn1_dist | 0.7711 | 0.5953 | no |
| density (kernel) | 0.7702 | 0.6043 | no |
| cons_ms_iptm_mean | **0.7292** | 0.6440 | **yes** |
| cons_ms_tcr_pep_iptm_mean | 0.7276 | 0.6373 | yes |
| complex_ipde | 0.7142 | 0.6259 | yes |

**Part 2 — classifiers, nested LOEO over epitopes.**

| features | model | top-k | singleton AUC | spliced macro |
|---|---|---|---:|---:|
| all | logistic L2 | 5 | **0.7949** | **0.8510** |
| all | logistic L2 (C=1) | 5 | 0.7943 | 0.8492 |
| all | boosting depth 2 | 5 | 0.7863 | 0.8430 |
| neighbour-free | logistic L2 | 15 | **0.7277** | 0.7995 |
| neighbour-free | logistic L2 | 5 | 0.7250 | 0.7990 |
| neighbour-free | boosting depth 2 | 15 | 0.7034 | 0.7869 |
| — | the frozen entry | — | 0.7941 | 0.8498 |

**Settled, with three readings in increasing order of interest.**

1. **The restriction buys +0.0012 on the macro and +0.0008 inside the class it
   was built for.** The best model's feature set contains the frozen score
   itself, and top-5 selection picks it first, so what the classifier learned is
   "use the score you were given". Fitting on four epitopes added nothing to it.

2. **Deleting the neighbour features costs −0.069 within singletons** (0.7277
   against 0.7941), and the fitted neighbour-free model does not even beat the
   best *single* neighbour-free feature (0.7292, multi-sample ipTM consistency).
   Fifteen features on four epitope folds is worse than picking one.

3. **The diagnosis, which is the point: "singleton" is not a neighbour-free
   population.** `nn3_dist` scores 0.7793 and `nn1_dist` 0.7711 *within*
   singletons -- higher than kernel density's 0.7702. The 0.12 radius creates a
   category, but the variable underneath it is continuous and keeps ranking
   below the threshold. **There is no subpopulation of this dataset for which the
   neighbourhood carries no information**, so "a signal that ranks TCRs with no
   relatives" was the wrong thing to look for: the TCRs in question do have
   relatives, just further away.

That closes the loop with E50 from the other side. Reweighting could not improve
`both_singleton` (oracle +0.0010) because the active signal inside that class is
still nearest-neighbour distance, which density already carries. And it explains
E51 mechanistically: candidate D's +0.0143 on both-singleton pairs comes from
its **top-8** members, and top-k aggregation beats a kernel at sigma=0.1
precisely where neighbours are far -- exactly what `nn3_dist > density` says
here. **E51 was already the correct exploitation of this finding, arrived at
from the other direction.**

Gain +0.0012 against the +0.02 bar. The entry is unchanged.

---

## E53 — Why a TCR is a singleton, and whether the radius is principled

Both questions are the user's, and they turn out to have one answer.
`scripts/50_singleton_cause_and_cutoff.py`.

### Part 1 — why a singleton?

Three causes make opposite predictions, and the target — *is this TCR a
singleton* — carries no label, so it can be predicted from features directly and
measured **within each epitope**, which holds cohort size fixed.

| feature | proposed cause | AUC for "is a singleton" | inverted |
|---|---|---:|---:|
| cdr3_len_total | chain atypicality | **0.635** | — |
| ms_iptm_sd | — | 0.370 | 0.630 |
| typicality | chain atypicality | 0.364 | 0.636 |
| complex_plddt | structurally odd | 0.363 | 0.637 |
| consist_tcr_rmsd_mean | — | 0.357 | 0.643 |
| log10_pgen_paired | hard to generate | 0.335 | **0.665** |
| pae_tcr_pep_min | structurally odd | 0.318 | 0.682 |
| ms_iptm_mean | structurally odd | 0.305 | 0.695 |
| iptm | structurally odd | 0.299 | **0.701** |

All nine together, in-sample within epitope: **0.850** — an upper bound, but it
settles the first question: **singleton status is substantially a property of the
TCR, not only of what was collected.**

**The structural reading is directionally right, and it is not a separate cause.**
Singletons have lower interface confidence (ipTM inverted 0.701), worse PAE
(0.682), lower Pgen (0.665), more atypical chain features (0.636) and longer
CDR3s (0.635). But these are facets of one quantity, not independent evidence:
E14 and E36 established that Boltz interface confidence *is* a canonicality
measure, and here `ms_iptm_sd` correlates with `ms_iptm_mean` at Spearman
**+0.843** and `consist_tcr_rmsd_mean` with it at **+0.676** — so the two
features that look like reproducibility measures are attenuated copies of the
confidence mean, not counter-evidence to it. A first reading of the output
treated them as contradicting the confidence features; they do not.

The unified reading: **a long CDR3 is improbable to generate, docks less
canonically, looks atypical in chain features, and has few convergent
relatives.** "Its structure is odd so it cannot make neighbours" and "it is hard
to generate so the relatives do not exist" are the same explanation from two
sides.

**Sampling is refuted at the epitope level.** Singleton rate against cohort
size: rho **−0.214** (n=8), and the extremes run the wrong way — GLCTLVAML has
the *smallest* cohort (192) and the *fewest* singletons (44.8%), while
CINGVCWTV has 243 rows and 97.5%.

| epitope | n | singleton rate | split |
|---|---:|---:|---|
| CINGVCWTV | 243 | **97.5%** | private |
| LLWNGPMAV | 229 | 95.6% | private |
| **TTDPSFLGRY** | 387 | **95.3%** | public |
| LTDEMIAQY | 130 | 85.4% | private |
| NLVPMVATV | 366 | 81.7% | public |
| YLQPRTFLL | 419 | 81.6% | public |
| GILGFVFTL | 670 | 52.4% | public |
| GLCTLVAML | 192 | 44.8% | public |

**Worth keeping: TTDPSFLGRY at 95.3% is the public epitope most like the private
ones.** It is the best single-epitope proxy for private transfer, which is
useful because LOEO weights all five equally and E30 showed LOEO is biased
upward precisely by family richness.

### Part 2 — is the radius principled? No.

**Smooth, no knee.** Validation rate by nn3 decile, taking each epitope's own
rate and averaging epitopes so the fourfold base-rate differences cannot mix:

| decile | nn3 mean | rate (per-epitope mean) |
|---:|---:|---:|
| 0 | 0.078 | 0.906 |
| 2 | 0.220 | 0.787 |
| 4 | 0.351 | 0.517 |
| 6 | 0.474 | 0.231 |
| 8 | 0.576 | 0.197 |
| 9 | 0.634 | 0.117 |

Largest single step **−0.216** against a mean absolute step of **0.107** — a
ratio of 2.0. A real threshold would put one step far above the rest.

**Against a background of other epitopes' TCRs**, within-epitope distance stops
being special at radius 0.45–0.49 for seven epitopes and at **0.21 for
CINGVCWTV**. At the radius in use, CIN has **zero** TCRs with a within-cohort
nn3 ≤ 0.12 while 0.41% have a cross-epitope neighbour that close — enrichment
**0.0**. A private epitope with no convergent structure at all at this radius,
published nonetheless at 0.80–0.83, because the continuous variable still ranks.

**The decisive sweep.** Binary clustered indicator at every radius 0.02–0.40,
against the continuous `nn3_dist` it is a threshold of:

| | mean per-epitope AUC |
|---|---:|
| **continuous nn3_dist** | **0.8334** |
| best binary, radius 0.36 (47.8% clustered) | 0.7544 |
| binary at the radius in use, 0.12 | 0.6725 |

**No radius recovers what the continuous variable carries — the best possible
binary loses 0.0790, and 0.12 loses 0.1609.** "Singleton" is a lossy summary of
`nn3_dist`, not a category.

This does not mean the model leaves 0.16 on the table: 0.12 is used for analysis
bookkeeping — defining families and pair classes — not as a predictor. But it
does mean every class-based analysis here (E22, E50, E52) rests on a lossy
summary, which is a limitation worth recording rather than a result to act on.

**A negative worth stating.** Enrichment is label-free, so unlike E49 this gives
all eight epitopes an x value against E47's published AUROC. It does not predict
difficulty: rho **+0.309, p 0.457**, broken by GLCTLVAML — second-highest
enrichment (2.020) and the lowest published band (0.67). Same outlier structure
as E49. The mechanism candidate is E9's: a tight cluster can be a knot of
co-annotated near-duplicates rather than a real motif family.

### Settled

Singleton status is intrinsic and the structural hypothesis is one face of it;
the radius is arbitrary and strictly lossy. **The practical consequence is
already banked**: candidate D's top-8 members read the continuous variable
directly, which is why E51 measured +0.0143 on both-singleton pairs where
reweighting (E50) and a singleton-only classifier (E52) each got about +0.001.
The entry is unchanged.

---

## E54 — Is family label purity more than the continuous density relation?

The design is the user's. E20 observed 45 families of ≥3 labelled members with
26 pure-positive, 3 pure-negative and the **20–60% purity band completely
empty**, and could not distinguish two readings:

- **A** families look pure because their members all sit at low `nn3_dist`,
  where the validation rate is high anyway, and families are small;
- **B** each family shares something extra — a recognition state, or a study or
  donor — that makes members agree beyond their position on the distance axis.

`scripts/52_family_purity_null.py`. The null permutes the label within
**(epitope × nn3 bin)**, preserving each epitope's label ratio, the
density-to-label relation at the bins' resolution, and every edge and family,
while destroying family-level agreement not already implied by position.

**Two design points that have bitten this project before.** The graph is built
on the **full cohort** — every row of the epitope including the ambiguous ones
that carry no label — and purity is measured on labelled members only; two
clustering bases already exist here (`scripts/20` labelled rows, `scripts/25`
the whole cohort) and mixing them is what made E30's accounting inconsistent.
Building it correctly **reproduces E20's 26 / 3 / 16 exactly**, so that headline
was not affected by the basis. And the bin count is swept rather than fixed.

**The density relation is flat inside this population.** Across nn3 deciles of
the 390 TCRs in families of ≥3, the validation rate runs 0.927, 0.925, 0.895,
0.949, 0.846, 0.949, 0.842, 0.947, 0.974, 0.825 — **no trend**, because family
members all sit at nn3 0.03–0.14 by construction. So conditioning on density
barely binds, and the question reduces cleanly to whether labels clump by family.

**The sweep, 2,000 draws per setting.**

| bins | strata | movable | pure⁺ obs | pure⁺ null | p | concordance obs | null | p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5 | 98.2% | 26 | 26.5 | 0.707 | 0.9092 | 0.8290 | **0.0005** |
| 3 | 15 | 90.5% | 26 | 26.7 | 0.742 | 0.9092 | 0.8272 | **0.0005** |
| 5 | 25 | 81.5% | 26 | 25.9 | 0.577 | 0.9092 | 0.8361 | **0.0005** |
| 10 | 47 | 64.6% | 26 | 24.5 | 0.270 | 0.9092 | 0.8541 | **0.0005** |
| 20 | 87 | 45.1% | 26 | 24.6 | 0.289 | 0.9092 | 0.8609 | **0.0020** |

**A split verdict, and the two halves do not conflict.**

**Result A for E20's headline.** The *number* of pure families is exactly what
the null produces at every density resolution — 26 observed against 24.5–26.7,
p 0.27 to 0.74. **The 26 pure families and the empty 20–60% band need no extra
mechanism**; they are what a high base rate plus small families gives.

**Result B, narrowed to six objects.** Pair-weighted agreement is above the null
everywhere (0.9092 against 0.827–0.861, p ≤ 0.002), and decomposing by family
size shows where all of it lives:

| family size | families | pairs | concordance | null | excess | p |
|---|---:|---:|---:|---:|---:|---:|
| 3–4 | 23 | 87 | 0.862 | 0.810 | +0.052 | 0.162 |
| 5–9 | 12 | 190 | 0.795 | 0.809 | −0.014 | 0.649 |
| 10–19 | 4 | 393 | 0.924 | 0.876 | +0.048 | 0.196 |
| **20+** | **6** | **2996** | **0.916** | 0.857 | **+0.059** | **0.006** |

Six families hold 2,996 of 3,666 pairs (82%) and are the only significant band.
A wider permutation pool — shuffling among *all* labelled rows of the stratum so
family members can trade with TCRs that have no family — gives the same answer
(concordance p 0.0005, pure-positive p 0.481), so this is not an artifact of
restricting the pool.

**The six, named, because a result resting on six objects must be inspectable.**

| family | labelled | validated | distinct references | shared by all | dominant TRBV |
|---|---:|---:|---:|---|---|
| GILGFVFTL#5 | 50 | 92.0% | 13 | none | TRBV19 (50/50) |
| GILGFVFTL#1 | 31 | 100.0% | 12 | none | TRBV19 (31/31) |
| GILGFVFTL#36 | 27 | 100.0% | 13 | none | TRBV19 (27/27) |
| GILGFVFTL#27 | 27 | 100.0% | 10 | none | TRBV19 (27/27) |
| GLCTLVAML#2 | 27 | 96.3% | 11 | none | TRBV20-1 (27/27) |
| GILGFVFTL#30 | 23 | 91.3% | 19 | none | TRBV19 (23/23) |

**The confounding check comes out clean.** No family shares a single reference
across all its members, and each is reported independently in **10–19** distinct
studies. Across all families, pure ones share a reference at 17.2% against mixed
ones at 12.5% — 5/29 against 2/16, no meaningful provenance effect. So the extra
agreement is not a study or donor artifact.

**Settled.** There is **no evidence for a per-family functional switch in
general**: E20's purity pattern is fully reproduced by the null. There *is*
extra agreement, it is real and not provenance, and it lives entirely in six
large germline-restricted public-motif families — five GILGFVFTL built on
TRBV19 and one GLCTLVAML on TRBV20-1, 96.2% validated against 85.9% in the other
39 families. The "something extra each family shares" is **public-motif
identity at the V-gene level**, a known and named thing rather than a hidden
per-cluster state.

This connects to E49 from a new direction. TRBV Simpson diversity rank-ordered
typicality's per-epitope AUC perfectly, and GILGFVFTL was the extreme (0.570);
here its four largest families are all TRBV19, which is what that concentration
*is*. But the connection does not make an epitope easy: GLCTLVAML's one large
family is uniformly TRBV20-1 at 96.3% validated, and GLCTLVAML still has the
lowest nn3-alone AUC (0.677) and the weakest published band.

Nothing here changes the predictor; it is a test of the biological reading.


**Amended by E58 (independence of the six families).** E58 compared the
families' CDR3β consensus sequences and found four of the five GILGFVFTL
families share length 13 and the prefix `CASSIRS`, differing by a median of 2
residues out of 13. They are sub-clusters of a single public motif that the
0.12 radius happened to split, not five independent objects. Everything above
that counts families — "six families hold 2,996 of 3,666 pairs", the p = 0.006
on the 20+ band — therefore rests on an **effective two or three independent
motifs**, not six. The direction of the result is unchanged; its weight is
lower than a count of six families suggests.

---

## E55 — How much of "singleton" is sample size?

E53 showed singleton status is predictable from the TCR itself — CDR3 length,
Pgen, interface confidence. That says *which* TCRs lose their neighbours first.
It does not say whether they would have neighbours in a larger sample, and the
AIRR literature is explicit that publicness tracks sample size, donor count and
sequencing depth. E21 already thinned the public cohorts but measured the AUC;
the statistic that answers this directly was never extracted.
`scripts/53_singleton_is_sampling.py`, 200 draws per size, graph rebuilt on the
thinned **full** cohort with the positive/negative split read off the labelled
survivors.

### The decisive comparison: public cohorts thinned to the private sizes

| private epitope | n | actual singleton rate | public thinned to that n | verdict |
|---|---:|---:|---|---|
| LTDEMIAQY | 130 | 85.4% | GLC 50.5%, GIL 69.7%, NLV 89.4%, YLQ 90.4%, TTD 97.6% | **inside** |
| LLWNGPMAV | 229 | 95.6% | GIL 63.2%, NLV 85.4%, YLQ 85.9%, TTD 96.4% | **inside** |
| CINGVCWTV | 243 | 97.5% | GIL 62.5%, NLV 85.3%, YLQ 85.8%, TTD 96.3% | outside, by **1.2 pp** |

**Two of three private cohorts are no sparser than a public cohort cut to the
same size, and the third exceeds the public range by 1.2 percentage points.**
The anchor throughout is TTDPSFLGRY, which at its own size is already 95.3%
singleton and at 243 reaches 96.3% — indistinguishable from CINGVCWTV's 97.5%.
And TTD's nn3-alone AUC is **0.860** (E53), so a cohort this sparse can still be
ranked well.

### But sample size is not the main thing

| | spread |
|---|---:|
| across epitopes at a fixed n = 130 | **47.0%** |
| within an epitope across the whole sweep (80 → full) | 12.1% |

**Which epitope it is moves the singleton rate 3.9× further than an eightfold
change in cohort size.** So "it is just sampling" is the wrong summary: sampling
fully accounts for the private cohorts' *level*, while the between-epitope
differences it cannot explain are far larger.

### Positives and negatives are structurally different

Going from full size to n = 130, the singleton rate rises **+0.144** for
positives and **+0.038** for negatives — positives start far more clustered and
so have more to lose. But conditioning on having *had* a neighbour reverses it:

| epitope | of positives that had a neighbour, lost it at n=130 | of negatives |
|---|---:|---:|
| GILGFVFTL | 0.357 | 0.597 |
| NLVPMVATV | 0.382 | 0.620 |
| YLQPRTFLL | 0.456 | 0.702 |
| GLCTLVAML | 0.093 | 0.119 |
| **mean** | **0.363** | **0.509** |

**Negatives lose their neighbour 1.40× more readily than positives.** A positive
sits in a family with redundant membership, so subsampling has to remove several
members before it is isolated; a negative's neighbour is more often a single
chance edge that one deletion cuts. That is a real structural asymmetry and it
supports the convergent-family reading over chance proximity — and it is the
same object E54 found, where the six large TRBV19/TRBV20-1 families carry all
the excess within-family agreement.

### Settled

Being a singleton *at the observed sample size* is largely a sampling statement
for the private cohorts — their levels sit inside or 1.2 pp outside what thinned
public cohorts produce. It is **not** evidence that those TCRs have no relatives
in a real repertoire. But the robustness of a TCR's neighbourhood to
subsampling is informative and differs by label, which is the continuous
`nn3_dist` story of E53 again rather than a new feature.

**Caveat, stated because the design cannot escape it:** this subsamples a finite
candidate library. It cannot reconstruct a real antigen-specific repertoire, so
it bounds the sampling explanation **from below**, not from above.

Figure: `reports/singleton_vs_cohort_size.html`. Nothing here changes the
predictor.

---

## E56 — Provenance audit: the gate on the study-aware plan (E56–E59)

Implements the first decision gate of the supplied plan
(`IMMREP26_E56-E59_Study-aware_TCR_Convergence_Analysis.md`), which states that
if sources are insufficient or inseparable, the later experiments must be
reported **not estimable** rather than run on defaults.
`scripts/54_provenance_audit.py`. Script numbering is offset from the plan's
(`54_*` implements E56.1); the predictor is untouched.

2,472 of 2,636 competition TCRs (93.8%) match a VDJdb record on
(cdr3 alpha, cdr3 beta, epitope), giving 3,638 records.

**The first version of this gate asked the wrong question** and failed the whole
plan because `meta.study.id` is only 26.4% populated. That is the wrong field.
The independence question the plan actually poses — *a clonotype in 12 papers is
not 12 rearrangements* — is answered by **donor**, which is stricter than a
study identifier, and `meta.subject.id` is well populated where it matters.

**Per epitope, can independent donors be resolved?**

| epitope | split | subject populated | donors | refs | largest ref share | |
|---|---|---:|---:|---:|---:|---|
| GILGFVFTL | public | 91.3% | **158** | 21 | 20.2% | yes |
| NLVPMVATV | public | 89.3% | 43 | 24 | 54.3% | yes |
| GLCTLVAML | public | 66.4% | 20 | 14 | 30.2% | yes |
| YLQPRTFLL | public | 38.4% | 39 | 10 | 42.9% | partial |
| TTDPSFLGRY | public | 35.4% | 8 | 2 | 61.8% | partial |
| CINGVCWTV | private | 30.4% | 5 | 2 | 69.6% | partial |
| LLWNGPMAV | private | **0.0%** | 0 | 3 | 77.8% | **no** |
| LTDEMIAQY | private | **0.0%** | 0 | 2 | **97.7%** | **no** |

**E54's six families, recounted in donors rather than papers.**

| family | refs | subject field | donors |
|---|---:|---:|---:|
| GILGFVFTL#30 | 17 | 96.1% | **110** |
| GILGFVFTL#36 | 14 | 97.5% | 68 |
| GILGFVFTL#5 | 13 | 92.0% | 67 |
| GILGFVFTL#27 | 10 | 96.0% | 52 |
| GILGFVFTL#1 | 12 | 88.5% | 50 |
| GLCTLVAML#2 | 11 | 60.7% | **12** |

**This strengthens E54 rather than qualifying it.** E54 argued the six families
are not a provenance artifact from 10–19 distinct references; the five
GILGFVFTL families in fact draw on **50–110 distinct donors**. Only 7.8% of all
clonotypes carry more than one reference, and 80.9% of those resolve to ≥2
donors, so reference counts are not systematically inflated by shared raw data.
GLCTLVAML#2 at 12 donors and 60.7% coverage is the one family where the
argument is thin, and it must be reported separately rather than pooled.

### The verdict, by object

| item | status |
|---|---|
| E57 cross-source persistence, **donor** as the unit | **estimable** on all six families; GLC#2 separately |
| E56.2 study-block resampling **as literally written** | **not estimable** — study id 26.4% populated; substitute donor blocks and say so |
| E59 matched control on **GILGFVFTL** | **estimable** — 278 labelled TRBV19 TCRs outside the big families |
| E59 matched control on **GLCTLVAML** | **not estimable** — 13 controls; loosening the matching to rescue it is what the plan forbids |
| any provenance-aware method on the **private split** | **not estimable** — LLW and LTD have zero donor coverage and LTD draws 97.7% of its records from one reference |

### A result that falls out of the feasibility probe

Decomposing GILGFVFTL's labelled TCRs three ways:

| group | n | validated |
|---|---:|---:|
| in a large TRBV19 family | 158 | **96.2%** |
| TRBV19, not in a large family | 278 | **74.5%** |
| other V gene, not in a large family | 139 | **15.8%** |

**Most of the work is done by the V gene, not by family membership.** Being
TRBV19 moves validation from 15.8% to 74.5%; being in one of the large TRBV19
families adds a further 21.7 points to 96.2% (Fisher p 6.4e-10, **uncontrolled
for density**, which the plan correctly identifies as a mediator rather than a
confounder to be regressed away blindly).

### Standing

The chain is runnable on GILGFVFTL and nowhere else. It cannot touch the
private split, so it carries no competitive value — which the plan states as
its intent. Nothing here changes the predictor.


**Amended by E58.** The donor counts above (50–110 per GILGFVFTL family) are
correct per family, but the five GILGFVFTL families are sub-clusters of one
`CASSIRS` motif. Donors are shared across those sub-clusters, so the five
rows must not be summed or treated as five independent replications of the
same finding.

---

## E57 — Do the six public families persist across independent donors?

E56 settled the unit: `meta.study.id` is unusable at 26.4% populated, while
`meta.subject.id` reaches 91.3% on GILGFVFTL. So the plan's leave-one-study-out
becomes **leave-one-donor-out**, which is stricter — two papers can share a
cohort, two donors cannot share a rearrangement event.
`scripts/55_cross_donor_persistence.py`.

Donor = `reference.id` + `meta.subject.id`. Two exclusions, both from the plan:
a neighbour sharing any donor with the query is not independent evidence, and
an identical αβ clonotype is the same receptor recorded twice rather than
convergence (6 such pairs exist in the cohort and never count). Every statistic
is computed **conservatively** (a TCR with no resolvable donor cannot supply
support) and **permissively** (missing donor counts as a different donor); the
truth is bracketed by the two.

**E57.1 — cross-donor support, against the controls that make it meaningful.**

| group | n | median neighbours | has cross-donor support (cons / perm) | median donors supporting |
|---|---:|---:|---|---:|
| GIL large TRBV19 families | 205 | 13 | **97.1% / 99.5%** | 19 |
| GIL TRBV19, no large family | 318 | 0 | 24.2% / 26.7% | 0 |
| GIL other V, no large family | 147 | 0 | 1.4% / 2.7% | 0 |
| GLCTLVAML#2 (12 donors, thin) | 32 | 8 | 90.6% / 100.0% | 4 |

The control does its job: this is not merely "families are dense". Holding the
epitope fixed, TRBV19 outside a large family reaches 24.2% against the families'
97.1%, and other V genes 1.4%. The V-gene gradient from E56 reappears.

**E57.2 — leave-one-donor-out. No family depends on any single donor.**

| family | donors | donors tested | TCRs keep a neighbour | worst donor |
|---|---:|---:|---:|---:|
| GILGFVFTL#30 | 112 | 112 | 100.0% | 100.0% |
| GILGFVFTL#36 | 69 | 69 | 100.0% | 100.0% |
| GILGFVFTL#5 | 68 | 68 | 100.0% | 100.0% |
| GILGFVFTL#27 | 53 | 53 | 100.0% | 100.0% |
| GILGFVFTL#1 | 51 | 51 | 99.8% | **90.9%** |
| GLCTLVAML#2 | 12 | 12 | 100.0% | 100.0% |

Removing every record contributed by any one donor leaves the family's TCRs with
a neighbour essentially always, and the single worst case across 365 donor
deletions still leaves 90.9%. **These motifs are reproduced by different people,
not propagated from one cohort.**

**E57.3 — and one genuinely new thing, after discarding a confounded version.**

The naive comparison — cross-donor support against validation, 94.6% vs 48.3%
on GIL — **is confounded and was not kept**: having cross-donor support and
having any neighbour at all are nearly the same condition (median 13 against 0),
so it restates density. Conditioning on *already having a neighbour* asks the
real question, whether that neighbour coming from another person adds anything:

| epitope | cross-donor neighbour | only same-donor neighbours | odds ratio | Fisher p |
|---|---|---|---:|---:|
| GILGFVFTL | **94.6%** (n=223) | **51.4%** (n=35) | 16.6 | 5.9e-10 |
| GLCTLVAML | **94.0%** (n=50) | **60.0%** (n=40) | 10.4 | 1.5e-04 |

**Among TCRs that already have relatives, whether those relatives come from a
different person separates validation at odds 10–17.** A TCR whose only
neighbours come from the same donor is not showing convergence — it is a clonal
expansion or a within-experiment duplicate — and density cannot tell the two
apart. This is the sharpest mechanistic result the project has on why density
works and where it misfires.

**Settled.** The six families are supported across 12–112 independent donors,
survive the deletion of any one of them, and the V-gene control shows the effect
is motif-level rather than germline-level alone. Combined with E56, the "these
are reproduced across studies" claim now rests on donors rather than on paper
counts.

**It cannot be used competitively**, and that is not a judgement call: E56
measured 0% donor coverage on LLWNGPMAV and LTDEMIAQY and 30.4% on CINGVCWTV,
with LTD drawing 97.7% of its records from a single reference. The quantity that
separates real convergence from clonal duplication is unavailable on exactly the
rows that decide the ranking. The entry is unchanged.


**Amended by E58.** "Five of six families show cross-donor support" is five
sub-clusters of about two motifs. The leave-one-donor-out result is unaffected,
because it is a within-family statement tested 365 times at the donor level,
and so is E57.3's odds ratio, which is computed over all GILGFVFTL TCRs rather
than per family. What must be discounted is any reading that treats the six
family-level rows as six independent confirmations.

---

## E58 — Are the public motifs enriched beyond what rearrangement explains?

E57 showed the six families are reproduced across 12–112 independent donors.
That says they are real observations, not one cohort propagated; it does not say
why they recur. V(D)J recombination makes some CDR3s far easier to produce, and
a motif that is simply easy to generate will be found repeatedly with no antigen
preferring it. `scripts/56_generation_adjusted_motif.py`.

**Two limits fixed before the numbers.** OLGA generates CDR3s without a V gene,
so a motif here is defined on **CDR3β only** and the families' TRBV19
restriction is reported separately rather than folded into the enrichment. And
*any* cohort motif beats a naive generation model, because the cohort is
antigen-selected to begin with — so the informative number is not the raw
enrichment but the comparison with the plan's required control: same epitope,
same TRBV19, matched size, not in a large family.

The motif definition is fixed before any enrichment is computed and uses no
labels: modal CDR3β length, position-wise consensus, radius = rounded-up median
within-family Hamming distance to that consensus.

**Part 1 — generation explains part of it.** Within GILGFVFTL:

| group | n | median log10 Pgen(β) | median CDR3β length |
|---|---:|---:|---:|
| in a large TRBV19 family | 205 | **−9.19** | 13 |
| TRBV19, no large family | 318 | −9.74 | 13 |
| other V, no large family | 147 | −10.26 | 14 |

Family CDR3βs are **+0.555 log10 units (3.6×) easier to generate** than matched
TRBV19 non-family ones at the same length. So "easy to produce" is part of why
they recur, and the enrichment below must be read against that.

**Part 2 — enrichment over a 200,000-sequence OLGA background.**

| family | L | r | consensus | cohort | OLGA | enrichment |
|---|---:|---:|---|---:|---:|---:|
| GILGFVFTL#36 | 13 | 1 | CASSIRSAYEQYF | 0.1597 | 3.5e-05 | **4563×** |
| GILGFVFTL#27 | 13 | 1 | CASSIRSTDTQYF | 0.0806 | 2.0e-05 | 4030× |
| GILGFVFTL#30 | 13 | 1 | CASSIRSSYEQYF | 0.2060 | 5.5e-05 | 3745× |
| GLCTLVAML#2 | 13 | 3 | CSARDGTGNGYTF | 0.1250 | 7.5e-05 | 1667× |
| GILGFVFTL#1 | 12 | 1 | CASSIGSYGYTF | 0.0418 | 4.0e-05 | 1045× |
| GILGFVFTL#5 | 13 | 3 | CASSIRSTGEQYF | 0.4090 | 5.5e-04 | 744× |

A first version scored **every** family against GILGFVFTL's cohort, which gave
the GLCTLVAML motif an observed frequency of exactly zero and an enrichment of
0.0. Fixed: a motif's observed frequency is measured in its own epitope's
cohort.

**Part 3 — against the matched control, which is the number that counts.**

| family | enrichment | control median | control p95 | p |
|---|---:|---:|---:|---:|
| GILGFVFTL#36 | 4563× | 34× | 130× | 0.005 |
| GILGFVFTL#27 | 4030× | 37× | 314× | 0.005 |
| GILGFVFTL#30 | 3745× | 34× | 130× | 0.005 |
| GILGFVFTL#1 | 1045× | 38× | 369× | 0.005 |
| GILGFVFTL#5 | 744× | 34× | 79× | 0.005 |
| GLCTLVAML#2 | 1667× | — | — | **not estimable** |

The control confirms the warning: a random same-size set of TRBV19 non-family
GIL sequences already shows ~34× enrichment over OLGA. The real motifs run
744–4563×, 20–130× above the control median, and all five estimable families
clear p = 0.005. **GLCTLVAML#2 is reported not estimable rather than forced**,
because E56 found only 13 matched controls and loosening the matching to rescue
it is what the plan forbids.

**The qualification that matters most, and it came out of the output itself.**
The five GILGFVFTL families are separate only because the 0.12 radius split
them. Their consensus sequences:

```
CASSIRSTGEQYF   .   L   2   2   2
CASSIGSYGYTF    L   .   L   L   L
CASSIRSTDTQYF   2   L   .   3   3
CASSIRSAYEQYF   2   L   3   .   1
CASSIRSSYEQYF   2   L   3   1   .
```

Four of the five share length 13 and the prefix `CASSIRS`, differing by a
**median of 2 residues out of 13**. They are sub-clusters of a single
CASSIRS-type public motif — the known TRBV19 influenza M1 response — not five
independent discoveries. **Five p-values of 0.005 are not five independent
pieces of evidence; the effective number of distinct GILGFVFTL motifs tested
here is about 2**, and every statement about "five families" in E54, E56 and
E57 carries the same qualification retroactively.

**Settled.** Enrichment beyond the rearrangement model is real and large
(20–130× above a matched control that already corrects for cohort selection),
and it is not explained by generation probability — though generation does
contribute, since family members are 3.6× more generatable than matched
controls. The evidence rests on roughly **two** independent motifs, one of them
the textbook CASSIRS/TRBV19 M1 response, plus one GLCTLVAML motif whose control
is not estimable. That is consistent with antigen-associated selection and does
not demonstrate it: the plan's own boundary, that a naive repertoire or
functional experiment is required, still stands. Nothing here changes the
predictor.

---

## E59 — Do the public motifs revalidate better than matched non-motif TCRs?

E58 measured how often the motifs occur. This asks whether occurring in one
predicts passing functional revalidation, against controls matched on the things
that could produce the same difference without any recognition advantage.
`scripts/57_functional_residual.py`. **GILGFVFTL only** — E56 found 13 matched
controls on GLCTLVAML, and the plan forbids loosening the matching to rescue it.

**The comparison**, within GILGFVFTL and restricted to TRBV19 so the V gene is
held fixed:

| group | n | validated |
|---|---:|---:|
| in a large family | 158 | **96.2%** |
| TRBV19, no large family | 278 | 74.5% |

Raw difference **+0.217**, odds ratio 8.7, Fisher p 6.4e-10.

**The matching ladder.**

| matching | difference | motif TCRs in a usable stratum | strata |
|---|---:|---:|---:|
| CDR3β length | +0.210 | 158/158 | 2 |
| length × Pgen tercile | **+0.230** | 158/158 | 6 |
| length × Pgen × **nn3** tercile | +0.082 | 137/158 | 9 |

Matching on generation probability does not reduce the advantage — it is
+0.230 after, against +0.217 before — so E58's finding that family members are
3.6× more generatable does **not** account for their validating more often.

The last row controls `nn3`, and it is reported because the plan asks for it,
**not as the causal estimate**: the motif's convergence is *why* its members
have close neighbours, so matching on density removes the mechanism along with
any confounding. That it still leaves **+0.082** is informative in the other
direction — even comparing TCRs with equally close neighbours, motif membership
carries something.

**E59.2 — the statistical unit, corrected by E58.** Leave-one-out at two levels,
because E58 showed four of the five families are sub-clusters of one motif:

| level | removals | range of the difference |
|---|---:|---|
| family | 5 | +0.208 to +0.237 |
| **motif** | 2 | +0.208 to +0.255 |

The two motif groups are `CASSIRS-13` (families #5, #27, #30, #36) and
`CASSIGSYGYTF` (family #1). Removing either leaves the advantage positive —
+0.208 without the small group, +0.255 without the CASSIRS cluster, where only
31 motif TCRs remain. So the effect is not carried by the textbook CASSIRS
cluster alone, which is the strongest statement available at an effective n of
two motifs.

**E59.3 — joint reading with E58, with the plan's boundary intact.**

Motifs that recur beyond the rearrangement model (E58: 744–4563× against a
matched control's ~34×) also revalidate more often than TCRs matched on epitope,
V gene, CDR3β length and generation probability (E59: +0.230). Both results
survive removing either motif group.

That is **consistent with antigen-associated selection and does not demonstrate
it.** The plan's boundary stands unchanged: distinguishing selection from
ascertainment quantitatively needs a naive or pre-exposure repertoire, or a
functional experiment, and the evidence here rests on about two independent
motifs within a single epitope. Nothing in E56–E59 changes the predictor, and
nothing in it can: E56 measured 0% donor coverage on two of three private
epitopes.

---

## E60 — Gate 0 of the recognition-profile roadmap: the external data audit

Implements Gate 0 of `IMMREP26_Recognition_Profile_to_Functional_Activation_Roadmap_2026-10-26.md`.
Both datasets downloaded and audited: PRP (Dryad 10.5061/dryad.h70rxwdzk, v5,
183,881,424 B) and SEEKER (10.5061/dryad.79cnp5jct, v3, 27,946,775 B) — both
**CC0-1.0**, both exactly the sizes the roadmap stated. 76 files, 273.5 MB,
manifest with sha256 at `data/derived/external_manifest.tsv`.
`scripts/58_gate0_external_audit.py`. Gitignored like TCRvdb.

### Two data quirks that would have corrupted a naive pipeline

**`TCR4.2`, `TCR4.3` and `TCR4.4` ship with no header row** — their first line
is already data — and they are *exactly* the three P2/P9-anchor files. Globbing
the directory with pandas defaults eats one peptide from each and gives them
different column names from the other eighteen.

**`TCR63` carries no `add_flow_high_copy` column.** The two SEEKER receptors do
not hold the same flow conditions, so the roadmap's F1 — fit single- against
high-copy antigen density and compare — is available on **one** receptor, not
two.

### What the audit found

| check | result |
|---|---|
| PRP anchor libraries | P2_P8 272,888 peptides, P2_P9 151,914, sharing 256 (0.17%) — genuinely different panels, kept apart |
| **PRP input coverage** | median **8.9%** of selected peptides appear in their own naive library; **0.2%** for the three P2/P9 files |
| SEEKER R0 | present for **100%** of candidates (19,322 TCR59, 16,369 TCR63) |
| SEEKER flow | 192 peptides per receptor, **~1%** of candidates, selected as top hits |
| top-10,000 tables | 6 of them; absence from one is **not** a measured negative |
| PRP ∩ SEEKER peptides | 19 of 6,313,097 and 35,374 |
| PRP ∩ SEEKER TCRs | **0** |

**The PRP input-coverage number is the one that bites.** The roadmap's
enrichment formula needs an input count in its denominator, and 91% of selected
peptides do not have one — the selections contain far more distinct peptides
(6.3M) than the two naive library files list (425k combined). For the P2/P9
receptors it is 99.8% missing. Input-corrected enrichment is therefore
computable on a small minority of the data, and any per-TCR profile built
without it is a within-receptor rank, not a calibrated quantity.

### Gate 0 verdict: **PROFILES ONLY**

| | |
|---|---|
| [FAIL] | PRP input library available for the enrichment denominator |
| [PASS] | SEEKER R0 exists |
| [PASS] | SEEKER flow measurements exist |
| [FAIL] | same TCR × same peptide × same HLA across the two sources |

Both studies use HLA-B\*27:05, so peptide overlap was possible and was measured:
**19 shared peptides and zero shared TCRs.** PRP measures yeast binding on 21
receptors (16 patient + 5 engineered variants of TCR19.2); SEEKER measures NFAT
activation on TCR59 and TCR63. **No pair carries both readouts**, so a
cross-source binding → activation map is not estimable, and the roadmap's rule
applies: report that rather than stitch the two into one label.

**What remains available**, and it is not nothing: per-receptor recognition
profiles from each source, and the within-SEEKER screen-to-flow calibration
(F1), where the same receptor has both an R1/R0 ratio and a flow readout — on
TCR59 with both copy conditions and TCR63 with single-copy only. The measured
flow subset is ~1% of candidates chosen as top hits, so any calibration fitted
on it describes ranked candidates, not the library, and must say so.

Nothing here touches the frozen predictor; the roadmap excludes it by design.

---

## E61 — F1: does a SEEKER reporter screen predict the flow activation it proxies?

E60 put Gate 0 at PROFILES ONLY, which left exactly one quantitative task open:
inside SEEKER, the same receptor carries both a screening enrichment and a flow
readout on the same peptides. `scripts/59_seeker_flow_calibration.py`.

**The conditional sampling, restated with numbers.** Flow was run on 192
peptides per receptor — 0.99% of TCR59's candidates and 1.17% of TCR63's —
chosen as top hits. Their median R1/R0 is **125.8 against 14.8** across all
TCR59 candidates, and **170.2 against 14.4** for TCR63. Everything below
therefore describes *ranked candidates*, and nothing is inferred about the
library's activation rate.

**Screening statistic against flow activation**, Spearman on the measured
subset plus 5-fold cross-validated isotonic and log-linear calibration:

| TCR | predictor | raw ρ | CV isotonic | CV log-linear |
|---|---|---:|---:|---:|
| TCR59 | **r2_r1_ratio** | **+0.687** | **+0.635** | +0.682 |
| TCR59 | count_pct_round1 | **−0.480** | +0.462 | +0.481 |
| TCR59 | r1_r0_ratio | **−0.450** | +0.448 | +0.448 |
| TCR59 | count_pct_round2 | +0.378 | +0.345 | +0.372 |
| TCR63 | count_pct_round2 | +0.537 | +0.488 | +0.521 |
| TCR63 | r2_r1_ratio | +0.324 | +0.427 | +0.320 |
| TCR63 | r1_r0_ratio | +0.143 | +0.066 | −0.187 |

**The signs are the result.** On TCR59, first-round abundance runs *negative*
against flow activation (`count_pct_round1` −0.480, `r1_r0_ratio` −0.450) while
the second-round gain runs strongly positive (`r2_r1_ratio` +0.687). **Early
enrichment is not activation strength; the later-round gain is what tracks it.**
The roadmap warned against reading a final-round abundance as binding strength,
and this measures the warning: within the top 1%, the round-one signal points
the wrong way.

A method bug was caught doing this. `sklearn.isotonic.IsotonicRegression`
defaults to an *increasing* fit, so on the negatively-correlated predictors it
returned near-zero (−0.059) while the linear model on the same data returned
+0.448. Fixed with `increasing="auto"`; the table above is after the fix.

**Single-copy against high-copy (TCR59 only — E60 found TCR63 ships no
high-copy column).**

| condition | n | median activation | ρ with r1_r0 |
|---|---:|---:|---:|
| single copy | 192 | 22.0% | −0.450 |
| high copy | 192 | 17.4% | −0.404 |

The two conditions agree at **ρ +0.923**, and **zero of 192** peptides are
inactive at single copy (<5%) and active at high copy (>25%). Median activation
is in fact slightly *lower* at high copy. **Within these top-ranked candidates,
presentation level is not what decides whether a peptide activates** — the
opposite of the risk the roadmap flagged. The claim is bounded to candidates
already selected for enrichment; it says nothing about weak or unselected
peptides, which is exactly where a presentation effect would be expected.

**Settled.** A SEEKER reporter screen does carry rank information about the flow
activation it proxies — out-of-fold ρ +0.635 on TCR59 — but only through the
*later-round gain*, and the measured subset is a top-1% slice selected on the
predictor itself, so the correlation is attenuated by range restriction and
cannot be extrapolated to the library. Antigen density does not reorder these
candidates. Both statements are about 192 peptides on two receptors under one
HLA, and neither is evidence about binding, which SEEKER does not measure.
Nothing here touches the frozen predictor.

---

## E62 — Substitution grammar: what explains how one residue change moves the response?

The question is whether a position-conditioned substitution model beats fixed
matrices (BLOSUM62, Grantham, physicochemical descriptors) at explaining
response change, under peptide-family holdout.
`scripts/60_substitution_grammar.py`.

**The substrate had to be chosen, not assumed.** Single mutants *relative to the
wild type* are too few: the Flu-M1 mutagenesis table holds **69**, covering 20
of 171 possible (position, residue) cells, with positions 2, 4 and 6 entirely
absent — and it is a top-10,000 table, so substitutions that killed the response
were never recorded. The plan's own fallback applies. What works instead is
peptide pairs *within* one screen differing at one residue:

| screen | peptides | Hamming-1 pairs | families |
|---|---:|---:|---:|
| Flu (4TCR) | 7,264 | 7,122 | 1,425 |
| EBV (4TCR) | 8,362 | 3,790 | 764 |
| 1G4 (4TCR) | 8,123 | 1,587 | 108 |
| DMF5 (4TCR) | 8,022 | 247 | 107 |
| TCR59 / TCR63 (B27) | — | **15 / 1** | unusable |

**A data trap worth recording.** These tables are ranked per **nt construct**,
not per peptide: 10,000 rows hold 8,545 distinct 9-mers and one peptide carries
up to **16 synonymous encodings**. Treating rows as peptides manufactures
duplicate pairs with a delta of exactly zero — a first run reported median
|Δ| 0.000 and 42,810 Flu pairs against the 7,122 a deduplicated count gives.
Reads are now summed over a peptide's encodings, and `r3_r2_ratio` is
recomputed from the summed counts with a declared pseudocount, because the
shipped column is NaN for the 2,785 rows unobserved in round 2.

**Part 1 — held-out Spearman by model** (grouped by Hamming-1 connected
component, response is the later-round gain per E61):

| model | 1G4 | DMF5 | EBV | Flu |
|---|---:|---:|---:|---:|
| BLOSUM62 | −0.079 | −0.006 | −0.024 | +0.013 |
| Grantham | −0.003 | +0.006 | −0.007 | −0.006 |
| physicochemical | +0.075 | −0.180 | −0.015 | +0.041 |
| position only | +0.089 | +0.191 | −0.029 | +0.009 |
| position + Grantham | +0.094 | +0.204 | −0.001 | +0.021 |
| position × chemistry | **+0.126** | +0.177 | −0.016 | **+0.062** |

**Fixed substitution matrices explain essentially nothing** — BLOSUM62 and
Grantham sit between −0.079 and +0.013 on every receptor.

**Part 2 — and the negative control removes the rest.** Holding the position
fixed and shuffling which chemical substitution it was:

| screen | position + Grantham | chemistry shuffled |
|---|---:|---:|
| DMF5 | +0.204 | **+0.209** |
| 1G4 | +0.094 | +0.080 |
| Flu | +0.021 | +0.012 |
| EBV | −0.001 | −0.032 |

DMF5's apparent +0.204 is **identical** after the chemistry is destroyed. What
the position-conditioned models learned is position, and nothing about which
residue replaced which.

**Part 3 — the cross-receptor split, and why it is not estimable.** Four
receptors on one NNK library allow the decomposition the plan asks for: an
effect shared across receptors is library or presentation, one that differs is
recognition. On the 80 (position, a, b) cells measured on all four, mean
cross-receptor ρ is **+0.169**, which looks like receptor-specificity.

It is not, and the check that settles it is the noise ceiling — does a receptor
agree with *itself*? Splitting each cell's pairs in half, 50 repeats:

| | split-half ρ |
|---|---:|
| DMF5 | +0.360 |
| 1G4 | +0.332 |
| Flu | +0.314 |
| EBV | +0.094 |
| **mean reliability** | **+0.275** |
| mean cross-receptor agreement | +0.169 |

**A receptor barely agrees with itself.** The cross-receptor number sits just
below a ceiling of 0.275, so the per-cell estimates are noise-dominated and
**nothing can be concluded about receptor-specificity**. Without this check the
result would have read as a clean finding.

**Settled.** The substitution grammar is not recoverable from these screens at
this noise level. That is a statement about the data, not about the models:
fixed matrices genuinely carry nothing here, position carries a little, and the
decomposition that would have been interesting cannot be made because single
(position, a, b) cells are estimated from too few pairs. Recovering it would
need replicate screens or far deeper sampling per cell, not a better model.
Nothing here touches the frozen predictor.


**SUPERSEDED IN PART by the 2026-10-27 code review.** The split-half ceiling
below (0.275) lacked a Spearman–Brown correction and the cross-receptor
correlation was compared against it raw, although that correlation is
attenuated by both receptors' noise. Corrected: ceiling **0.419**, and
disattenuating per pair while dropping EBV (self-reliability 0.171) gives a
shared component of **46%**, so about **54% of a substitution's effect is
receptor-specific**. **"Nothing can be concluded about receptor-specificity"
is withdrawn** — the honest statement is roughly half shared, half specific,
with wide uncertainty. Everything else in this entry stands: fixed matrices
explain nothing, and the position-conditioned models do not survive the
chemistry-shuffled control.

---

## E63 — Recognition distance against sequence distance (PRP)

PRP's Dryad package carries peptide counts and nothing else. The TCR sequences
are in the companion GitHub repo as `data/tcr_cdr3b.csv`, 487 bytes, and
reading it produced the finding this experiment is built on:

```
CASSVAYYSTDTQYF  ->  135.7, 19.2C1, 28.1
CASSVATYSTDTQYF  ->  19.2C2, 4.3
CASSVGLYSTDTQYF  ->  19.2C5, 4.1
CASSVGTYSTDTQYF  ->  19.2C4, 8.4
```

**Twenty-one screened receptors carry sixteen distinct CDR3β.** The engineered
19.2 variants were built carrying another patient TCR's β, so five comparable
pairs have a CDR3β distance of exactly **zero** while being different
receptors. `scripts/61_recognition_vs_sequence.py`.

**The limitation that bounds everything below: only CDR3β exists here** — no
alpha, no V/J. "Sequence distance zero" means zero *on the axis that can be
measured*.

**Noise floor first**, learnt from E62. Splitting each receptor's round-4 reads
by a binomial draw and rebuilding the profile from each half gives a within-
receptor reliability of **0.687** (range 0.178–0.990). Profiles here are
reliable, unlike E62's per-cell estimates.

**The five identical-CDR3β pairs**, against 151 others with ≥5,000 common
peptides:

| a | b | common peptides | profile ρ | ceiling | ρ / ceiling |
|---|---|---:|---:|---:|---:|
| 135.7 | 28.1 | 46,940 | 0.304 | 0.788 | 0.386 |
| 19.2C5 | 4.1 | 53,566 | 0.269 | 0.671 | 0.401 |
| 19.2C4 | 8.4 | 38,547 | 0.259 | 0.546 | 0.475 |
| 135.7 | 19.2C1 | 57,972 | 0.246 | 0.656 | 0.376 |
| 19.2C1 | 28.1 | 43,804 | 0.218 | 0.580 | 0.375 |

| | median profile ρ |
|---|---:|
| identical CDR3β (5 pairs) | **0.259** |
| everything else (151 pairs) | **0.282** (range 0.139–0.482) |

**Identical CDR3β buys −0.023.** Receptors carrying literally the same β chain
recognise no more similarly than a random comparable pair.

**And sequence distance explains nothing across the board**: CDR3β Levenshtein
against profile distance gives ρ **+0.015** (n=156); Hamming gives +0.010.

Counterexamples exist in both directions. Sequence close, recognition far:
19.2C3 vs 26.1 at Levenshtein 1 and ρ 0.145. Sequence far, recognition close:
4.2 vs 4.4 at Levenshtein 5 and ρ **0.482**, the highest in the set and 0.70 of
the ceiling.

**Settled.** On this dataset CDR3β distance carries no information about
recognition-profile similarity, and even a zero distance does not predict
similar recognition. That confirms the roadmap's H0 in the direction it
predicted — but the correct statement is **"CDR3β alone is insufficient"**, not
"sequence distance fails". The divergence between receptors sharing a β could
be entirely an alpha effect, and alpha is unobserved. Twenty-one receptors, one
HLA, one disease context.

---

## E64 — E61's calibration under family holdout, and does model shape help?

E61 used a plain KFold, sound only if nearby peptides cannot land on both sides
of a split. `scripts/62_flow_calibration_v2.py`.

**Measuring the concern before fixing it.** Among the 192 flow-measured
peptides there is **1 Hamming-1 pair on TCR59 and 4 on TCR63**, and only 68 and
80 are even 9-mers — these are real proteome and virome peptides, not an NNK
library. **Sequence-family leakage was never live for E61.** The groupings that
do exist are `Gene` (160 and 157 distinct) and `origin` (human against virus,
which makes a real domain shift rather than a fold).

**Held-out Spearman of `r2_r1_ratio` against single-copy flow:**

| TCR | model | plain KFold | gene holdout | origin holdout |
|---|---|---:|---:|---:|
| TCR59 | **linear (log x)** | +0.682 | **+0.683** | +0.680 |
| TCR59 | two-piece linear, learnable knot | +0.682 | +0.683 | +0.680 |
| TCR59 | isotonic, direction from data | +0.635 | +0.639 | +0.665 |
| TCR59 | saturating sigmoid | +0.578 | +0.588 | +0.657 |
| TCR63 | **isotonic, direction from data** | +0.427 | +0.424 | **+0.514** |
| TCR63 | saturating sigmoid | +0.386 | +0.405 | +0.368 |
| TCR63 | linear (log x) | +0.320 | +0.309 | +0.319 |

**E61's number stands.** TCR59 goes +0.682 → +0.683 → +0.680 across the three
schemes; TCR63 goes +0.427 → +0.424 → +0.514. Holding out whole source genes
changes nothing, and training on human peptides to test on viral ones changes
nothing on TCR59 and *helps* TCR63.

**Model shape buys nothing.** The two-piece model with a learnable knot returns
*exactly* the linear numbers on TCR59 — its second slope fits to about zero, so
no knot is supported by the data. The saturating sigmoid is worse than linear
everywhere on TCR59. On TCR63 isotonic wins, but linear wins on TCR59, so there
is no consistent shape. 192 points selected on the predictor itself cannot
identify a four-parameter curve, which the roadmap said in advance.

**Settled.** The screen-to-flow relationship is a monotone trend with no
evidence of a threshold or saturation in the measured range, and E61's estimate
was not inflated by leakage. Getting further would need flow on peptides that
were *not* chosen for being top-ranked, which no public dataset here provides.
Nothing in E62–E64 touches the frozen predictor.


**SUPERSEDED IN PART by the 2026-10-27 code review.** The two-piece knot model
was never fitted: its `position` parameter started at 2.0, outside the observed
range of log10(r2_r1_ratio) ([−0.44, 1.13]), so `clip()` was zero everywhere,
`slope2` had no gradient and stayed at its initial −10 in all five folds, and
the model was a straight line by construction. **"Its second slope fits to
about zero, so no knot is supported" is withdrawn** — the slope never fitted.

With the knot bounded inside the data, TCR59 stays linear (+0.683 against the
knot's +0.658) but **TCR63 has a real threshold**: knot **+0.581** under gene
holdout and **+0.632** across the human/virus split, against linear's +0.309,
with the knot at log10(r2_r1) ≈ 0.33 in every fold. TCR63's best held-out ρ is
therefore **0.632**, not the 0.427 tabled below. The leakage finding — that
E61's plain KFold was not inflating anything — is unaffected.

---

## Code review of E62–E64 (2026-10-27)

Three bugs, all in the same part of the machinery — the noise models. Two of
them changed a conclusion.

### 1. E64: the two-piece knot was never fitted — **conclusion reversed**

`knot(x, intercept, slope1, slope2, position)` was handed `p0 = [20, 10, -10,
2.0]` with no bounds. The observed range of `log10(r2_r1_ratio)` is
**[−0.44, 1.13]** on TCR59, so a knot at 2.0 sits outside the data,
`np.clip(x − 2.0, 0, None)` is zero everywhere, `slope2` has no gradient, and
it stayed at exactly −10.0 in all five folds. **The model was a straight line
by construction**, which is why it returned the linear numbers to three
decimals — and why the write-up said "its second slope fits to about zero, so
no knot is supported". The slope never fitted.

Refitting with `position` bounded to the 10th–90th percentile of x:

| TCR | linear | knot, bounded | gain |
|---|---:|---:|---:|
| TCR59 | **+0.683** | +0.658 | −0.026 |
| TCR63 | +0.309 | **+0.581** | **+0.272** |

Under the human↔virus split TCR63's knot reaches **+0.632**, and the knot lands
at log10(r2_r1) ≈ **0.33** in every fold. **The corrected conclusion is the
opposite of the one logged**: shape buys nothing on TCR59, which is linear, and
a great deal on TCR63, which has a genuine threshold the broken model hid. The
E64 entry above is superseded on this point; TCR63's best held-out ρ rises from
0.427 to **0.632**.

### 2. E62: the ceiling was uncorrected and the comparison unadjusted — **conclusion softened**

Two faults in one comparison. The split-half reliability used half the pairs
per cell with no Spearman–Brown correction, understating full-depth
reliability; and the cross-receptor correlation was compared against it raw,
although that correlation is attenuated by **both** receptors' noise.

| | as logged | corrected |
|---|---:|---:|
| mean within-receptor reliability | 0.275 | **0.419** |
| cross-receptor agreement (raw) | 0.169 | 0.169 |
| disattenuated, per pair | — | 1G4–Flu 0.641, DMF5–Flu 0.381, 1G4–DMF5 0.363 |

EBV's own reliability is 0.171 even after correction, so it is **dropped**
rather than disattenuated — dividing by 0.17 manufactures numbers from noise,
and it would otherwise have contributed a spurious DMF5–EBV value of 0.704.

Among the three usable receptors the shared component is **46%**, so roughly
**54% of a substitution's effect is receptor-specific**. E62 concluded "nothing
can be concluded" — **too pessimistic**. The honest statement is that about half
the effect is shared and half is receptor-specific, with wide uncertainty,
since dividing by a reliability near 0.5 doubles the estimate and its error
bars together. What stands unchanged: fixed matrices explain nothing, and the
position-conditioned models survive only position, not chemistry.

### 3. E63: the ceiling shared its denominator — **conclusions unaffected**

Both halves of the split divided by the *same* Round-1 term, which is perfectly
correlated between them. Measured across eight receptors that inflates the
ceiling by **+0.221** (0.652 against 0.431 when Round 1 is split too). The
missing Spearman–Brown correction pulled the other way, so the reported 0.687
was wrong in method but close in value; the corrected ceiling is **0.527**.

**None of E63's conclusions use the ceiling.** Identical CDR3β still gives
median ρ 0.259 against 0.282 for all other pairs, and CDR3β Levenshtein against
profile distance is still ρ **+0.015** (n=156). Only the `ρ / ceiling` column
and one narrative sentence change. One receptor has a *negative* split-half
reliability after correction, i.e. its profile is noise, which the inflated
version had hidden.

### What this says about the method, not the data

All three are the same mistake in different clothes: **a noise model that was
not itself checked.** E62 was nearly misread in one direction and E64 in the
other, both because a reliability number was taken at face value. The rule
worth keeping: when a conclusion rests on comparing an observed correlation
with a ceiling, the ceiling needs its own validation — split the right thing,
correct for the split, and drop any unit whose self-agreement is too low to
divide by.

---

## Full code review of E49–E64 (2026-10-27)

Prompted by the three bugs found in E62–E64: if that class were widespread,
many conclusions would be in doubt. Every experiment back to E49 was scanned
by bug class rather than by reading each script end to end, since the three
bugs shared a shape.

| bug class | scripts with the pattern | real bugs | conclusions changed |
|---|---:|---:|---:|
| unbounded optimizer parameter | 3 | **1** (E64) | 1 |
| noise ceiling not validated | 5 | **2** (E62, E63) | 1 |
| `IsotonicRegression` default direction | 2 | 0 remaining | 0 |
| permutation null on a different unit than the statistic | 10 | 0 new | 0 |
| asymptotic p at small n | 5 | 1 immaterial (E53) | 0 |
| resamples silently dropped, biasing a null | 1 checked | 0 | 0 |
| independence inflation of the analysis unit | 2 checked | 0 new | 0 |
| paired vs unpaired resampling | 1 checked | 0 | 0 |

**What was checked and cleared, with the reason:**

- **The two other `minimize` hits are not this bug.** `scripts/32` and
  `scripts/34` call OpenMM's `simulation.minimizeEnergy(maxIterations=...)` —
  a physics minimiser with an iteration cap, not a statistical curve fit with
  free parameters that can stick at their initial values.
- **E50's "ceiling" is an oracle bound, not a noise model.** It is the best
  score achievable *given the labels*, computed exactly, so there is nothing to
  disattenuate. The same is true of E49's and E51's oracle numbers.
- **E54's permutation is same-unit.** Observed and null both go through
  `purity_counts(·, codes, n_families)` with identical codes, which is exactly
  the E41→E46 failure avoided.
- **E45's stratified shuffle** permutes assignment indices within each stratum,
  the correct construction.
- **E58's control is the same code path.** Pseudo-motifs go through
  `motif_from` → `match_rate(pool)` → `match_rate(cohort)` identically to the
  real motifs. And the `isfinite` filter that could have biased the control
  distribution downward by dropping zero-background draws **never fires**:
  measured 0/200 dropped at every motif size.
- **E49's cluster bootstrap is paired.** `raw` and `cond` are computed on the
  same resample before differencing, so the CI is on the drop rather than on
  two independently resampled quantities.
- **E59's weighting is standard Mantel–Haenszel**, `n_a·n_b/(n_a+n_b)`.
- **E53's n=8 correlation** used scipy's asymptotic p with only three distinct
  values in y. Recomputed exactly over all 40,320 permutations: **0.510 against
  the reported 0.457**. Both null, so the "enrichment does not predict
  difficulty" conclusion is unaffected; the script now reports the exact p.
- **E63's five identical-CDR3β pairs are not five independent observations** —
  {135.7, 19.2C1, 28.1} is a triple contributing three of them. Collapsing it
  to one value leaves the median at **0.259**, unchanged to three decimals.
- **E39's Mann–Whitney** runs at n in the hundreds, where the asymptotic p is
  appropriate.

**Why the three bugs clustered where they did.** They all live in machinery
that only appears in E62–E64: reliability ceilings and free-form curve fitting,
both introduced for the external datasets. The IMMREP26 experiments (E49–E59)
are built from direct AUC measurement, oracle upper bounds and same-unit
permutations — a smaller and more robust toolkit, and one whose central output
has been validated against the leaderboard **four times** (0.8396→0.839,
0.8498→0.849, 0.8601→0.860, 0.8557→0.855). The competition conclusions are the
least exposed to this bug class, not the most.

**The three bugs were real but not the tip of an iceberg.** Their shared
lesson stands as a standing rule: *a noise model is a result and needs its own
validation.* Split the right thing, correct the split for its own halving, and
drop any unit whose self-agreement is too low to be a denominator.

---

## E67 / E68 — Configuration D's error anatomy, and whether any of it is label-free

Commissioned by `IMMREP26_E62-E64_Audit_and_E67-E70_Outlier_Roadmap.md`. No
model is trained and the submitted file is untouched.
`scripts/63_outlier_anatomy.py`. D's four components are rebuilt and verified
against the submitted CSV at **r = 0.999999** before anything is measured.

**The hazard the brief flagged is real and is honoured.** This project uses two
different 0.12 thresholds — the single-linkage *edge* radius that defines
families, and `nn3_dist`, the mean distance to the three nearest members.
Cluster state here is the edge radius on the **full cohort** including
unlabelled rows, as in E54 and E57; `nn3_dist` appears only as a continuous
covariate and is never thresholded.

### E67 — where D's remaining loss sits

| class | share of pairs | mean error rate | share of loss |
|---|---:|---:|---:|
| positive singleton × negative singleton | 49.4% | 0.192 | **64.0%** |
| positive singleton × negative cluster | 3.9% | **0.634** | 19.3% |
| cluster × cluster | 5.6% | 0.215 | 10.2% |
| positive cluster × negative singleton | 41.1% | 0.026 | 6.6% |

**D did not change the shape of the error.** E22 put `both_singleton` at 62% of
loss on the frozen entry and E50 at 59.8% under the label-free partition; D
gives 64.0%. The adversarial class is still ~19% of loss at a 63% error rate.
Buying +0.0102 of macro AUC moved the level, not the structure.

**Pair weight against independent TCRs**, which the brief rightly insists on
separating:

| epitope | misranked pair-slots | TCRs | carried by the worst 5% of TCRs |
|---|---:|---:|---:|
| GILGFVFTL | 16,820 | 575 | 34.7% |
| YLQPRTFLL | 4,602 | 382 | **43.2%** |
| NLVPMVATV | 5,220 | 317 | 32.7% |
| TTDPSFLGRY | 8,806 | 342 | 23.0% |
| GLCTLVAML | 2,996 | 166 | 17.7% |

The loss is concentrated in few receptors on four of five epitopes — but **not
on GLCTLVAML**, where it is the most evenly spread. GLC is not failing because
of a handful of outliers; it is failing broadly.

Casebook of the worst 40 at `reports/38_outlier_casebook.tsv`: **26 high-scoring
negatives against 14 low-scoring positives**, and 26 of the 30 GILGFVFTL cases
are negatives.

### E68 — is the failure visible without labels? Four gates, and it dies on the last two

`U = SD(rank_density, rank_pool, rank_Vdomain, rank_structure)`, computable on
the private split.

**Gate A** — misranked pairs rise with U on GLC (12.5 → 24.1 across quartiles)
and NLV, but not monotonically on GIL or TTD, and flat on YLQ.

**Gate B — passes.** Poisson regression of the misranked-pair count on D's own
score, `nn3`, cluster state and epitope, with standard errors clustered by
family: **U coefficient +1.148, p = 0.0047**, pseudo-R² 0.0854 → 0.0904. So U
does carry residual information about *where D fails*.

**Gate C — fails, and this is the one that matters.** A misranked pair means a
positive was scored too low *or* a negative too high, so a correction needs U
to say which. Fitted separately:

| side | n | U coefficient | p |
|---|---:|---:|---:|
| positives (scored too low) | 910 | **+1.122** | <0.0001 |
| negatives (scored too high) | 872 | **+1.573** | <0.0001 |

**Same sign on both.** U marks risk in both directions at once. A shift
`s' = s + δ(U)` moves high-risk positives and high-risk negatives the same way.

**Gate D — priced anyway, and it is worth nothing.** Over a grid of
`s' = s + c·U_rank`:

| | macro AUC | vs D |
|---|---:|---:|
| D as submitted | 0.8600 | — |
| best single global c (+0.025) | 0.8601 | **+0.0000** |
| oracle per-epitope c | 0.8614 | **+0.0014** |

And the per-epitope oracle directions **contradict each other** — GIL +0.100,
TTD +0.150, NLV 0.000, GLC −0.050, YLQ −0.150 — which is exactly what a
label-blind risk flag should do.

### Settled

**U is a real uncertainty statistic with no usable direction.** It survives the
controls the brief specified and is worth +0.0014 to an oracle that is allowed
to choose its sign per epitope. No correction is built; the brief's own warning
that E68 "must not become E49/E50/E52 under a new name" is honoured by stopping
here rather than by searching for a form that works.

The substantive finding is E67's: **D changed the level of the error and not
its shape**, 64% of the loss is still in pairs where neither member has
relatives, and GLCTLVAML's failure is diffuse rather than concentrated in
outliers — which argues against the "fix a handful of GLC cases" reading that
the macro-AUC sensitivity arithmetic invites.

---

## E70 — Validating TCR63's hinge, because it appeared only after a bug fix

The code review found E64's two-piece model never fitted — its knot started
outside the data, so the hinge term was zero everywhere. Bounding it reversed
the conclusion. **A result that appears only after a bug fix has to be
validated harder than one that was there all along.**
`scripts/64_hinge_validation.py`, six checks, the first of which can kill it.

**Check 1 — synthetic recovery, both halves.** x resampled from TCR63's own
values so the leverage structure matches:

| scenario | recovered slope change | held-out gain over linear |
|---|---:|---:|
| true hinge, change −60 | **−61.6** | +0.148 |
| true hinge, change −25 | −30.5 | **−0.021** |
| **pure line, change 0** | **+0.04** | **−0.006** |

**The fitter does not manufacture a knot**: on pure lines it recovers a change
of +0.04 and *loses* 0.006 out of sample. And the middle row is a useful
calibration — a moderate hinge is recovered as a parameter but does not pay off
out of sample, so a positive held-out gain requires a large one.

**Checks 2–5 on the real data.**

| | TCR59 | TCR63 |
|---|---|---|
| knot | 0.327 | 0.321 |
| slopes | +131.0 → +12.9 | **+248.9 → −65.5** |
| points each side | 98 / 94 | 73 / 119 |
| knot sitting at a bound | no | no |
| folds failing identifiability | 0/7 | 0/7 |
| knot after dropping 5/10/20 extreme-x points | 0.322 / 0.322 / 0.288 | 0.301 / 0.289 / 0.276 |
| bootstrap knot, 95% CI | [0.230, 0.577] | **[0.270, 0.447]** |
| gene holdout gain | **−0.020** | **+0.267** |
| origin holdout gain | −0.038 | **+0.310** |

**Check 6 — and it changes what the result means.** Saturation would put the
post-knot slope at zero; TCR63's goes *below* it:

| | slope after knot | 95% CI | P(slope<0) | Spearman above the knot |
|---|---:|---|---:|---:|
| TCR59 | +12.9 | [−55.0, +54.5] | 0.368 | +0.050 (n=94) |
| **TCR63** | **−65.5** | **[−129.0, −20.7]** | **1.000** | **−0.360** (n=119) |

**TCR63's relationship is not saturating, it reverses.** Above
`r2_r1_ratio ≈ 10^0.32 ≈ 2.1`, higher screen enrichment goes with *lower* flow
activation, and the model-free rank correlation in that region is **−0.360** on
119 peptides. That is why isotonic, which is constrained to be monotone,
reached only 0.424 while the hinge reaches 0.576–0.629 — the data are
non-monotone and a monotone model cannot represent them.

TCR59 flattens but does not reverse: its post-knot slope CI spans zero and the
rank correlation above the knot is +0.050, so a line captures its ordering just
as well and the hinge costs 0.020–0.038.

### Settled

**The hinge survives every check the brief specified**, and the finding is
sharper than first logged: TCR63 shows a **reproducible non-monotone reversal**
between screen enrichment and flow activation within the assayed range, stable
under gene holdout and the human↔virus split, not driven by extreme points, and
with a knot bootstrap CI of [0.270, 0.447].

It joins E61 from the other side. There, first-round abundance ran *negative*
against flow on TCR59 while the later-round gain ran +0.687. Here the
later-round gain itself reverses above a threshold on TCR63. Both say the same
thing: **beyond some point the screen stops tracking function**, which is what
to expect if the top of a selection is shaped by amplification, bottlenecking
or non-functional enrichment rather than by activation strength.

**What still cannot be claimed**, exactly as the brief states: this is not a
measured T-cell activation threshold, not a TCR–pMHC binding threshold and not
kinetic proofreading. SEEKER's two readouts are both functional-screen
quantities, and the 192 measured peptides were selected for being top-ranked.
Establishing a biological threshold needs peptides chosen *without* reference
to enrichment, across a dose gradient, with expression and presentation
controls.

---

## E69 — Does GLCTLVAML's public motif explain where D fails?

Groups fixed from the known GLC biology **before any error was examined**, as
the brief requires: TRAV5, TRBV20-1, and the single single-linkage family with
≥20 members (which E54 found is 27/27 TRBV20-1).
`scripts/65_glc_motif_decomposition.py`. The design is not vacuous — TRAV5
covers 96 of 192 cohort rows, TRBV20-1 covers 50, the family has 32.

**The groups, with D's behaviour inside each** (166 labelled rows):

| group | n | validated | mean D rank | AUC within group | share of misranking |
|---|---:|---:|---:|---:|---:|
| public motif family | 27 | **0.963** | 0.834 | **0.923** | 3.0% |
| TRAV5 only | 46 | 0.826 | 0.642 | 0.651 | 23.3% |
| **neither** | **80** | 0.512 | 0.315 | 0.627 | **61.2%** |
| TRAV5+TRBV20-1, no family | 8 | 0.500 | 0.569 | — | 8.7% |
| TRBV20-1 only | 5 | 0.200 | 0.308 | — | 3.8% |

**D is excellent inside the public motif (AUC 0.923) and mediocre outside it**
(0.627 and 0.651), and 61.2% of the misranking sits in the 80 rows that carry
no public motif at all. The whole-epitope AUC of 0.7568 is what those two
regimes average to.

**The residual-error test, conditioning on D's own score and nn3**, with
standard errors clustered by family:

| added term | coefficient | p | pseudo-R² |
|---|---:|---:|---|
| TRAV5 | +0.187 | 0.376 | 0.1675 → 0.1708 |
| TRBV20-1 | −0.044 | 0.896 | 0.1675 → 0.1677 |
| **in_big_family** | **−1.369** | **<0.0001** | 0.1675 → **0.2288** |

**TRAV5 and TRBV20-1 add nothing on their own** — E42's boundary holding
exactly as the brief predicted: V-gene identity is already inside D through the
V-domain term. Family membership does survive, and survives a provenance
control too (adding donor count leaves the coefficient at −1.349, p<0.0001,
with donors resolvable for 106/166 rows).

### The two checks that close it anyway

Surviving a residual regression is not enough. A correction must vary on the
private split and must say which way to move.

**(a) The family does not exist where it would be needed.** Rows in a family of
≥20, by epitope:

| | GIL | GLC | NLV | TTD | YLQ | **CIN** | **LLW** | **LTD** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| share | 30.6% | 16.7% | 0.0% | 0.0% | 0.0% | **0.0%** | **0.0%** | **0.0%** |

`in_big_family` is **identically false on all three private epitopes** — it
fails E3's `private_variance_gate` by construction, so a correction built on it
contributes exactly nothing where the ranking is decided. It is also false on
three of the five public epitopes, so it could not be validated across folds
either.

**(b) The direction is not estimable, and that explains the survival.**

| side | n | in family | coefficient | p |
|---|---:|---:|---:|---:|
| positives (scored too low) | 110 | 26 | −1.288 | <0.0001 |
| negatives (scored too high) | 56 | **1** | — | **not estimable** |

Only **1 of 56** GLC negatives is in the family. The family is 26/27 positive,
so conditioning on D's score does not stop it being very nearly the label
itself — which is why it survived the regression, and why nothing can be built
on it.

### Settled

**This line closes.** The pre-fixed hypothesis was partly right about the
biology: GLC does hold a reliable TRAV5/TRBV20-1 public motif, D ranks inside
it at 0.923, and the epitope's low overall AUC comes from the 80 rows outside
any motif. But the motif identity that explains it is already in the model
(E42), and the one term that adds something is label-like, directionless and
absent from every private epitope.

This also confirms E67's reading against the macro-AUC sensitivity arithmetic.
GLC's loss is diffuse — spread over 80 non-motif rows, not concentrated in
outliers — so there is no handful of cases to fix.

---

## E71 — Is there systematic sequence–function decoupling in GLC outside the public motif?

Two minimal analyses as asked, plus the control that makes the second one
readable. `scripts/66_glc_decoupling.py`, figure at
`reports/glc_score_distribution.html`.

**Why the second needs a control.** Removing the family removes 26 positives
and 1 negative, so the base rate moves and what remains is by construction the
harder half. An AUC drop is expected with no decoupling at all.

### 1. D's score distribution on GLCTLVAML

| | n | D rank median | IQR | in family |
|---|---:|---:|---|---:|
| validated | 110 | 0.606 | [0.383, 0.830] | 26 |
| not validated | 56 | 0.331 | [0.174, 0.476] | **1** |

24.5% of positives sit inside the negatives' interquartile range. The figure
shows why the whole-epitope AUC is 0.757: the validated distribution has a
sharp right peak at 0.8–0.95 that **is** the public-motif family, and removing
it leaves two heavily overlapping distributions.

### 2. AUC with the family removed, against three references

| epitope | removed | AUC all | **AUC without family** | 95% CI | nn3 alone on remainder |
|---|---:|---:|---:|---|---:|
| GILGFVFTL | 158 | 0.886 | **0.851** | [0.809, 0.890] | 0.812 |
| **GLCTLVAML** | 27 | 0.757 | **0.695** | **[0.603, 0.777]** | 0.614 |
| NLVPMVATV | 0 | 0.873 | 0.873 | [0.825, 0.914] | 0.845 |
| TTDPSFLGRY | 0 | 0.849 | 0.849 | [0.807, 0.887] | 0.860 |
| YLQPRTFLL | 0 | 0.935 | 0.935 | [0.908, 0.960] | 0.927 |

### The answer: **no, there is no large-scale decoupling**

GLC's remainder scores **0.695 with a 95% CI of [0.603, 0.777], which excludes
0.5**, and `nn3_dist` alone still reaches 0.614 there. The sequence
neighbourhood has not stopped predicting function outside the public motif —
it predicts it **weakly**. The control says part of the drop is simply the cost
of removing easy cases: GILGFVFTL loses 0.035 to the same operation and lands
at 0.851, while GLC loses 0.062 and lands at 0.695.

**And the weakening is a gradient, not a population.** Inside the remainder,
by `nn3` tercile:

| tercile | n | validated | AUC |
|---|---:|---:|---:|
| near | 46 | 0.717 | **0.774** |
| mid | 47 | 0.617 | 0.665 |
| far | 46 | 0.478 | **0.606** |

Validation rate and D's accuracy decline together and continuously, and even
the far tercile stays above chance. There is no sub-population where the
relationship breaks; there is a smooth loss of signal with distance.

### One specific, suggestive finding — and it is about calibration, not biology

GLC is the only epitope where restricting to **within-TRBV pairs improves** D
on the remainder:

| epitope | overall | within TRBV | delta | 95% CI | pairs kept |
|---|---:|---:|---:|---|---:|
| **GLCTLVAML** | 0.695 | **0.793** | **+0.098** | **[+0.005, +0.198]** | 11.8% |
| NLVPMVATV | 0.873 | 0.895 | +0.023 | [−0.046, +0.080] | 3.4% |
| TTDPSFLGRY | 0.849 | 0.850 | +0.001 | [−0.052, +0.052] | 9.0% |
| YLQPRTFLL | 0.935 | 0.915 | −0.020 | [−0.057, +0.012] | 6.3% |
| GILGFVFTL | 0.851 | 0.700 | **−0.152** | [−0.211, −0.093] | 34.4% |

If this holds, it says D **ranks correctly inside each TRBV group on GLC's
remainder and misorders the groups against each other** — a between-group
calibration failure rather than a failure on any TCR population.

**It is reported as suggestive, not established.** The CI's lower bound is
+0.005, only 11.8% of pairs survive the conditioning, and five epitopes were
tested without a multiplicity correction. GILGFVFTL moves the opposite way at
−0.152, which is itself informative: on an epitope whose response is
TRBV19-dominated, between-group ordering is where most of D's information
lives.

### Settled

Outside its public motif GLC is **harder, not decoupled**. The remaining
question that this analysis can pose but not answer is whether the between-TRBV
ordering is systematically wrong there — and even if it is, E69 already showed
the route is closed for the competition: V-gene identity is inside D through
the V-domain term (E42), and any fix would need a label-free way to order V
genes that works on three private cohorts with no families at all.

**Addendum to E71 — is GLC's optimal top-k larger, because it is saturated?**
E20's mechanism for GLC is that *everything is dense*, so the natural reading
is that top-3 saturates there and a larger k would restore resolution. E17
killed the analogous idea for kernel bandwidth, but its direction argument does
not obviously carry: for top-k, dense→large and sparse→small *is* geometric.

Tested on the density term alone, k ∈ {1,2,3,5,8,12,20,30,50,80}:

| epitope | close-neighbour fraction | oracle k | gain over k=8 |
|---|---:|---:|---:|
| GLCTLVAML | 0.552 | **8** | **+0.000** |
| GILGFVFTL | 0.476 | 12 | +0.002 |
| YLQPRTFLL | 0.184 | 20 | +0.001 |
| NLVPMVATV | 0.183 | **1** | +0.022 |
| TTDPSFLGRY | 0.047 | 3 | +0.009 |

**The idea dies on GLC specifically: its oracle k is already 8, the value D
uses, and larger k buys exactly nothing.** The oracle gain (+0.0069 overall)
comes entirely from the *sparse* cohorts wanting *small* k — E17's mechanism
again, where the few genuinely close pairs must dominate. Spearman between the
close-neighbour fraction and the oracle k is +0.500 on five points with
YLQPRTFLL breaking it outright. Not a usable rule, and no headroom on GLC.

---

## E72 — Why is GLCTLVAML over-converged? The proximate cause, measured

E20 established *that* GLC's density fails because "everything is dense". This
asks *why* its cohort is dense, against four competing explanations: a narrow
recognition solution space, a database artifact of deeply-sampled donors,
generation bias, or cohort composition. Six of the eight epitopes are
HLA-A\*02:01, so **HLA is held fixed by the design** and any explanation must
be peptide-specific.

**The shape of the convergence.** GLC is not one big motif:

| epitope | n | has a close neighbour | clusters | largest | in 2–4 member clusters |
|---|---:|---:|---:|---:|---:|
| **GLCTLVAML** | 192 | **55.2%** | 22 | 32 | **41.5%** |
| GILGFVFTL | 670 | 47.6% | 39 | 71 | 21.6% |
| NLVPMVATV | 366 | 18.3% | 21 | 17 | 56.7% |
| YLQPRTFLL | 419 | 18.4% | 19 | 12 | 51.9% |

**V-gene restriction is not the explanation.** Rarefied to n=192, TRBV Simpson
diversity is **0.390 on GILGFVFTL** against GLC's 0.853 and 0.886–0.955 for the
rest. GIL is far more V-restricted than GLC, yet density works there.

**Generation bias is not the explanation either**, and this was already
measured: E33's R08 found GIL's clustering-to-validation gap collapses from
0.400 to 0.067 once Pgen and length are matched, while **GLC's survives, 0.281
→ 0.375**. GLC's clusters are genuinely selected, not easy-to-generate.

**Donor concentration alone is not sufficient.** GLC draws on 20 resolvable
donors with the top three contributing 47%; NLVPMVATV is *more* concentrated
(43 donors, top three 69%) and shows none of the pathology.

### What the data does say

| epitope | donor resolution | close pairs | cross-donor | **same-donor** |
|---|---:|---:|---:|---:|
| GILGFVFTL | 86% | 1304 | 72% | 28% |
| NLVPMVATV | 89% | 53 | 68% | 32% |
| YLQPRTFLL | 37% | 13 | 69% | 31% |
| **GLCTLVAML** | 61% | 145 | **53%** | **47%** |

And the keystone, restricted to donors contributing ≥5 TCRs — **how convergent
is one person's response to this epitope?**

| epitope | qualifying donors | within-donor pairs ≤0.12 |
|---|---:|---:|
| **GLCTLVAML** | 9 | **8.7%** |
| GILGFVFTL | 74 | 4.2% |
| YLQPRTFLL | 11 | 0.6% |
| NLVPMVATV | 12 | 0.1% |
| TTDPSFLGRY | 8 | 0.0% |
| CINGVCWTV | 4 | 0.0% |

**Inside a single donor, GLC's TCRs sit within the clustering radius of each
other at twice GILGFVFTL's rate and 15–90× everything else.**

### The mechanism, as a chain

1. One person's GLC response is **oligoclonal** — 8.7% of within-donor pairs
   are inside the radius, the highest of all eight epitopes.
2. GLC's cohort draws on only **20 donors**, so those within-donor pairs are a
   large share of all close pairs: **47% same-donor**, against GIL's 28%.
3. E57 measured what that costs: on GLC, conditional on *having* a neighbour,
   a cross-donor neighbour predicts validation at **94.0%** and a same-donor-only
   neighbour at **60.0%** (odds ratio 10.4, p 1.5e-04).
4. So GLC's density is high but **spent on the uninformative kind of
   neighbour** — within-donor clonal expansion carries no convergence evidence.
5. Hence density collapses exactly where E20 found it: 0.551 on GLC's
   singletons, 0.677 for nn3 alone, the lowest of the eight.

**This replaces "everything is dense" with "the wrong things are dense"**,
which is both more specific and testable.

### The biology, separated from the measurement

What is measured above is the *proximate* cause. The *ultimate* cause — why one
person's GLC response is oligoclonal — this dataset cannot settle, and the
obvious candidate fails:

GLCTLVAML is the immunodominant HLA-A\*02:01 epitope of EBV's lytic BMLF1
protein, and EBV persists for life with periodic reactivation, so decades of
repeated restimulation focusing the response is the natural story. **But CMV's
NLVPMVATV is also a lifelong herpesvirus epitope with famously large responses,
and it shows 0.1% within-donor close pairs — essentially none.** "Persistent
herpesvirus" does not separate them. Something more specific to the GLC
response is responsible, and distinguishing the candidates (episodic lytic
reactivation versus continuous low-level stimulation; a narrower peptide
contact surface; founder effects in the acute infectious-mononucleosis
response) needs per-donor clone sizes and longitudinal sampling, neither of
which VDJdb carries here.

**Caveat on the keystone number.** GLC qualifies only 9 donors at 61% donor
resolution, so 8.7% rests on thin data. The contrast with NLVPMVATV — 12
donors at 89% resolution and 0.1% — is large enough that noise is an unlikely
explanation, but the number itself should not be quoted precisely.

---

## E73 — Should GLC be excluded from validation because the private side has no oligoclonality?

The proposal, and it is the *opposite* of what E20 forbade: E20 said do not tune
*for* GLC; this asks whether to drop GLC *from the validation set*, on the
grounds that its pathology (E72: within-donor oligoclonality inflating density)
does not exist on the private split.

**The premise holds, and extends further than proposed.** Close-neighbour
fraction, which is label-free and does not need donor information (LLW and LTD
have 0% donor coverage, E56):

| epitope | split | close-neighbour fraction |
|---|---|---:|
| **GLCTLVAML** | public | **55.2%** |
| **GILGFVFTL** | public | **47.6%** |
| YLQPRTFLL | public | 18.4% |
| NLVPMVATV | public | 18.3% |
| LTDEMIAQY | private | 14.6% |
| TTDPSFLGRY | public | 4.7% |
| LLWNGPMAV | private | 4.4% |
| CINGVCWTV | private | 2.5% |

The private range is 2.5–14.6%, matched by YLQ, NLV and TTD. **Both GLC and GIL
sit outside it** — so the proposal, taken seriously, excludes two epitopes.

### The direct test, and it refutes the proposal

E47's published tcrdist3 bands give the answer the proposal is trying to
approximate: the private macro is **0.767**. So ask which public subset
estimates it best.

| subset used | tcrdist3 estimate | error |
|---|---:|---:|
| all five public | 0.825 | +0.058 |
| the "private-like" three (YLQ, NLV, TTD) | 0.854 | **+0.088** |
| **the two excluded (GLC, GIL)** | 0.780 | **+0.014** |

**Excluding GLC makes the private estimate worse, not better.** tcrdist3's
measured per-epitope values explain why: GIL 0.871, **GLC 0.690**, NLV 0.809,
TTD 0.826, YLQ 0.928, against published private bands of 0.815 / 0.815 /
**0.67**. **GLC is the only public epitope in the weak band, and therefore the
only public proxy for LTDEMIAQY.** Dropping it removes the sole representative
of one third of the private set.

### And it would not have changed anything anyway

| configuration | all five | private-like three | rank |
|---|---:|---:|---|
| D | 0.8600 | 0.8857 | 1 → 1 |
| F | 0.8557 | 0.8838 | 2 → 2 |
| frozen | 0.8498 | 0.8795 | 3 → 3 |

The ordering is identical. The optimal structure weight moves **0.20 → 0.15**
under the restriction — *less* structure, not more, and both below the 0.30 in
use, whose justification was never a public optimum but E20's dominance
argument about sparse cohorts.

### Settled, and the general lesson

**Cohort geometry and difficulty are different axes, and they disagree on GLC.**

| | by density geometry | by published difficulty |
|---|---|---|
| is GLC like the private epitopes? | **no** (55.2% against 2.5–14.6%) | **it is the closest** (0.690 against LTD's 0.67 band) |

This is E21's result again — "sparsity does not determine performance;
intrinsic separability does" — and E49's, where GLC broke the `trbv_simpson`
ordering. Of the two axes, the published per-epitope AUROC is the more direct
evidence, because it *is* the quantity being predicted, while density is a
proxy for it.

The idea was sound and worth testing; the test says keep GLC. Nothing changes.

---

## E74 — The margin over tcrdist3 is not uniform, and the better model is validated

E47 anchored the private estimate on the paper's published tcrdist3 bands and
carried **our margin over tcrdist3 uniformly** across the three private
epitopes. E47 named that as the single remaining assumption. E73 made it
testable by establishing that the published per-epitope AUROC is the right axis
for comparing epitopes. This tests the assumption.

**It is wrong, and systematically.** Configuration D's margin against
tcrdist3's measured per-epitope AUC:

| epitope | tcrdist3 | D's margin | remaining headroom | fraction captured |
|---|---:|---:|---:|---:|
| GLCTLVAML | 0.690 | **+0.067** | 0.310 | 0.215 |
| NLVPMVATV | 0.809 | +0.064 | 0.191 | 0.335 |
| TTDPSFLGRY | 0.826 | +0.023 | 0.174 | 0.133 |
| GILGFVFTL | 0.871 | +0.015 | 0.129 | 0.120 |
| YLQPRTFLL | 0.928 | **+0.007** | 0.072 | 0.099 |

Spearman between tcrdist3's level and our margin is **−1.000, exact p 0.017**
(exact over all 120 permutations, per E49's lesson).

**Part of that is mechanical and must be removed.** AUC is bounded at 1, so a
baseline of 0.928 leaves only 0.072 to win. On the headroom-fraction scale the
relationship weakens to Spearman −0.900, **exact p 0.083** — so the ceiling
explains much of the slope but not all of it.

### The leave-one-epitope-out test, which decides between the models

Each model is fitted on four public epitopes and asked to predict the fifth's
margin — the discipline E40 demanded of E30's projection, which failed it at
MAE 0.083.

| margin model | held-out MAE |
|---|---:|
| uniform (what E47 assumed) | 0.0301 |
| linear in tcrdist3's level | 0.0161 |
| **constant fraction of remaining headroom** | **0.0156** |

**Both difficulty-aware models roughly halve the held-out error**, and the
headroom model — which is also the one with a mechanical justification — wins.
At MAE 0.0156 it is five times more accurate than E30's withdrawn projection.

### Revised private estimates

| configuration | public | held-out MAE | uniform (logged) | **headroom model** | band interval |
|---|---:|---:|---:|---:|---|
| frozen | 0.8498 | 0.0185 | 0.792 | **0.795** | [0.783, 0.807] |
| **D** | 0.8600 | 0.0156 | 0.802 | **0.809** | [0.798, 0.820] |
| F | 0.8557 | 0.0130 | 0.798 | **0.804** | [0.793, 0.815] |

D per-epitope under the headroom model: **CIN 0.848, LLW 0.848, LTD 0.730.**

The intervals come from the published bands' own width (0.80–0.83 and
0.66–0.68), not from sampling error, so they describe what the paper states
rather than what we measured.

### What changes and what does not

The estimates rise by +0.003 to +0.007 and the **ordering is unchanged** —
D > F > frozen under every margin model. So nothing about the submission
follows from this.

What does follow is that **E47's last stated assumption has now been tested
rather than assumed**, and replaced with a model that is validated on held-out
epitopes. The per-epitope breakdown is also no longer an assertion: LTDEMIAQY
is projected at 0.730 against CIN and LLW at 0.848, which is where the weak
band and the headroom argument both put it.

**The honest limits.** Five epitopes, one held-out fold each. The headroom
fraction itself ranges 0.099–0.335 across them, a 3.4× spread, so "constant
fraction" is true only on average and MAE 0.0156 is the price of that. And
LTDEMIAQY's published band (0.66–0.68) sits just below the public range's
minimum (GLC at 0.690), so its projection is a short extrapolation rather than
an interpolation.

---

## E75 — Structurally, how does GLCTLVAML differ? It does not.

Every GLC analysis so far has been distributional. This asks the structural
question from the 2,636 Boltz-2 complexes, with a control already built into
the data: **six of eight epitopes are HLA-A\*02:01 9-mers**, so HLA and peptide
length are fixed and anything left is peptide-specific.
`scripts/67_structural_comparison.py`.

The hypothesis from E72 was that GLC's over-convergence reflects a narrow set
of docking solutions, which should show as a tighter ensemble, a narrower
contact footprint and less variable geometry.

**1. The docking geometry is unexceptional.** Across six spread measures
(crossing angle, TCR height, α and β contact counts, CDR3 shares), GLCTLVAML's
mean rank is **4.3 of 6** where 1 is tightest. Its crossing-angle spread is
15.13°, second *highest* after YLQPRTFLL's 15.75°.

**2. The contact footprint is nearly identical across all six**, median
α+β contacts per peptide position:

| epitope | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GLCTLVAML | 0 | 0 | 7 | **34** | 1 | 0 | 27 | 11 | 26 |
| GILGFVFTL | 0 | 0 | 8 | **35** | 1 | 0 | 16 | 10 | 26 |
| NLVPMVATV | 2 | 0 | 8 | **34** | 1 | 0 | 22 | 12 | 26 |
| YLQPRTFLL | 3 | 1 | 7 | **32** | 2 | 0 | 30 | 16 | 26 |
| LLWNGPMAV | 9 | 0 | 7 | **35** | 0 | 0 | 18 | 7 | 27 |
| CINGVCWTV | 0 | 1 | 6 | 33 | 1 | 0 | **35** | 12 | 26 |

P4 dominates everywhere but CINGVCWTV, P9 sits at 26–27 for all six, P3 at
7–8 for all six. The only position that genuinely varies is **P7 (16–35)**.

**3. "Ensemble tightness" is a Boltz artifact, not biology.** The pairwise
fingerprint distance within an epitope correlates with cohort size at
**ρ 0.000** and with interface ipTM at **ρ −0.810** — an epitope Boltz folds
confidently produces similar structures whatever the biology. After regressing
tightness on ipTM, **GLCTLVAML's residual is +0.0000, rank 5 of 8: exactly
average.** LLWNGPMAV looked tightest (0.0783) purely because it has the highest
ipTM (0.9472).

**4. The one apparent relationship does not survive.** Footprint width
(effective number of contacted positions) against close-neighbour fraction
gives ρ −0.486 — the two narrowest footprints, GLCTLVAML (4.48) and GILGFVFTL
(4.49), are exactly the two most sequence-convergent epitopes. But the **exact
p over all 720 permutations is 0.356**, and removing those two flips the
remaining four to **+0.800**. Two points out of six carry it and the rest
reverse.

### Settled, and it is partly a statement about Boltz

**GLC is structurally indistinguishable from the other A\*02:01 9-mers in
anything these structures resolve.** Its sequence convergence (55.2% with a
close neighbour, the highest of eight) has no structural counterpart.

That joins E36 from a second direction. There, folding each TCR against
anchor-matched decoy peptides gave a margin of **+0.0061** in ipTM — Boltz
assigns a TCR's true target and an unrelated A2 nonamer essentially the same
confidence. Here the contact profiles say the same thing in a different
currency: **Boltz builds nearly the same interface for every A\*02:01 9-mer**,
P4-dominant with a P9 secondary contact and a variable P7. The negative is
therefore at least as much about what the folding model can resolve as about
what distinguishes the epitopes, and separating those two would need
experimental structures across these epitopes, not more folds.

---

## E76 — Is GLC's over-convergence a sampling artifact that more donors would erase?

The proposal: if E72 is right that GLC's density comes from within-donor
oligoclonality and a 20-donor cohort, then the epitope itself is unremarkable
and a larger sample would bring its TCR distribution into line with the others.

**The literal version is backwards.** GLC is the *smallest* cohort (192) with
the *highest* close-neighbour fraction (55.2%), while GILGFVFTL has 670 rows at
47.6%. Adding TCRs raises density; it cannot lower it. The testable version is
the mechanism, not the count: **sample more donors, not more TCRs** — which is
measurable now by counting only cross-donor neighbours, since that removes
within-donor expansion without shrinking the cohort.

### Removing within-donor neighbours

| epitope | resolvable n | all neighbours | **cross-donor only** | drop |
|---|---:|---:|---:|---:|
| **GLCTLVAML** | 118 | 56.8% | **28.8%** | **−28.0 pp** |
| GILGFVFTL | 576 | 47.2% | 42.0% | −5.2 |
| NLVPMVATV | 326 | 14.7% | 8.6% | −6.1 |
| YLQPRTFLL | 155 | 11.6% | 7.7% | −3.9 |
| CINGVCWTV | 188 | 2.1% | 0.0% | −2.1 |
| TTDPSFLGRY | 141 | 1.4% | 1.4% | −0.0 |

**GLC loses more than half its density, a 28-point drop against 3.9–6.1 for
everything else — five times larger.** The proposal is substantially right.

**But not entirely.** What survives is 28.8%, still second only to GILGFVFTL's
42.0% and three to four times NLVPMVATV and YLQPRTFLL. More donors would move
GLC from first place to second, not into the pack.

### Is the epitope itself low quality? No — the composition is

Classifying each labelled TCR by what kind of neighbour it has:

| epitope | cross-donor neighbour | only same-donor | singleton |
|---|---|---|---|
| **GLCTLVAML** | n=**32**, 94% validated | n=**33**, 55% | n=76, 51% |
| GILGFVFTL | n=**195**, 94% validated | n=**34**, 50% | n=317, 48% |
| NLVPMVATV | n=28, 75% | n=23, 74% | n=252, 15% |

**GLC and GILGFVFTL behave identically group for group** — a cross-donor
neighbour means 94% validated on both, and a same-donor-only neighbour means
50–55%, which is the base rate. What differs is the **mix**: GLC runs 32 : 33,
a one-to-one ratio of informative to uninformative clustered TCRs, while GIL
runs 195 : 34, six to one.

So the epitope's biology is not mediocre and its labels are not noisier. **Its
cohort's donor composition is what is unusual**, and it dilutes the one feature
the model depends on. This also disposes of a competing reading: E4 found
VDJdb's own confidence score anti-correlates with validation overall (pooled
0.455) but reaches 0.581 on GLCTLVAML — the only epitope where it is above
chance, which argues against "GLC's metadata is worse".

### The prediction this makes

If the mix is the problem, then a GLCTLVAML cohort assembled from **more donors
with fewer TCRs each** should have a density feature that behaves like
GILGFVFTL's, and its AUC should rise from 0.757 toward the 0.85–0.89 the other
epitopes reach. That is a prediction about data collection, not about method,
and it cannot be tested on this cohort — capping TCRs per donor here leaves
too few rows to measure density on.

**Caveat.** GLCTLVAML's donor resolution is 61%, so 51 of its 192 rows are
unclassified and every number above rests on the other 141. GIL and NLV are at
86% and 89%, so the comparison is between a well-resolved and a
partly-resolved cohort.

---

## E77 — De-duplicating the top-8 neighbourhood by family. The hypothesis is backwards.

Proposed externally (as "E73" in a collaborator's plan; renumbered here because
E73–E76 are taken). The idea: a TCR's eight nearest neighbours may represent
only two or three independent recognition modes, so capping each family's
contribution should sharpen the ranking. It is motivated by our own E57, E72
and E76, and it is genuinely new — E11 tested graph topology, E17 adaptive
bandwidth, E50 class-conditioned weights, but nobody had de-duplicated the
top-k set. It also uses *family* rather than donor, so unlike E57's machinery
it is computable on the private split.

**The premise holds on two epitopes and fails on six.** Effective number of
distinct families among the top 8:

| epitope | split | effective families | top-8 all one family |
|---|---|---:|---:|
| GILGFVFTL | public | **4.34** | 28% |
| GLCTLVAML | public | **4.06** | 28% |
| YLQPRTFLL | public | 6.50 | 3% |
| NLVPMVATV | public | 6.71 | 5% |
| **LTDEMIAQY** | private | **7.27** | 0% |
| **LLWNGPMAV** | private | **7.79** | 0% |
| **CINGVCWTV** | private | **7.91** | 0% |

**On all three private epitopes the top-8 already spans 7.3–7.9 distinct
families**, so there is essentially nothing to de-duplicate where the ranking
is decided. That alone closes it for the competition.

**And de-duplication makes the density worse**, macro 0.8414 → **0.8381**
(−0.0033). It helps only GLCTLVAML (+0.0170) and costs GILGFVFTL (−0.0122) and
NLVPMVATV (−0.0186).

**The proposal's own stopping conditions both fire.** The within-epitope
Spearman between the standard and de-duplicated score is **+0.9667** ("stop if
the new feature is nearly identical to D's ranking"), and the only epitope that
gains is the public-motif one ("stop if it only moves public motif scores").

**The hypothesis is backwards, and the mechanism is explicit.** Scoring the
*effective family count* on its own gives **0.2957** — far below chance, so
more independent families means *less* likely to validate. A top-8 that
collapses into one family is not redundancy to be discounted; it is the
strongest convergence evidence available, and E54 measured that large families
are 96.2% validated. De-duplication removes exactly that.

**A prediction for the companion proposal** (background-corrected neighbourhood
enrichment, matched on V/J, length and Pgen): the same failure, for a reason
already measured. E8 ran the core version — dividing cohort density by density
against the other epitopes' TCRs — and got 0.8396 → 0.8262, with **background
density alone scoring 0.60 in the positive direction**. Its settled mechanism
was that "crowded regions of TCR space are where genuinely public clonotypes
live, so normalising subtracts signal rather than noise". E77 now says the same
thing from a third angle: local redundancy *is* the signal. A finer matched
background is a refinement of a design whose sign is already known.

---

## E78 — Matched-background neighbourhood enrichment. E8's mechanism survives the refinement.

Proposed externally (as "E74"; renumbered). The question: a TCR has many
references nearby, but **given its V gene, CDR3 length and generation
background, how unusual is that?** E8 ran the coarse version and lost; this is
the fine-grained one it did not run, with the background matched per query on
**TRBV stem × CDR3β length bin × paired-Pgen tercile**, relaxing through
`trbv+len`, `len+pgen`, `len only` and finally the whole foreign pool.
`scripts/68_matched_background_enrichment.py`. Median background pool 238.

**The prediction was recorded before running**, because a prediction is not a
measurement — the same discipline E36, E50 and E58 applied to their own
controls.

| feature | macro | GIL | GLC | NLV | TTD | YLQ |
|---|---:|---:|---:|---:|---:|---:|
| own-cohort top-8 (raw) | **0.8414** | 0.880 | 0.714 | 0.833 | 0.851 | 0.930 |
| **matched background alone** | **0.5881** | 0.596 | 0.694 | 0.574 | 0.601 | 0.476 |
| enrichment = own − background | 0.8154 | 0.883 | 0.675 | 0.812 | 0.793 | 0.914 |

**Dividing the background out costs 0.026.** And into configuration D: the
enrichment correlates with D at within-epitope ρ +0.826, an **oracle allowed to
choose a single global weight picks w = 0.00** for a gain of +0.0000, and
nested LOEO gives −0.0014.

**E8's mechanism survives the refinement intact.** The matched background
scores **0.5881 on its own** — against E8's 0.60 for the coarse version — so it
is positively informative, and subtracting it removes signal however finely it
is matched. Crowded regions of TCR space are where genuinely public clonotypes
live, and conditioning on germline and generation class does not change that.

This is the third independent route to the same conclusion. E8 measured it with
a coarse background, E77 measured it by de-duplicating the top-8 (−0.0033, with
effective family count scoring 0.2957 on its own), and E78 measures it with a
matched background. **Local crowding is the signal, not a nuisance to be
normalised away.**

---

## E79 — Weighting a neighbour's contribution by its annotation quality

Proposed externally (as "E75"). Genuinely different from what is closed: E4
scored `vdjdb_score` as a feature **of the query** (pooled 0.455, below chance)
and E33's S05 did the same for provenance (0.534–0.557, constant within private
epitopes). This weights **the neighbour**, using only VDJdb metadata and no
validation labels.

**E3's gate fires before the experiment**, on the fields the weighting needs:

| epitope | split | distinct `n_refs` | distinct `score` | distinct `n_methods` |
|---|---|---:|---:|---:|
| GILGFVFTL | public | 11 | 3 | 7 |
| GLCTLVAML | public | 7 | 4 | 4 |
| NLVPMVATV | public | 5 | 4 | 4 |
| CINGVCWTV | private | 2 | 2 | **1** |
| LLWNGPMAV | private | 2 | **1** | 2 |
| **LTDEMIAQY** | private | **1** | **1** | **1** |
| TTDPSFLGRY | public | **1** | **1** | **1** |

**On LTDEMIAQY all three weights are constant**, so the weighting is a no-op
there, and TTDPSFLGRY is the same, leaving four usable validation folds.

Run anyway, weights floored at 1 so a `score = 0` neighbour is not erased:

| weighting | macro | vs unweighted |
|---|---:|---:|
| unweighted | 0.8414 | — |
| by `n_refs` | **0.8421** | +0.0007 |
| by `n_methods` | 0.8409 | −0.0005 |
| by `vdjdb.score` | 0.8388 | −0.0026 |

Into configuration D, an oracle allowed a single global weight picks **w = 0.00**
for **+0.0000**.

**A gate limitation worth recording.** `private_variance_gate` *passes* every
one of these (243/229/130 distinct values on the private epitopes) even though
the weight itself is constant on LTDEMIAQY. The weighted density inherits its
variance from the distances, not from the weighting, so **the gate can pass a
feature whose novel component is constant**. Checking the variance of the new
*ingredient*, not just of the resulting score, is the stricter test and it is
what caught this.

---

## E80 — Cross-epitope reference expansion constrained by peptide similarity

Proposed externally (as "E76"). Genuinely untried: every pool experiment here
(E12, E27, E33, E35) used **same-epitope** TCRs. This uses TCRs annotated to
*different but similar* peptides, on the reasoning that a receptor specific to
a near neighbour of the target may be a legitimate reference.

**Availability first.** Paired-αβ VDJdb TCRs on peptides of the same HLA class
and length within Hamming 4 of each target:

| target | split | nearest similar peptide | records available |
|---|---|---|---:|
| NLVPMVATV | public | ALVPMVATV (d=1) | 399 |
| **CINGVCWTV** | private | CVNGSCFTV (d=3) | **14** |
| **LLWNGPMAV** | private | LLEWLAMAV (d=4) | **3** |
| **LTDEMIAQY** | private | AMDEFIEQY (d=4) | **3** |
| GILGFVFTL | public | GILEFVFTL (d=1) | 2 |
| GLCTLVAML | public | SLYNLVATL (d=4) | 1 |
| YLQPRTFLL / TTDPSFLGRY | public | — | 0 |

**The three private epitopes can gain 20 records between them**, against 602
rows — a multiplier of 1.03, where E27 measured that even 1.834 yields +0.0013
realisable.

**NLVPMVATV is the only case where the idea can be evaluated at all**, and
testing it exposes why the counts are misleading. Its 126 "similar peptides"
include ALVPMVATV, CLVPMVATV, DLVPMVATV and ELVPMVATV all at Hamming 1 —
**a positional scanning mutagenesis panel of NLVPMVATV itself**, with the same
receptors tested against every variant. After deduplicating by paired clonotype
and removing TCRs already in the cohort, the 399 records collapse to **13 new
clonotypes**, worth **+0.0011** (0.8330 → 0.8341).

**Settled, with a data-structure finding worth keeping.** Hamming-based peptide
similarity in VDJdb retrieves **mutagenesis variants of the same epitope**, not
biologically distinct similar epitopes, so the TCRs it returns are largely the
ones already present. Any future cross-epitope expansion has to filter the
peptide list by provenance — a variant panel from one study is not a second
epitope — before counting what it adds.

---

## E81 — Twelve representations, including three trained on corrupted data

Proposed externally (as "E77": a new paired-αβ representation). No
`transformers` or `esm` locally, but SCEPTR ships **fourteen variants** and E6
tested only seven. The ten untested ones include three published **ablation
controls**, which answer a question this project has circled for weeks: how
much of the representation's value is learned antigen biology?

A gene-name fix was needed first — every variant rejects `TRAV14_DV4*01`; the
working path replaces `_` with `/` and drops genes outside the model vocabulary,
exactly as `scripts/33`'s embedder does.

**Top-8 density, per epitope:**

| variant | GIL | GLC | NLV | TTD | YLQ | macro |
|---|---:|---:|---:|---:|---:|---:|
| **blosum** | 0.882 | 0.711 | 0.836 | **0.875** | 0.925 | **0.8460** |
| left_aligned | 0.883 | 0.707 | 0.841 | 0.858 | 0.922 | 0.8422 |
| default (used by D) | 0.880 | 0.714 | 0.833 | 0.851 | 0.930 | 0.8414 |
| **dropout_noise_only** | 0.870 | 0.709 | 0.820 | 0.861 | 0.931 | 0.8381 |
| **shuffled_data** | 0.878 | 0.705 | 0.827 | 0.853 | 0.925 | **0.8377** |
| tiny | 0.884 | 0.741 | 0.815 | 0.803 | 0.933 | 0.8352 |
| small | 0.884 | 0.720 | 0.821 | 0.832 | 0.919 | 0.8349 |
| average_pooling | 0.877 | 0.709 | 0.830 | 0.849 | 0.909 | 0.8347 |
| **synthetic_data** | 0.879 | 0.700 | 0.830 | 0.840 | 0.920 | **0.8337** |
| mlm_only | 0.882 | **0.750** | 0.819 | 0.804 | 0.911 | 0.8333 |
| large | 0.876 | 0.687 | 0.837 | 0.837 | 0.913 | 0.8300 |
| finetuned | 0.872 | 0.720 | 0.789 | 0.795 | 0.929 | 0.8208 |

### The controls are the result

- **`synthetic_data`** was trained on a size-matched set of **OLGA-generated
  sequences** — a probabilistic VDJ recombination model with **no antigen
  selection whatever**. It scores **0.8337** against the real model's 0.8414, a
  gap of **0.0077**.
- **`shuffled_data`** was trained on the same real repertoire with **α-β chain
  pairing randomised**. It scores **0.8377**, a gap of **0.0037**.

**A representation that has never seen a correctly paired, antigen-selected TCR
reaches within 0.008 of one that has.** The entire span across twelve
architecturally and data-wise very different models is 0.8208–0.8460, 0.025
wide. That is the sharpest version yet of what E6 found ("the metric space does
not matter"), and it joins E7, E33's R08, E77 and E78: what these
representations encode that matters here is recombination-space geometry, not
learned antigen specificity. `shuffled_data`'s near-parity also agrees with E45
and E63, where αβ pairing carried nothing beyond per-chain typicality.

### And `blosum` is not a justified swap, for a measured reason

It is the best single representation ever tested here, +0.0046 over the default
D uses. But choosing it *is* a selection. Nested over the twelve, picking the
best on four epitopes and applying it to the fifth:

| | macro |
|---|---:|
| fixed `blosum`, chosen post hoc | 0.8460 |
| fixed `default`, what D uses | 0.8414 |
| **nested variant selection** | **0.8316** |

**Selection costs −0.0144** — worse than either fixed choice. The mechanism is
visible in the folds: `blosum` wins four of five, but the TTDPSFLGRY fold picks
`tiny`, which scores 0.803 there against `blosum`'s 0.875, and that one pick
erases the gain. This is E15's −0.0048 again at three times the magnitude,
because the variants differ more than the configurations E15 swept.

Blending `blosum` into D rather than swapping gives +0.0006 at an oracle weight
of 0.10–0.15 and **−0.0035 under nested LOEO**; it correlates with D at
within-epitope ρ +0.905, so it is measuring what D already measures.

**The entry is unchanged.**


**Amended by the 2026-10-27 review.** The table above puts all twelve variants
at a **common top-8 reduction**, so the comparison between them is exact. E6's
published numbers used each representation's *own best* reduction, so the two
sets are not comparable row by row: `default` reads **0.8414** at top-8 here
against E6's **0.8396** at kernel σ=0.1. Nothing in the conclusion depends on
it — the controls, the span and the selection cost are all computed within this
table — but the earlier text invited the wrong comparison.

---

## E82 — The two-task control: is antigen specificity in the pretraining or in the retrieval?

A collaborator reframed E81 correctly and more sharply than "protein language
models have not learned immunology". The system decomposes into three layers —
a pretrained sequence representation, an **antigen-conditioned reference set**
of TCRs known to bind the target peptide, and a ranking of functional
revalidation against that reference. On that reading, an antigen-agnostic
embedding can still support good prediction because **the retrieval, not the
representation, carries the antigen information**:

> antigen-conditioned retrieval over an antigen-agnostic sequence representation

The proposed decisive control is two tasks with the same representations:
**A**, IMMREP26-style functional revalidation (done, E81), and **B**, independent
recognition generalisation measured experimentally rather than from VDJdb.

### Task B is runnable, in a reduced form

What each external source can and cannot supply turned out to matter:

| source | paired sequences | shared peptide panel |
|---|---|---|
| PRP | **CDR3β only** (GitHub `tcr_cdr3b.csv`) | yes — 38k–100k common peptides per pair |
| SEEKER PacBio | **yes**, 101 full αβ with V/J | **no** — 59,812 rows hold 59,812 distinct peptides, each assigned to one TCR, so every pair shares zero |
| SEEKER 4TCR | no sequences in the package | yes — 236 peptides across all four, but only 4 receptors |
| SEEKER B27 | — | 2 receptors |

**The pieces exist but never together.** What is available is PRP: 21 receptors
over a shared library, β sequence only. Every SCEPTR variant will embed from
CDR3β alone, so the comparison can be made β-only.

### Result: nothing works, including the real-pretrained model

156 comparable pairs, recognition distance = 1 − Spearman on ≥5,000 common
peptides within one anchor library:

| representation | pretrained on | ρ with recognition distance | p |
|---|---|---:|---:|
| default | real repertoire | **+0.015** | 0.850 |
| blosum | real, BLOSUM embedding | +0.009 | 0.914 |
| shuffled_data | real, αβ pairing randomised | +0.028 | 0.725 |
| **synthetic_data** | **OLGA synthetic only** | **+0.039** | 0.628 |
| cdr3_only | real, CDR3 only | +0.034 | 0.669 |
| b_sceptr | real, β chain only | +0.043 | 0.598 |
| Levenshtein (E63's baseline) | — | +0.015 | 0.854 |

**Every value is indistinguishable from zero**, and the OLGA-pretrained model
is numerically *above* the real-pretrained one.

### Which branch of the proposal this lands in — neither, exactly

The framework offered two outcomes: OLGA matching on A but losing on B (the
pretraining carries task-specific value), or matching on both (recombination
geometry already approximates the biophysics). **The actual outcome is a third:
on task B nothing works at all**, so there is no difference to attribute. The
conclusion is not "recombination structure suffices" but "**β sequence carries
no recognition-profile similarity in any representation tested**".

**The bound that matters.** This is β-only, because PRP has no α. E63 found the
same dataset's five identical-CDR3β pairs recognise no more similarly than
random pairs — and those pairs differ precisely in α. So α is doing work that β
cannot, and this experiment does not test whether a *paired* representation
would succeed.

**What would settle it**, stated precisely because the proposal deserves it:
**≥10 fully paired αβ TCRs profiled over a shared peptide panel.** Neither PRP
(β only) nor SEEKER (pairs without a shared panel, or a shared panel without
sequences) provides it, and that gap is worth more as a specification for the
next dataset than as a reason to weaken the claim.

---

## Code review of E77–E82 (2026-10-27)

**The first finding is a process failure of mine.** Five of the six experiments
were run as inline heredocs and left **no retained script** — E77, E79, E80,
E81 and E82. That is exactly the gap E60 criticised in a collaborator's ledger
("M70 and pool v3 have conclusions but no code, and cannot be re-run"). Fixed:
`scripts/69` through `scripts/73` now implement them, and **every one
reproduces its logged numbers exactly** —

| experiment | script | check |
|---|---|---|
| E77 family de-duplication | `69_family_dedup_topk.py` | 0.8414 → 0.8381, effective count 0.2957, ρ +0.9667 ✓ |
| E79 reference quality weighting | `70_reference_quality_weighting.py` | n_refs +0.0007, into D +0.0000 ✓ |
| E80 cross-epitope expansion | `71_cross_epitope_expansion.py` | 20 private records, 13 clonotypes, 0.8330 → 0.8341 ✓ |
| E81 representation ablations | `72_representation_ablations.py` | nested 0.8316, cost −0.0144 ✓ |
| E82 task B | `73_task_b_recognition.py` | all ρ ≤ +0.043, all p ≥ 0.598 ✓ |

**No computational bugs were found in the six experiments.** The bug classes
that bit E62–E64 were checked again and are absent here: no free-parameter
curve fitting, no isotonic direction assumption, no reliability ceiling
compared raw, and every permutation or resample uses the same unit as its
observed statistic.

**Two reporting defects, neither changing a conclusion.**

1. **E81 mixed reductions when citing E6.** All twelve variants there use a
   common top-8, so the comparison between them is exact, but E6's published
   figures used each representation's own best reduction. `default` is 0.8414
   at top-8 and 0.8396 at kernel σ=0.1. The entry is amended in place.
2. **E82's 156 pairs are not independent.** Twenty-one receptors carrying
   sixteen distinct CDR3β generate 156 pairs, so each receptor appears in about
   seventeen and the p-values are anticonservative. Every result is null at
   p ≥ 0.598, so the conclusion is unaffected; the limit is now in the script's
   docstring alongside the β-only bound.

**The standing rule this re-earns:** an experiment without a retained script is
not a result, whoever ran it. The five scripts were written after the fact and
checked against the logged numbers rather than trusted.

---

## E83 — The contact fingerprint as a density space, not a confidence score

Every structural feature tried before (E14, E19, E28, E36) is a **scalar
summarising one interface**, and E36 settled that those carry no peptide-
specific information — the anchor-matched decoy margin was +0.0061 against a
decoy-choice noise of 0.036. This is a different axis: the 60-dimensional
contact fingerprint from `scripts/24` is indexed **by peptide position**, and
within an epitope the peptide is identical, so two TCRs can be compared by
*where they touch the antigen*. Not "how confident is this interface" but
"**do TCRs that engage the peptide the same way validate the same way**".
`scripts/74_structural_neighbourhood.py`.

Both scalings are reported rather than one being chosen, since the 60 columns
mix contact counts, capped distances and CDR3 shares:

| feature | macro | GIL | GLC | NLV | TTD | YLQ |
|---|---:|---:|---:|---:|---:|---:|
| `pae_tcr_pep_min` (D's term) | **0.7565** | 0.779 | 0.736 | 0.763 | 0.629 | 0.875 |
| fingerprint, column-standardised | 0.7510 | 0.827 | 0.730 | 0.764 | 0.581 | 0.853 |
| fingerprint, raw cosine | 0.7428 | 0.817 | 0.736 | 0.723 | 0.553 | 0.884 |

**The structural neighbourhood is slightly worse than the structural scalar it
would replace** — but it is genuinely decorrelated, at within-epitope ρ **+0.444**
with `pae_tcr_pep_min` and **+0.535** with D as a whole, against blosum's +0.905
(E81) and the matched-background enrichment's +0.826 (E78). It also clears the
typicality gate (0.7510 against 0.6378, ρ +0.216).

**Adding it to D gives +0.0007 at an oracle weight and −0.0024 nested.** But
*combining the two structural views inside D's existing structure slot* is
different:

| D's structure term | macro | vs D |
|---|---:|---:|
| `pae` alone (D as submitted) | 0.8600 | — |
| fingerprint alone | 0.8557 | −0.0043 |
| **pae 0.7 + fingerprint 0.3** | **0.8617** | **+0.0017** |
| pae 0.5 + fingerprint 0.5 | 0.8616 | +0.0015 |
| **nested LOEO over the three forms** | **0.8616** | **+0.0015** |

**This is the first candidate in thirty-odd experiments whose nested LOEO is
positive**, which is worth recording even though +0.0015 is far below the +0.02
bar and the entry does not change. The reason it behaves differently from the
other near-misses is visible in the correlations: it is not re-encoding what D
already has. E28's rule was "decorrelation without margin over chance buys
nothing" — this has both and still buys only +0.0015, which sets the price of
the one structural axis that was still unexplored.

**What it does not do** is overturn E36. That result is about the
TCR-versus-peptide axis, and this is about TCR-versus-TCR similarity in a
structural coordinate; the fingerprint being mildly useful does not mean Boltz
resolves peptide specificity, and E75 measured that it builds nearly the same
interface for every A\*02:01 nonamer.

---

## E84 — What are the three submissions actually disagreeing about?

Three entries exist — frozen (public 0.8498), D (0.8600), F (0.8557) — and every
quantity separating them comes from the published bands rather than from the
private rows. So this asks something the private labels are not needed for:
*what are they disagreeing about?* `scripts/75_submission_disagreement.py`.

ROC AUC **is** a pairwise statistic, so the natural decomposition is pairwise.
Within an epitope, every (positive, negative) pair falls in one of four cells —
both right, neither right, only A, only B — and the first two cancel exactly.
The macro difference is therefore *identically* the per-epitope mean of
(only-A − only-B); the script asserts this to 1e-9, which is also a check that
all three files are aligned on the same rows.

| epitope | pairs | both | neither | only D | only frozen | ΔAUC |
|---|---:|---:|---:|---:|---:|---:|
| GILGFVFTL | 73914 | 0.863 | 0.108 | 0.023 | 0.006 | +0.0169 |
| GLCTLVAML | 6160 | 0.723 | 0.225 | 0.034 | 0.018 | +0.0154 |
| NLVPMVATV | 20566 | 0.846 | 0.108 | 0.027 | 0.019 | +0.0080 |
| TTDPSFLGRY | 29240 | 0.819 | 0.129 | 0.031 | 0.022 | +0.0092 |
| YLQPRTFLL | 35392 | 0.923 | 0.055 | 0.012 | 0.010 | +0.0015 |
| **MACRO** | | | | | | **+0.0102** |

**The whole difference between the entries is 5,745 discordant pairs out of
165,272 — 3.5%**, of which D wins 66.0%. Ninety-six point five per cent of the
ranking is common to all three. The disagreement is also not close on either
side: D's median winning gap is 0.0384 and frozen's 0.0445, so these are not
tie-breaks, they are genuine reorderings in the middle of the list.

### The pool question, and why the obvious answer to it is wrong

D holds a pool-expanded density term and E27 priced its public-to-private
transfer at +0.0112 → +0.0013. So: are D's winning pairs on TCRs that actually
acquired pool neighbours? Counting, per TCR, the share of its top-8 density
neighbours that are pool rows:

| pool share of top-8 | pairs | D wins | net pairs to D |
|---|---:|---:|---:|
| neither member (0) | 885 | 59.5% | +169 |
| 0 < max ≤ 0.25 | 1264 | 52.4% | +60 |
| 0.25 < max ≤ 0.5 | 404 | 67.8% | +144 |
| max > 0.5 | 3192 | 72.9% | **+1460** |

That reads as 91% of D's margin sitting on pool-affected pairs, which would be
bad news. **It is confounded, and the confound is visible as soon as the epitope
is held fixed:**

| epitope | pool-free pairs | D wins | pooled pairs | D wins |
|---|---:|---:|---:|---:|
| TTDPSFLGRY | 696 | 57.0% | 836 | 60.3% |
| YLQPRTFLL | 189 | 68.8% | 584 | 48.6% |

Within TTD the two strata are indistinguishable and within YLQ the sign
*reverses*. The pool-affected TCRs are simply the ones in the large cohorts
(GIL 3,005 pool rows, NLV 9,881) and those are also where D wins for unrelated
reasons. GIL, GLC and NLV contribute no pool-free pairs at all, so the table
above is mostly an epitope effect wearing a pool label.

**The decomposition that is a control rather than a correlation already
exists**, because F *is* D with the pool removed everywhere:

- `F − frozen` = **+0.0059**, 57% of D's margin, pool-free **by construction**
- `D − F` = **+0.0044**, 43% of D's margin, the pool's own contribution

and the pool term shows a clean dose–response against how much pool each cohort
has:

| epitope | pool rows | top-8 pool share | D − F |
|---|---:|---:|---:|
| NLVPMVATV | 9881 | 0.781 | +0.0058 |
| GILGFVFTL | 3005 | 0.702 | +0.0044 |
| GLCTLVAML | 459 | 0.556 | +0.0115 |
| YLQPRTFLL | 54 | 0.100 | +0.0003 |
| TTDPSFLGRY | 25 | 0.067 | −0.0002 |

Spearman +0.700. **Below ~60 pool rows the term does nothing at all** — TTD and
YLQ move by 0.0003 — and the three private epitopes hold **11, 191 and 5** pool
rows. So D's pool half should contribute approximately nothing on CINGVCWTV and
LTDEMIAQY and something on LLWNGPMAV alone, which is an independent measurement
of E27's discount rather than a restatement of it: D's realisable advantage over
F is about one third of +0.0044, and D's advantage over frozen on the private
side is carried by the pool-free +0.0059.

**The ordering is unchanged** — D > F > frozen — but the margin between D and F
is now measured as mostly unavailable rather than assumed to transfer.

---

## E85 — Inverting the difficulty model into a private risk ranking

E74 fitted the margin model forwards to get a private point estimate. This runs
it backwards: *given a configuration, where is it most likely to fail?* No score
follows. `scripts/76_private_risk_ranking.py`.

The headroom fraction, per entry, over the five public epitopes:

| epitope | tcrdist3 | headroom | frozen f | D f | F f |
|---|---:|---:|---:|---:|---:|
| GLCTLVAML | 0.690 | 0.310 | 0.166 | 0.216 | 0.178 |
| NLVPMVATV | 0.809 | 0.191 | 0.293 | **0.335** | 0.304 |
| TTDPSFLGRY | 0.826 | 0.174 | 0.081 | 0.134 | 0.135 |
| GILGFVFTL | 0.871 | 0.129 | **−0.013** | 0.118 | 0.084 |
| YLQPRTFLL | 0.928 | 0.072 | 0.075 | 0.097 | 0.093 |
| **mean** | | | 0.120 | **0.180** | 0.159 |

The frozen entry's fraction on GILGFVFTL is **negative** — it scores below
tcrdist3's published 0.871 there — which is worth recording as the one place the
frozen entry is beaten by the published baseline.

**Band risk and fraction risk are separate failure modes** and the second is
larger. Varying tcrdist3 across its published interval moves D's LTDEMIAQY
prediction over [0.721, 0.738]; varying the fraction across its five observed
values moves it over [0.702, 0.780], three times wider. LTDEMIAQY carries the
widest total interval of the three (0.079 against 0.044), as both its band
position and its headroom make it.

**Sampling risk**, Hanley–McNeil at the private cohort sizes:

| epitope | rows | SE @20% pos | SE @40% | SE @50% |
|---|---:|---:|---:|---:|
| CINGVCWTV | 243 | 0.0360 | 0.0270 | 0.0251 |
| LLWNGPMAV | 229 | 0.0372 | 0.0278 | 0.0259 |
| LTDEMIAQY | **130** | **0.0603** | 0.0465 | 0.0441 |

LTDEMIAQY is first on every axis: worst band, widest fraction interval, smallest
cohort.

### The private cohorts are the opposite of GLCTLVAML

The mechanism risk was expected to be GLC-like over-convergence, since E20 and
E72 established that failure mode and the family diagnostics are label-free and
therefore legitimate on the private rows. It is the other way round:

| epitope | rows | singleton % | largest family | close pairs % | nn3 IQR | split |
|---|---:|---:|---:|---:|---:|---|
| CINGVCWTV | 243 | **97.5** | **2** | 0.01 | 0.1677 | private |
| LLWNGPMAV | 229 | 95.6 | 4 | 0.03 | 0.2282 | private |
| LTDEMIAQY | 130 | 85.4 | 7 | 0.23 | **0.1239** | private |
| TTDPSFLGRY | 387 | 95.3 | 5 | 0.02 | 0.1360 | public |
| NLVPMVATV | 366 | 81.7 | 17 | 0.15 | 0.2253 | public |
| YLQPRTFLL | 419 | 81.6 | 12 | 0.11 | 0.3852 | public |
| GILGFVFTL | 670 | 52.4 | 71 | 0.70 | 0.3445 | public |
| GLCTLVAML | 192 | 44.8 | 32 | 1.27 | 0.3540 | public |

**The private split has almost no family structure at all.** CINGVCWTV's largest
family is two TCRs. There is no over-convergence risk there; the risk, if any, is
the opposite — nothing to be typical *of*. That TTDPSFLGRY sits in the same
regime (95.3% singleton, nn3 IQR 0.1360) and still scores 0.8493 is the
reassurance, and it also makes TTD the natural analogue for the private split
rather than GIL or GLC.

**None of the four diagnostics predicts per-epitope AUC** over five points
(|ρ| ≤ 0.400), so this table is descriptive. It says what the private cohorts
look like, not what they will score.

### The structure-matched fraction, which fails its own holdout

The table above suggests a sharper model than E74's: if the private cohorts
resemble TTDPSFLGRY, take the fraction from the structurally *nearest* public
epitope rather than averaging over all five, which mixes in two family-rich
cohorts the private split does not resemble. Matching on six standardised
structural features gives CIN → TTD and LLW → TTD (f = 0.134) and LTD → NLV
(f = 0.335), and raises D's private estimate from 0.809 to **0.820**.

**It is rejected.** Under E74's own leave-one-epitope-out test:

| margin model | held-out MAE |
|---|---:|
| uniform | 0.0301 |
| **headroom (E74)** | **0.0157** |
| structure-matched | 0.0261 |

The failure is visible: holding out TTDPSFLGRY matches it to NLVPMVATV and
predicts +0.0583 against an actual +0.0233, a 2.5× overshoot. Nearest-neighbour
over four candidates is simply unstable, and the appeal of the story is not
evidence. **E74's 0.809 stands** — this is the same shape as E30's withdrawn
projection and E49's third part, a mechanism that reads well and does not
survive its holdout.

What does survive from the matching exercise is a risk flag: LTDEMIAQY's nearest
public analogue sits at standardised distance **3.01**, against 1.68 for the
other two, in a space where the public epitopes are 1.91–5.82 apart. **LTDEMIAQY
has no good analogue in the public set on any axis** — band, cohort size or
cohort structure.

### Does the ordering survive an adversarial fraction?

| entry | mean f | best f | worst f | macro @mean | @best | @worst |
|---|---:|---:|---:|---:|---:|---:|
| frozen | 0.120 | 0.293 | −0.013 | 0.7948 | 0.8350 | 0.7636 |
| **D** | 0.180 | 0.335 | 0.097 | **0.8086** | 0.8448 | **0.7892** |
| F | 0.159 | 0.304 | 0.084 | 0.8037 | 0.8377 | 0.7862 |

D at its **worst** observed fraction (0.7892) still exceeds the best rival at its
**best** (F, 0.7862). The ordering therefore does not depend on which fraction is
the right one — though the margin is +0.0030, well inside every uncertainty on
this page, so this is robustness of the ordering and not of the gap.

---

## E86 — A similarity map of the eight cohorts: which public epitope is each private one like?

E85 needed one number from this — each private epitope's nearest public analogue
— and turned up something worth drawing in full. `scripts/77_epitope_similarity_map.py`,
figure at `reports/figures/epitope_similarity_map.png`.

"Alike" is not one thing, so three axes are computed, and the first result is
that **they are mutually orthogonal**:

| axis | axis | ρ over the 28 epitope pairs |
|---|---|---:|
| cohort structure | TCR repertoire | **−0.073** |
| cohort structure | peptide | +0.032 |
| TCR repertoire | peptide | +0.133 |

So there is no single answer to "which epitope resembles which", and in
particular **cohorts whose TCRs sit in similar regions of sequence space are not
organised alike**. That is a point in favour of the sampling reading of E55 and
E76: family structure is a property of how a cohort was collected, not of the
receptors in it.

### 1. Cohort structure — the axis that matters, and the only one that groups

Standardised over all eight (E85 used the public five, since it was predicting
public folds; the public–public distances shift by at most 0.45 and no ordering
changes). Distance to each public epitope, nearest first:

| private | 1st | 2nd | 3rd | 4th | 5th |
|---|---|---|---|---|---|
| CINGVCWTV | **TTDP 1.40** | NLVP 3.09 | YLQP 4.04 | GILG 5.82 | GLCT 5.92 |
| LLWNGPMAV | **TTDP 1.47** | NLVP 2.17 | YLQP 2.98 | GILG 4.98 | GLCT 5.07 |
| LTDEMIAQY | NLVP 2.56 | TTDP 2.62 | YLQP 3.79 | GLCT 4.54 | GILG 5.19 |

The dendrogram splits cleanly in two, and **all three private cohorts land in
the same branch as exactly one public epitope**:

```
  singleton-dominated   LTDEMIAQY · TTDPSFLGRY · CINGVCWTV · LLWNGPMAV
  family-rich           GLCTLVAML · GILGFVFTL · NLVPMVATV · YLQPRTFLL
```

CINGVCWTV and LLWNGPMAV are closer to TTDPSFLGRY (1.40, 1.47) than any two
public epitopes are to each other except NLVP–YLQP (1.80). **TTDPSFLGRY is the
private split's public twin**, and it scores 0.8493 under D.

The clustering is **not** recovering HLA. LTDEMIAQY and TTDPSFLGRY are the only
two HLA-A\*01:01 cohorts, and LTDEMIAQY's nearest neighbour is nonetheless the
A\*02:01 NLVPMVATV, at 2.56 against 2.62.

**LTDEMIAQY is the one with no twin.** Its nearest analogue sits at 2.56 in a
space where the public epitopes span 1.80–5.59, and E85 already put it first on
every other risk axis — worst published band, widest fraction interval, smallest
cohort at 130 rows. Every axis independently nominates the same epitope.

### 2. TCR repertoire — flat, which is itself the result

Symmetrised Chamfer distance in SCEPTR space (for every TCR in A, the cosine
distance to its nearest neighbour in B, averaged, and the same back). All 28
pairs fall in **0.485–0.600**, a span of 0.115 on an axis where neighbours
*inside* a cohort sit at 0.12–0.40. Between cohorts, every epitope is about
equally far from every other. The two mild exceptions are NLVP–GILG (0.485) and
CINGVCWTV–GLCTLVAML (0.513).

This is what makes the orthogonality above interpretable rather than puzzling:
there is barely any repertoire-level structure *to* correlate with.

### 3. Antigen — carries nothing, as expected

BLOSUM62 similarity normalised by self-score, range −0.42 to +0.27. E75 already
established that six of the eight are HLA-A\*02:01 nonamers, so this axis is
nearly flat by construction; it is drawn to show that it supplies no grouping
rather than because it was expected to.

### The inference the figure invites, and which does not hold

The bar panel labels each public epitope with tcrdist3's published AUROC, which
invites reading *closer to the private cohorts* as *easier*. It does not:

| public epitope | mean distance to the three private | tcrdist3 |
|---|---:|---:|
| TTDPSFLGRY | 1.83 | 0.826 |
| NLVPMVATV | 2.61 | 0.809 |
| YLQPRTFLL | 3.60 | 0.928 |
| GLCTLVAML | 5.18 | 0.690 |
| GILGFVFTL | 5.33 | 0.871 |

Spearman **+0.100**. The two most distant cohorts are the easiest and nearly the
hardest. Structural resemblance says what the private cohorts *look like*, not
what they will score — which is the same conclusion E85 reached when the
structure-matched margin model failed its leave-one-epitope-out test at MAE
0.0261 against 0.0157. **This map must not be used to re-weight anything.**

---

## E87 — Should the public epitopes be weighted by resemblance to the private ones?

E86 measured that CINGVCWTV and LLWNGPMAV are structural near-twins of
TTDPSFLGRY (1.40, 1.47) while GILGFVFTL and GLCTLVAML are the two public
cohorts least like the private split (mean distance 5.33, 5.18). The move that
follows is textbook covariate-shift correction: when choosing a configuration,
stop averaging the five public epitopes equally and upweight the ones the
private split actually resembles. `scripts/78_proximity_weighted_selection.py`.

E85 rejected a *different* use of the same distances — predicting the private
AUC **level** — and that says nothing about this. Predicting a level and
predicting *which configuration wins* are different questions, and only the
second is what selection needs.

The prior is against, and quantitatively so: E15 measured that cross-validated
selection over this same 323-configuration bank costs **−0.0048** against fixing
the configuration in advance, because with five folds the selection step has
worse signal-to-noise than the thing being selected. Proximity weighting reduces
the effective number of folds further, so it must beat equal weighting by enough
to pay for that as well.

### 1. Configuration transfer does not track structural distance

For each pair of epitopes, the Spearman correlation between their
323-vectors of configuration AUC — does the same configuration win on both?

| pair | structural distance | config transfer ρ |
|---|---:|---:|
| NLVP–YLQP | **1.80** | **+0.014** |
| NLVP–TTDP | 2.18 | +0.542 |
| GILG–YLQP | 2.92 | +0.652 |
| GILG–GLCT | 3.07 | +0.063 |
| GILG–NLVP | 3.10 | +0.223 |
| TTDP–YLQP | 3.44 | +0.191 |
| GLCT–YLQP | 3.90 | +0.484 |
| GLCT–NLVP | 3.96 | −0.068 |
| GILG–TTDP | **4.92** | **+0.520** |
| GLCT–TTDP | 5.59 | −0.247 |

Spearman **−0.345, exact p 0.331** over the ten pairs (enumerated over all
3,628,800 permutations, per E49). The sign is the one the hypothesis wants and
the magnitude is nothing. The counterexamples are what matter: the structurally
*closest* public pair, NLVP–YLQP at 1.80, has essentially **zero** configuration
transfer (+0.014), while GILG–TTDP at 4.92 — nearly the most distant pair —
transfers at +0.520.

**And it is not that the wrong distance was picked.** Every label-free
similarity available was tested against the same transfer correlations:

| label-free similarity | ρ with config transfer | exact p |
|---|---:|---:|
| repertoire Chamfer (E86) | −0.479 | 0.165 |
| \|log cohort size ratio\| | −0.406 | 0.248 |
| \|nn3 IQR difference\| | −0.394 | 0.261 |
| structural distance (E86) | −0.345 | 0.331 |
| \|singleton % difference\| | −0.139 | 0.706 |
| \|log largest-family ratio\| | +0.442 | 0.204 |
| HLA mismatch | +0.071 | 0.915 |

Nothing reaches significance across seven tests at n=10. Five of the seven carry
the hypothesised sign, which is weak directional evidence that the idea is not
absurd — but the strongest candidate is the repertoire Chamfer distance, which
E86 showed is the *flattest* axis of the three, so even that is not a mechanism.

**The obstacle is structural, not a matter of finding a better metric.** The
quantity that would predict configuration transfer is configuration transfer
itself, and computing it needs labels. The private split has none.

### 2. The payoff: weighting makes nested LOEO monotonically worse

Nested leave-one-epitope-out with the *inner* selection objective weighted by
`exp(-d/τ)` toward the held-out epitope — exactly the procedure that would be
run against the private split, simulated where the answer is known. `τ=∞` is
E15's equal-weight nested LOEO; `τ=0` selects on the single nearest epitope.

| τ | nested LOEO | vs equal weight | distinct configs chosen |
|---|---:|---:|---:|
| **∞ (equal)** | **0.8460** | — | 3 |
| 8 | 0.8463 | +0.0003 | 2 |
| 4 | 0.8444 | −0.0016 | 3 |
| 2 | 0.8444 | −0.0016 | 3 |
| 1 | 0.8412 | −0.0048 | 4 |
| 0.5 | 0.8417 | −0.0043 | 3 |
| 0.25 | 0.8386 | **−0.0074** | 4 |
| 0 (nearest only) | 0.8410 | −0.0050 | 3 |

Oracle, for scale: 0.8508.

The only non-negative point is τ=8, whose weights are 0.19–0.22 — equal
weighting with rounding error, and +0.0003 is noise. Everything with real
weighting loses, and **τ=1 costs −0.0048, which is E15's selection cost over
again**: weighting pays the selection penalty a second time, which is what
reducing the effective fold count buys. Per fold, nearest-neighbour-only loses
on four of the five.

### What this closes, and what it does not

**Closed:** using cohort resemblance to reweight the selection objective. It was
measured across eight temperatures and seven similarity metrics and is negative
at every setting that actually changes the weights.

**Not closed, and still the right use of E86:** cohort resemblance as a *risk
description*. LTDEMIAQY being first on every independent axis — published band,
headroom-fraction interval, cohort size, analogue distance — remains worth
knowing. The distances say what the private cohorts look like. They do not say
which configuration will win there, and now that has been tested rather than
assumed in either direction.
