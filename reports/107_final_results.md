# Final results — the private leaderboard, and what it did and did not resolve

Written 2026-11-19, after the private leaderboard was published. Everything
below was read from the Kaggle API, not from the web page.

## 1. Outcome

The competition closed 2026-11-02 15:00 UTC with 163 teams and 165 scored
entries. Our six entries:

| Kaggle team | public | private | method |
|---|---:|---:|---|
| HeidAI05 | 0.860 | **0.819** | D — four-term unsupervised rank average |
| Jielin2026 | — | 0.813 | not identified in this repository |
| saarufer | 0.882 | 0.811 | C1 — rank consensus of S4 and FHB30 |
| HeidAI | 0.871 | 0.810 | S4 — D plus nine refinements |
| Jielin Wei | 0.861 | 0.809 | S09 — frozen anchor plus SCEPTR reference support |
| HeidAI04 | 0.875 | 0.794 | FHB30 — bounded pairwise correction on a frozen anchor |

D, S4 and C1 are this repository's; S09 and FHB30 come from collaborators and
are described here only as far as needed to read the table.

**Our best private score is 0.819, from D.** The pre-registered estimate
recorded in D's submission description on 2026-10-26 was 0.798–0.802; the
realised value was 0.819, so the estimate was conservative by about 0.02.
C1's pre-registered estimate was ~0.79 against a realised 0.811, the same
direction and the same magnitude.

The organisers' tcrdist3 baseline scored 0.768 (rank 92). D beats it by
+0.051.

## 2. The public score did not predict the private score

| | public | private | public rank among our entries | private rank |
|---|---:|---:|---:|---:|
| C1 | 0.882 | 0.811 | 1 | 3 |
| FHB30 | 0.875 | 0.794 | 2 | 6 |
| S4 | 0.871 | 0.810 | 3 | 4 |
| S09 | 0.861 | 0.809 | 4 | 5 |
| D | 0.860 | **0.819** | 5 | **1** |

The ordering is close to reversed. The two most elaborate entries — S4's nine
refinements and C1's cross-method consensus — bought public score and lost
private score relative to the plainest member of the family.

The consensus result deserves a precise statement rather than a dramatic one.
On a bootstrap of the public rows shaped like the private split (three cohorts
of 243/229/130 rows), C1's expected margin over D was **+0.022 ± 0.006**. The
realised private margin was **−0.008**. The swing is about five standard
deviations of the paired difference, so *the relationship between the two
genuinely changed between the public and private peptides*. What the private
split alone cannot establish is the weaker claim that D is better than C1
there: that 0.008 gap is only 1.3 paired standard deviations.

This is the cleanest confirmation the project has that the standing rule —
**prefer the configuration with fewer fitted quantities when the evaluation
set is three unseen cohorts** — was the right one, and that the +0.0086
eighteen-fold gain recorded for C1 was in-sample with respect to its members'
own fitted weights, exactly as its submission description warned.

## 3. What the leaderboard can and cannot resolve

Measured by bootstrapping the public rows into a private-shaped split
(`scripts/129_leaderboard_resolution.py`):

- A single team's macro-AUC has a row-sampling standard deviation of
  **0.017**. That wobble is common to every team, since all are scored on the
  same 602 rows, so it does not by itself reorder the board.
- What reorders the board is the standard deviation of the *difference*
  between two teams, which depends on how correlated their predictions are:

  | two entries | ρ (within-peptide) | sd(difference) | gap needed for 95% |
  |---|---:|---:|---:|
  | near-identical (D vs F) | 0.998 | 0.002 | 0.004 |
  | same family (D vs C1) | 0.959 | 0.006 | 0.012 |
  | related (D vs density alone) | 0.929 | 0.009 | 0.018 |
  | unrelated teams | ~0 | 0.024 | **0.047** |

Applying the last row: the best score among teams that did not recover labels
is 0.846, ours is 0.819, a gap of 0.027 — **below** the 0.047 two unrelated
methods would need to be distinguished. Counting by team, after collapsing
multi-account clusters and removing the label-recovery entries, **35 of 107
teams sit within 0.047 of the top**, and 16 within one standard deviation. The
top third of the clean field is one large tie.

Our position in that field, taking the best of our six accounts: **19th of 107
teams**, tied with two others at 0.819; 25th if each multi-account group is
represented by its median rather than its best.

## 4. The top of the board is label recovery, not prediction

Seven entries scored 0.953–0.999, then the field drops to 0.846 and continues
smoothly down to 0.467. A score of 0.999 over 602 rows and three unseen
peptides is not attainable by prediction: the IMMREP25 unseen-peptide winner
scored 0.60, and this project's own 150-experiment search saturated near 0.86
on the public peptides.

**No leak exists inside the released files.** Checked on 2026-11-09:

- `tcr_id` decomposes as plate/row/well/index/peptide (2 × 192 wells, matching
  the well- and plate-specific primer tables of the source preprint). Scoring
  each field by its leave-one-peptide-out positive rate gives AUC 0.481–0.506
  — pure noise.
- All 17 metadata columns of `test.csv`, plus missingness indicators: best is
  `TRAV_IMGT` at 0.559.
- Cross-peptide clonotype overlap: `clonotype_aa_screened` shares **zero**
  keys between peptides; CDR3β links only 3 private rows to public rows.
- The released TCRvdb file is 2,034 rows × 24 columns, exactly the public
  half, with a single `pos`/`neg`/`ambiguous` label triple and no cross-pool
  reactivity profile.

**The private labels are derivable from public data.** NCBI BioProject
[PRJNA1464244](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1464244),
submitted by "MATCHMAKERS Grand Challenge" and registered **2026-05-10** — five
months before the competition opened — holds 51 SRA samples of full-length
unique TCR molecule ONT sequencing from the pooled VDJdb-10 screens. The
library names include:

```
CINGVCWTV_CD69_pos_rep_1 / _2 / _3
LLWNGPMAV_CD69_pos_rep_1 / _2 / _3
LTDEMIAQY_CD69_pos_rep_1 / _2 / _3
TTDPSFLGRY_CD69_pos_rep_1..3 · GILGFVFTL_CD69_pos_rep_1..3
NLVPMVATV_CD69_pos / _neg · GLCTLVAML_CD69_pos / _neg / _pMHC_multimer
YLQPRTFLL_CD69_pos / _pMHC_multimer
VDJdb-10_Jurkat_baseline_for_LTD · _for_LTD_and_TTD · _for_GIL · _for_GIL_and_NLV
VDJdb-10_plasmid_1 / _2
```

All three private peptides are present, with matched baselines. The recovery
path is complete without any private information: `test.csv` supplies
`full_seq_alpha_aa` and `full_seq_beta_aa` to map reads against, the released
TCRvdb labels for the five public peptides calibrate the enrichment threshold,
and the same pipeline then applies to CIN, LLW and LTD. A score of 0.999 is
consistent with a join missing one or two ambiguous rows.

**This corrects E2.** Its follow-up concluded "no deposit containing them was
found… the private leaderboard is sound". That was wrong: only the TCRvdb
GitHub repository and its changelog were checked, not SRA. The deposit had
been public since May. We have not attempted any reconstruction and attribute
nothing to any team; the point is that the benchmark's private split was not
in fact hidden.

## 5. Method summary

Two lines, meeting in C1.

### D — unsupervised cohort density (public 0.860, private 0.819)

Premise: receptors genuinely reactive to a peptide converge within their
cohort, whereas entries that enter a database through technical artefacts are
effectively random draws from the repertoire. Each row is therefore scored by
how typical its receptor is of the cohort annotated to the same peptide.

Four components, each converted to a within-peptide percentile rank and
averaged with **equal weights and no fitted parameters**:

| component | definition |
|---|---|
| `density3` | mean of within-peptide ranks of exponential kernel densities in default-SCEPTR (bandwidth 0.1), small-SCEPTR (0.1) and paired α+β tcrdist3 (20) distance spaces; self-distances excluded |
| `pooled_support` | negative unweighted mean of the eight smallest default-SCEPTR cosine distances in a same-peptide reference pool: the competition cohort plus paired human VDJdb records (release 2026-06-03, exact peptide, deduplicated by paired CDR3, own row excluded) |
| `full_sequence_support` | negative unweighted mean of the eight smallest combined α+β full-V-domain alignment distances within the same-peptide cohort |
| `structure_support` | negative minimum TCR–peptide predicted aligned error from Boltz-2 complexes (TCR V domains, peptide, reduced class-I α1/α2 platform, empty MSAs) — model confidence, not measured affinity |

No binding labels, no peptide supervision, and no VDJdb-trained binding
predictor: the last was ruled out in E121 as circular for a task auditing
VDJdb (the strongest public candidate overlapped the test set 100% and carried
872 inverted labels).

A single SCEPTR density alone scores 0.8396 of D's 0.8601 — **94.3% of the
above-chance signal**. The ladder is 0.8396 → +structure 0.8498 → +V-domain
0.8557 → +pool 0.8601. Three of the four terms are distance measures and they
are the three strongest; structure is the weakest at 0.7565 standalone.

The premise fails where within-cohort convergence is absent. Density alone
reaches AUC 0.927 on YLQPRTFLL and 0.677 on GLCTLVAML, where the validated and
non-validated distributions sit on top of each other. Nothing tested across
150 experiments recovered that cohort.

### S4 — D plus nine refinements (public 0.871, private 0.810)

D is re-ranked, then two structural corrections are applied — a
peptide-adaptive kernel bandwidth, and a 60-value contact fingerprint (2
chains × 10 peptide positions × 3 descriptors, L2-normalised, scored as the
negative mean cosine distance to eight same-peptide neighbours) — followed by
seven sequential rank blends at 5–15% each: interface polar fraction,
interface charge product, non-CDR3 interface polarity, paired CDR3-length
frequency, interface packing, TwoNN local intrinsic dimensionality, and
mutual-nearest-neighbour support. Prescribed rank fusion with fixed export
weights, not a trained classifier; the residual effective weight on the D
baseline is 0.5018.

### C1 — rank consensus (public 0.882, private 0.811)

Equal-weight average of the within-peptide percentile ranks of S4 and a
collaborator's FHB30, which fits a bounded supervised correction
`s = a + 0.05·tanh(xᵀw)` on a frozen sequence-similarity anchor. No fitted
parameters in the consensus; the weights were fixed at one half a priori. The
rationale was variance reduction rather than added information: at a
within-peptide rank correlation of 0.861, an equal-weight average retains
(1+ρ)/2 = 0.931 of a single member's variance. Section 2 records what happened
to that argument on the private split.

## 6. What this changes for the next benchmark of this kind

1. **A public leaderboard whose labels are published measures a join.** 19
   teams reached ≥0.984 here. Standing on it carried no information, and the
   project dropped it from consideration after E2.
2. **A private split is only hidden if its upstream data is.** Check the
   sequence archives, not just the curated release.
3. **Report the spread.** With three cohorts and 602 rows, differences below
   0.047 between unrelated methods are not resolvable. A ranked table without
   that caveat invites over-reading.
4. **Prefer fewer fitted quantities.** Within one family, evaluated on three
   unseen cohorts, every elaboration we added helped the public score and hurt
   the private one.
