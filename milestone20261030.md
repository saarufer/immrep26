# IMMREP26 — Milestone, 2026-10-30

**Competition:** Kaggle *IMMREP26: Validation of public TCR-pMHC data*
**Deadline:** 2026-11-02 15:00 UTC
**Task:** given a VDJdb TCR–pMHC annotation, predict whether it survives
experimental revalidation in TCRvdb
**Metric:** mean per-epitope ROC AUC, truncated by Kaggle to 3 decimals
(granularity **0.0005**)

**Data:** 2,636 rows / 8 epitopes. 2,034 public (1,782 scorable, 51.1%
positive), **602 private across 3 epitopes** — CINGVCWTV 243, LLWNGPMAV 229,
LTDEMIAQY 130. The public rows carry visible labels, so the public leaderboard
is in-sample for anything fitted on them; only the 602 private rows decide the
ranking. Our harness reproduces the leaderboard exactly (E47, reconfirmed by
E133: harness +0.0116 against the leaderboard's +0.011).

**Status: 135 logged experiments (E1–E145). The entry is D, unchanged since
2026-10-27.**

---

## Part 1 — The submissions

### Ranked by public macro, which is the only measured number

| # | submission | public macro | sha256 (16) | principle |
|---:|---|---:|---|---|
| 1 | `S4_D_stack9` | **0.8716** | `d6b588c4c69f57d7` | D restructured (density3 bandwidth 0→0.5, structure slot `pae`→0.7·pae+0.3·fingerprint) then **7 components mixed sequentially**: polar_frac 0.15, charge_product 0.15, pf_germline 0.10, lenfreq 0.05, packing 0.05, TwoNN 0.05, mutual-kNN 0.10. **Only ~50% of its rank is D**, and 0.20 of the 0.65 added weight is `charge_product` (residual 0.5542) and `packing` (0.5494), both of which **failed** the residual-on-D gate. |
| 2 | `S2_D_stack5` | 0.8675 | `04afccd0b0584fe5` | D restructured + 3 components: lenfreq 0.05, TwoNN 0.05, polar_frac 0.15. The only added term that cleared the residual gate is `polar_frac` (0.5586 against E83's 0.5570 benchmark). |
| 3 | `S5_D_ttd_tuned` | 0.8664 | `be42ea85dc0adb1f` | **Same three components as S2**, weights fitted on **TTDPSFLGRY alone** — the public cohort nearest the private ones by pool size (25 rows). Measured transfer of one-epitope fitting: **+0.0007 ± 0.0088** over 20 (fit, holdout) pairs. |
| 4 | `S1_D_plus_polar` | 0.8643 | `4f088170ec34874d` | D + `polar_frac` at w=0.15. **One component, the least selected** — the fraction of N/O atoms among TCR atoms within 9 Å of pMHC, read off Boltz-2 structures. E109 localised it to the **germline loops facing the MHC**, not CDR3 facing the peptide; E111 showed it is a V-gene property but not the one D already carries; E115 replicated it on an independent Boltz run (ρ +0.856, reliability 0.812). |
| — | **`D` = `d_density3_m30pool_m50_struct`** | **0.8600** (7th) | `64b3dd4076a51fbd` | **The entry.** `(density3 + m30pool + m50 + structure) / 4`, uniform weights. `density3` = mean rank of kernel densities in three spaces (SCEPTR default σ=0.1, SCEPTR small σ=0.1, tcrdist3 σ=20); `m30pool` = top-30 neighbour distance against the epitope's reference pool; `m50` = whole-V-domain top-k alignment similarity; `structure` = Boltz-2 `pae_tcr_pep_min`. |

Per-epitope, public:

| submission | GILGFVFTL | GLCTLVAML | NLVPMVATV | TTDPSFLGRY | YLQPRTFLL |
|---|---:|---:|---:|---:|---:|
| S4 | 0.8906 | **0.7929** | 0.8824 | **0.8695** | 0.9229 |
| S2 | 0.8864 | 0.7721 | 0.8883 | 0.8588 | 0.9318 |
| S5 | 0.8858 | 0.7679 | 0.8866 | 0.8595 | 0.9322 |
| S1 | 0.8850 | 0.7706 | 0.8814 | 0.8526 | 0.9320 |
| **D** | 0.8862 | 0.7568 | 0.8729 | 0.8493 | **0.9350** |

Remaining submissions: `S3_F_stack5` 0.8638, `S7_D_atlas_pool` 0.8601,
`f_..._nopool` 0.8557, `density4_struct03` 0.8525, `density3_struct03` 0.8498
(the earlier frozen entry), `density_sceptr_default_sigma0.1` 0.8396.

### Ranked by expected private score, which is what matters

**The private level is not identifiable.** The class-composition projection
behind the often-quoted 0.799 fails its own leave-one-epitope-out test — MAE
0.0833, worst error 0.2129 (GLCTLVAML predicted 0.954 against an actual 0.741)
— and its sensitivity surface spans **0.829–0.857**. Only *differences* are
estimable, and the level is shared by the whole family.

| submission | public vs D | 5-fold OOF | **18-fold OOF** | loses on |
|---|---:|---:|---:|---:|
| S4 | +0.0116 | +0.0116 | **−0.0017** | 1/5 |
| S2 | +0.0074 | +0.0048 | **−0.0045** | 1/5 |
| S5 | +0.0063 | **+0.0007** | — | 2/5 |
| S1 | +0.0043 | +0.0016 | **−0.0048** | 2/5 |
| **D** | — | — | — | **0/5** |

**The ranking flips with the estimator.** Under five epitope folds every fitted
candidate is positive; under eighteen pseudo-epitope folds — E118's better
estimator, null sd 0.0045 against 0.0059 — **every one goes negative**, and S4,
the most heavily fitted, comes out least negative. All values sit inside ±1
null sd, so **the top seven are a statistical tie and the ordering among them
is not determined by the data.**

The bottom four differ from D by **architecture** rather than fitted weights,
so their −0.0044 to −0.0205 carry across unshrunk and those rankings are solid.

### Why D is the entry

1. **18-fold out-of-fold**: D is the reference and every fitted candidate is
   below it (−0.0017 to −0.0080).
2. **Only submission that loses on none of the five public epitopes** (0/5;
   S4 loses 1, the rest lose 2).
3. **Only audited artefact**: E123 ran 14 checks, reproduced D from source to
   **1.11e-16**, and confirmed all four components carry variance in all three
   private cohorts.
4. **S4's public +0.0116 is in-sample for its own weight selection.** E133
   tested the one argument for S4 — ρ(D's per-epitope AUC, S4's gain) = −1.000,
   suggesting S4 helps exactly where D is weak, which would favour it since two
   private cohorts have pools of 11 and 5. Sixty stacks of **seven pure-noise
   components**, greedily fitted the same way, reproduce ρ = −1.000 in 5% of
   trials with a mean of −0.638: **the pattern is mechanical.** S4's gain is
   nevertheless real (3.5× the maximum noise reached) — real, but with its
   per-epitope allocation fitted, across-epitope sd 0.0179 = 1.5× its own mean.

### Portfolio, if several entries are available

| portfolio | E[max] vs D alone |
|---|---:|
| D alone | +0.00000 |
| **D + S4** | **+0.00054** |
| **D + S4 + S2** | **+0.00063** |
| D + S4 + S2 + S5 | +0.00069 |

**Everything is below Kaggle's 0.0005 truncation step.** Recommended: **D, S4,
S2** and stop; the fourth slot is worth +0.00003.

**S2 and S5 must not be split across two entrants** — private-row ρ **0.996**
and **0.0%** of private rows move more than 10 percentile points. That is one
submission entered twice.

All three of D, S4, S2 pass the full integrity audit (E141): exact columns,
2,636 rows, no duplicate ids, id set identical to `test.csv`, numeric with no
NaN or inf, no epitope constant.

**Action required: D must be selected manually.** It is 7th of 11 by public
score, so Kaggle's default — best public — will not pick it at any allowed
count. Entering D versus not is E134's unresolvable **[−0.0041, +0.0273]**,
which is **50× the entire portfolio effect**.

---

## Part 2 — What the 135 experiments established

### The governing fact

> **Real orthogonal information at AUC 0.54–0.56 exists and does not improve
> the ranking, because D is at 0.86 while every such channel is at 0.52–0.76
> about the same rows.**

Established from four independent directions:

| channel | residual on D | nested |
|---|---:|---:|
| `polar_frac` — structure (E104) | 0.5586 | +0.0016 |
| Minervina exo-TCR margin — wet-lab negatives (E126) | **0.5609** | −0.0003 |
| `pdockq1` — docking quality (E145) | 0.5406 | +0.0001 |
| excess convergence — composition null (E135) | 0.5070 | — |

The floor is **0.5098** (a term already inside D) and the record is 0.5609.

### Why the search is bounded, quantitatively

| route | bound |
|---|---|
| **add information** (E96–E140, E145) | clearing a measured +0.005 needs a **true** effect of +0.0086; clearing +0.02 needs **+0.0309**; the largest channel ever found is **+0.0016** (E131) |
| **remove estimation error — bandwidth** (E142) | response surface is flat: spread 0.0063 across the entire grid, D's fixed point near its top and beating fold-based reselection. +0.0001 ± 0.0013 |
| **remove estimation error — aggregation** (E143) | EV-neutral; break-even at an unobservable 0.036 flip rate against 0/20 observed. All values ≤ 0.0005 |
| **remove estimation error — distillation** (E144) | D is **not smooth** in sequence space (student ρ ceiling 0.71), so there is nothing to smooth |
| **portfolio / entry selection** (E141) | best 3-entry portfolio +0.00063, below the rounding step |

### Measurement instruments built, and what they caught

- **residual-on-D** (E103): residualise a candidate on D within epitope and
  rescore. Floor 0.5098.
- **the E112 tie check**: the residual is only readable when
  |ρ(residual, D)| < 0.1. Small-integer or heavily-missing features tie under
  percentile rank and the residual comes back monotone in D. **This caught five
  spurious "records"**: E112's 0.6626, E139's 0.6715/0.5900/0.5691/0.5591,
  E145's 0.6345/0.5822/0.5364, and two `is_zero` controls at 0.8422/0.6781.
  Without it, at least three log entries would have claimed a new record.
  For tied features the correct instrument is **stratification within
  D-quintiles**, not residual regression.
- **permutation null** (E106): shuffle labels within epitope and re-run the
  *entire* procedure including sign selection. Honest null sd 0.0076, max
  +0.0203. Nesting the sign step moved p from 0.082 to **0.200**.
- **18 pseudo-epitope folds** (E118): split each cohort by source study. Null
  sd 0.0045 against the 5-fold 0.0059. Turned the stack's fitted +0.0103 into
  **−0.0017, p 0.683**.
- **harness power** (E131): MDE +0.0031 for nested > 0; the binding constraint
  is the **between-epitope sd of D's own AUC, 0.0587**, which does not shrink
  with more rows.
- **the four admission questions** (E131/E132, `reports/87_candidate_admission.md`):
  name the channel not the functional, with a *pre-registered* residual;
  show per-cohort sd on the 602 private rows before fitting; prove independence
  from the labelling process; say which constraint is relieved. Applied
  retrospectively it kills 7 of 11 candidates before any code — all 7 returned
  ≤ +0.0019 — and admits all 4 that were worth running.

### The science, independent of the leaderboard

**Both competition tasks are antigen-naive.** An OLGA-trained SCEPTR separates
the cohorts at 2.67× against the real model's 2.60× (E89); CDR3 lengths alone
recover 47% (E95); E36 folded 300 TCRs against anchor-matched decoy peptides
and every target-minus-decoy margin collapsed to chance.

**In both tasks the density's power is low-order positional composition.**
E127: destroying residue co-occurrence while holding V/J, length and
per-position composition fixed makes epitope *attribution* **better** —
A−C = **−0.0175 ± 0.0038**, negative in 5/5 seeds. E135: on *revalidation*, the
positional-composition null recovers **95.6%** of the real density's
above-chance signal, and the excess sits **below the floor** (0.5070 against
0.5098). A collaborator's two-layer hypothesis — low-order composition decides
which cohort you resemble, high-order convergence decides whether you are
reliable evidence — is therefore **confirmed in layer 1 and refuted in layer 2**.

**`polar_frac` is a structurally-read, epitope-independent V-gene property**
on the MHC-facing germline loops (E104/E109/E111), replicating on an
independent Boltz run (E115: ρ +0.856, reliability 0.812, disattenuated ceiling
0.5650). E116 then showed that ceiling is **not** a noise ceiling: averaging 6
structures raised measured reliability 0.812→0.968, tracking Spearman–Brown to
within 0.031, while the residual moved **−0.0005**. All of its information is
readable from **one** structure.

**Any VDJdb-trained binding predictor is disqualified by construction** for a
task that audits VDJdb (E121). Hi-TPH's Level-III training set contains **100%**
of the 2,636 competition receptors, and **872 of 1,782** scorable rows are
IMMREP26 negatives it was trained to call binders. This generalises past one
model and was reconfirmed by the 2026-10-30 literature sweep, where Hi-TPH is
the headline new resource.

**Source-study false-positive rates run 0.084–0.923** (E90), with the worst
case explained mechanistically; replication across ≥2 independent studies is
OR 17 for validation but **absent on the private rows**.

**αβ pairing authenticity is not the problem** (E139/E140). `test.csv`'s unused
`clonotype_aa_screened` column records what TCRvdb actually screened, and it
equals the stored pair in **2,635/2,636** α and **2,636/2,636** β. Cell-level
dual-α is real at **9.05%** of cells in the 10x application note, and shows no
association with revalidation failure (+0.0143, 95% CI [−0.079, +0.098],
Fisher p = 1.000). The structural obstruction: the rows whose pairing can be
audited validate at **0.957** against 0.511 overall — the error class can only
be detected where the errors are not.

**Reusable negatives for the field:** OTS cannot supply a healthy naive human
repertoire (all 384 `includes_naive` units are COVID patients; the healthy,
untreated, no-cancer subset has zero naive-labelled units). VDJdb 2026's own
clustering requires ≥5 members per cluster and therefore covers **0 of 243**
CINGVCWTV rows.

### Dead ends, by category

- **Provenance / metadata**: dead on the private split (E3, E4, E90–E94, E92).
- **Different functionals of the same distance matrix**: hubness (E99), TwoNN
  (E100), cross-epitope LID (E101), per-epitope bandwidth (E97), lenfreq (E96)
  — all ≤ +0.0019, and E118 put the whole stack at −0.0017.
- **Pairwise / reweighting formulations**: all three closed; 60% of the loss
  sits in relative-less pairs (E22, E23, E50).
- **Cohort-similarity reweighting**: measured **negative** (E87) — it pays
  E15's −0.0048 selection cost twice.
- **Pool expansion**: the pool term is dead below ~60 same-epitope rows;
  CINGVCWTV has 11 and LTDEMIAQY has 5. The TCR Motif Atlas added 128 rows to
  LTDEMIAQY and scored **+0.0000** macro, with the one checkable analogue
  (YLQPRTFLL 54→121) at −0.0008 (E122).
- **Epitope-conditioned anything**: dead by construction — the private cohorts
  are singletons with no families, and TTDPSFLGRY is the public analogue, not
  GILGFVFTL or GLCTLVAML.
- **Label-representation change** (GuanRank-style): blocked upstream of its
  data — **D is not a trained model**, so a richer label representation has no
  objective to act on.
- **Generation-to-mature-repertoire selection** (E138): `Q_paired` residual
  0.5126 against a floor of 0.5098, pre-registration hit, and it **loses to
  the Pgen control it was meant to improve on** (`log10_pgen_beta` 0.5792
  standalone against `Q_beta`'s 0.5474).
- **Cross-reactivity breadth** (E128): direction favours promiscuity
  consistently across three panels and two species, all nine estimates
  positive, sign test p ≈ 0.125, **nothing resolved** — at n=16–20 the smallest
  detectable |ρ| is 0.651/0.591.

---

## Part 3 — Known risks on the private side

1. **The public mean is unrepresentative in the one dimension known to matter.**
   It is dominated by cohorts with 3,005 and 9,881 pool rows, while two of the
   three private cohorts have **11 and 5**. Splitting the one accepted channel
   by pool richness: `polar_frac` earns **+0.0070 on pool-rich** against
   **+0.0002 on pool-poor** — its gain comes from the cohorts least like the
   private ones. Not acted on: cohort-similarity reweighting is negative (E87)
   and pool-share attribution is a standing barred inference (E84).
2. **`structure` is the component most worth watching.** If one of D's four
   terms is inverted on an unseen cohort, the median-aggregation payout is
   +0.0063 for `structure` against +0.0007–0.0020 for the others, and
   `structure` has the lowest standalone AUC of the four (0.6294 on
   TTDPSFLGRY, 0.7270 across folds). First place to look if D underperforms.
3. **D carries 22 private ties** (CINGVCWTV 8, LLWNGPMAV 10, LTDEMIAQY 4) where
   S4 and S2 carry zero. Quantified: tied pairs are 0.0750% of comparable
   pairs, the ceiling for breaking every one correctly is +0.000375, and
   `polar_frac` resolves 11.7% of the pairs D cannot order — a realistic
   **+0.000044**, which is **11× smaller than one rounding step**. The only
   clean non-fitted argument available, and still too small to act on.
4. **"We measured no effect" is not "there is no effect."** Candidates whose
   true effect sat near +0.0005–0.0007 faced a procedure that reports a
   positive only **47%** of the time (E131). The decisions stay correct *under
   the bar* — a true +0.0007 is 29× below what +0.02 demands — but the two
   claims were conflated in earlier write-ups and are not any more.

---

## Part 4 — Artefacts

| | |
|---|---|
| master log | `reports/00_experiment_log.md` (135 experiments, ~10k lines) |
| admission standard | `reports/87_candidate_admission.md` |
| per-experiment reports | `reports/01`–`reports/101_*.md` |
| scripts | `scripts/01`–`scripts/123_*.py` |
| submissions | `submissions/*.csv`, sha256 recorded above |
| git | `saarufer/immrep26` (private), branch `master` |

**Not redistributable:** TCRvdb (`data/external/TCRvdb`, `data/external/tcrvdb`)
— "Redistribution or public sharing of the dataset, in whole or in part, is not
permitted". Gitignored and never committed. Cite doi:10.1016/j.immuni.2026.09.007.
