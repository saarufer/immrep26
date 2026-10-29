# The admission standard: what must be answered before a candidate is built

> 我希望下一次提出的不是"再试一个模型"，而是能够先回答：它究竟凭什么掌握了
> D 尚未掌握、同时又可以在 Private 上获得的信息？

Every gate this project runs is **post hoc**: the typicality gate, the private
variance gate, residual-on-D, the permutation null, the 16 pseudo-epitope
folds. All of them require the candidate to exist first. Between E96 and E122
that cost eleven experiments, of which the honest yield was one feature.

This document is the **prior** gate. Three questions, answered in writing with
numbers, before any code is written.

---

## Q1 — Name the channel. Not the functional.

> Which physical or informational channel does this candidate read that is not
> already read by D **or by any feature already blended into it**?

A channel is a source of measurement. A functional is a formula applied to a
channel. **D's ceiling has never been set by the choice of formula**, and the
project has four independent proofs of this:

- **E98** — three candidate mechanisms were orthogonal (ρ +0.006 / −0.016 /
  +0.155) and summed at **98% of additive**, reaching +0.0009. The ceiling is
  not overlap, so it cannot be relieved by de-overlapping.
- **E103** — residual-on-D's floor is **0.5098**: that is what a term *already
  inside D* scores after residualisation. Every reformulation of the same
  channel lands near it.
- **E106** — the permutation null's **max is +0.0203**. Any fitted gain under
  that is reachable by shuffled labels.
- **E118** — 16 pseudo-epitope folds turned the fitted +0.0103 into
  **−0.0017, p 0.683**. Formula search transfers at zero.

So "a different functional of the same distance matrix" is **not an answer**.
E99 (hubness), E100 (TwoNN), E101 (cross-epitope LID) each proposed one, each
was defensible, each returned ≤ +0.0019 and none survived E118.

**Evidence required:** name the channel, and state the *predicted* residual-on-D
as a number before running. The floor is 0.5098 and the best any candidate has
ever reached is 0.5606. A proposal that cannot predict where it will land in
that interval does not have a mechanism — it has a hope. Run it anyway if it is
cheap, but do not call the result a discovery.

**Validity condition (E112):** the residual is only readable when
`|ρ(residual, D)| ≈ 0`. A feature with few distinct values ties under
percentile rank and the residual comes back monotone in D — this produced a
spurious 0.6626 and was caught only because a 0.6626 residual and a −0.0005
blend cannot both be true. All gated features have since been re-audited;
all leak ≤ 0.041.

---

## Q2 — Show the variance, on the three private cohorts, before fitting.

> Is the channel present **and variable** on CINGVCWTV (243), LLWNGPMAV (229)
> and LTDEMIAQY (130)?

Availability is not the test; *within-cohort variance* is. The metric is mean
per-epitope AUC, so a channel that is constant inside an epitope contributes
exactly zero to that epitope regardless of how well it separates globally.

- **E91** — platform looked discriminative until the strata were checked.
- **E122** — the Atlas raised LTDEMIAQY's pool from 5 to 133 rows and left
  CINGVCWTV at 0; macro gain **+0.0000**, and the one analogue that could be
  checked (YLQPRTFLL, 54→121) scored **−0.0008**.
- **Pool rows below ~60 are dead.** CINGVCWTV has 11, LTDEMIAQY has 5.
- **The private cohorts are singletons.** There are no families there, so
  anything that needs to fit per-epitope has one cohort to fit on, and
  TTDPSFLGRY — not GILGFVFTL or GLCTLVAML — is the public analogue.

**Evidence required:** the per-cohort standard deviation, computed on the
private rows, printed before any blend weight is chosen. `polar_frac` passes
this (sd 0.010–0.014 in all three). E113's V-gene prior failed it: 54% pair
coverage.

---

## Q3 — Prove independence from the process that generated the labels.

> The task is *which VDJdb annotations are wrong*. Is the candidate downstream
> of the annotation process it is being asked to audit?

This is the sharpest lesson in the project and it was nearly missed.

**E121** — Hi-TpH's Level-III training set contains **100% of the 2,636
competition receptors**, and **872 of the 1,782 scorable rows are IMMREP26
negatives it was explicitly trained to call binders.** Its predictions are not
independent evidence about those rows; they are a readout of the labelling it
is supposed to be checking, with the sign inverted on half the set.

This generalises past one model. **Any binding predictor trained on VDJdb —
UpTCR, TCRBinder, and most of the published field — is disqualified for this
task by construction**, no matter how well it validates elsewhere. The
disqualification is not about quality. It is about circularity.

**Evidence required:** the candidate's training corpus or derivation path,
checked against the 2,636 rows by sequence. Boltz-2 passes (trained on PDB
coordinates, not binding assays). E113 failed: it aggregated the competition's
own labels by V gene, then asked them to predict those labels.

---

## The fourth failure mode: relieving a constraint that was never binding

A candidate can pass all three questions and still hold nothing. **E116 is the
worked example, and it is worth more than the eleven failures combined because
it failed while being right.**

The proposal: average `polar_frac` over all five Boltz models instead of one.
Q1 — arguably fine, it is the same channel sampled better. Q2 — fine, present
everywhere. Q3 — fine, no labels involved. And measurement theory made a
falsifiable quantitative prediction: Spearman–Brown, from a single-structure
reliability of 0.812.

| m | reliability | S–B predicts | residual on D | blend gain |
|---:|---:|---:|---:|---:|
| 1 | 0.812 | 0.812 | 0.5586 | +0.0016 |
| 3 | 0.945 | 0.928 | **0.5606** | +0.0036 |
| 6 | 0.968 | 0.963 | 0.5581 | +0.0047 |

**The prediction was correct — within 0.031 at every m — and the candidate
still held nothing.** Reliability rose 19%. The residual moved **−0.0005**,
crossing **3.6%** of the floor-to-ceiling span, with ρ(reliability, residual)
= **−0.700**. The blend gain tripled, which looks like success until it is
divided by E118's measured null sd for this estimator: **+0.0047 = 1.04 σ.**

Noise was never the binding constraint. `polar_frac` carries ≈0.56 of
epitope-independent AUC that D lacks, and **all of it was readable from a
single structure.** The gap to E115's 0.5650 ceiling is not sampling error and
cannot be bought with compute.

**Evidence required, added as a standing question:** which constraint does this
relieve — information, or noise? And is that constraint the one that binds?

---

## Retrospective validation

Applied to every candidate since E96, scoring only on what was knowable
*before* each was built:

| candidate | Q1 channel | Q2 private | Q3 independent | would be | actually returned |
|---|---|---|---|---|---|
| E96 lenfreq | **fail** — length is already in SCEPTR | pass | pass | **killed** | +0.0005 |
| E97 per-epitope bandwidth | **fail** — reparameterisation | partial | pass | **killed** | +0.0007, 1/3 transfers |
| E99 hubness | **fail** — same matrix | pass | pass | **killed** | −0.0007 |
| E100 TwoNN LID | **fail** — same matrix | pass | pass | **killed** | +0.0019 |
| E101 cross-epitope LID | **fail** — same matrix | pass | pass | **killed** | +0.0007, unstable |
| E113 V-gene prior | pass | **fail** — 54% coverage | **fail** — label source | **killed** | below chance every fold |
| E121 Hi-TpH | pass | pass | **fail** — 100% overlap | **killed before download** | disqualified |
| E103 germline docking | pass — coordinates | pass | pass | run | ≤0, reducible to m50 |
| E104 `polar_frac` | pass — coordinates + elements | pass — sd 0.010–0.014 | pass | run | **+0.0016, the one survivor** |
| E112 pMHC side | pass — the other interface | pass | pass | run | nothing cleared |
| E122 Motif Atlas | pass — new receptors | partial — LTD only | pass | run | +0.0000 |

**Seven of eleven killed before any code. All seven returned ≤ +0.0019, and
E118 put the whole stack at −0.0017.** The four it would have admitted include
every experiment that was worth running, and `polar_frac` — the only feature in
the project carrying information not reducible to something D already had.

The gate is not conservative. It is *correct* on the sample it was built from,
which is weak evidence, so it is stated here with its derivation so it can be
argued with rather than obeyed.

---

## Current channel inventory

What D reads: receptor sequence (SCEPTR), reference-pool composition
(m30pool, m50), structural contact fingerprint (struct03).

| channel | status |
|---|---|
| Atomic coordinates beyond the contact fingerprint | **mined.** E103/E104/E112; `polar_frac` is the survivor and E116 shows it is fully extracted |
| Folding-model confidence (ipTM, pLDDT, PAE — 35 columns, all 2,636 rows) | **inside the typicality well.** 0.7387 on task A, **0.5362** on the real task; `iptm_mhc_pep` — an interface with no TCR — scores 0.647; as a residual **−0.0193** |
| Provenance / metadata / platform | **dead on private** |
| External receptor repertoires | **measured empty** on CINGVCWTV (E122); needs ~60 pool rows |
| VDJdb-supervised binding models | **Q3-disqualified by construction** (E121) |
| Epitope-conditioned anything | **dead by construction** — private cohorts are singletons |
| Competition labels aggregated by any receptor property | **Q3-disqualified** (E113) |

No channel currently in hand answers all three questions. That is the output of
the gate today, and it is a statement about the inventory, not about the time
remaining — the deadline is 2026-11-02 15:00 UTC. The productive use of what is
left is to **look for a channel**, not another functional on the channels above.
