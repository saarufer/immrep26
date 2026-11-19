# Which source studies are reliable, and why

Compiled 2026-11-16. Sources: `reports/56_study_false_positive_rate.md` (E90,
`scripts/81_study_false_positive_rate.py`), NCBI BioProject PRJNA1464244 and
PRJNA1068303, Europe PMC / PubMed records for each cited study.

The competition label is the outcome of **re-testing** a published annotation
under uniform conditions, so `pos = 0` is, to a first approximation, a false
positive in the source literature. That makes the dataset a direct audit of the
papers that fed VDJdb.

---

## 0. What the Schumacher group actually validated

**Scope: exactly the eight competition epitopes, and nothing else.**

NCBI BioProject **PRJNA1464244** (registered 2026-05-10, "MATCHMAKERS Grand
Challenge", 51 SRA samples, ONT full-length unique-TCR-molecule sequencing)
contains CD69-sorted reactivity screens for:

| screened | assays deposited |
|---|---|
| GILGFVFTL, NLVPMVATV, TTDPSFLGRY, CINGVCWTV, LLWNGPMAV, LTDEMIAQY | CD69-positive, 3 replicates each |
| GLCTLVAML, YLQPRTFLL | CD69-positive and CD69-negative, plus pMHC-multimer-positive sorts |
| — | matched Jurkat baselines (`..._baseline_for_GIL`, `_for_GIL_and_NLV`, `_for_LTD`, `_for_LTD_and_TTD`, `_for_pMHC_multimer_sort`) and `VDJdb-10_plasmid_1/2` library controls |

The VDJdb-10 library was assembled from the ten most densely sampled VDJdb
epitopes. **The two remaining members, SPRWYFYYL and ATDALMTGF, were not
screened.** So there is no functional revalidation of any VDJdb target outside
the competition's eight.

A separate, earlier NKI deposit, **PRJNA1068303** (2024-01-24, "Functional
discovery of tumor-reactive T cell receptors by massively parallel library
synthesis and screening: NKIRTIL063 neoantigen screen"), covers tumour
neoantigens. That is TCR *discovery*, not revalidation of published annotations,
and shares no targets with this benchmark.

**What the validation axis actually is.** T-RAP synthesises each receptor at
uniform expression in a Jurkat reporter line; the readout is CD69 upregulation
on peptide challenge. So "validated" means *this αβ pair signals against this
pMHC in a standardised reporter*. It is not primary-T-cell function, not
affinity, and not a measurement in the donor's own genetic background.
Approximately half of tested VDJdb entries did not validate; the deposit
description gives 56% validating.

---

## 1. The measurement

23 (study, epitope) cells with n ≥ 10, covering **15 studies and 1,559 TCRs**.
1,615 of the 1,782 scorable rows carry exactly one VDJdb reference and form the
primary analysis; the 162 multi-reference rows are analysed separately because
replication is itself a candidate quality signal.

The dominant confound is the epitope — base validation rates run 0.262
(NLVPMVATV) to 0.632 (GLCTLVAML), and groups specialise. Every cell is therefore
compared against its epitope's base rate computed **with that study removed**.

Dispersion statistic **383.1**; within-epitope label-permutation null mean 57.3,
95th percentile 75.5, maximum 94.7 over 2,000 draws. **p < 0.0005.** On a single
epitope, study validation rates span **0.084 to 0.923**.

---

## 2. Reliable and unreliable, named

`excess` = validation rate minus the leave-this-study-out epitope base rate.

### Unreliable

| study | year, venue | epitope | n | validated | expected | excess | z |
|---|---|---|---:|---:|---:|---:|---:|
| Chour et al., `doi:10.1101/2020.05.04.20085779` | 2020 medRxiv, never published | NLVPMVATV | 143 | **0.084** | 0.424 | −0.340 | **−8.23** |
| Sant et al., PMID:29997621 | 2018 Front Immunol | GILGFVFTL | 81 | 0.358 | 0.626 | −0.268 | −4.99 |
| Rowntree et al., PMID:35750048 | 2022 Immunity | TTDPSFLGRY | 129 | 0.380 | 0.568 | −0.188 | −4.32 |
| Nguyen et al., `doi:10.1016/j.xcrm.2023.101017` | 2023 Cell Rep Med | YLQPRTFLL | 39 | 0.154 | 0.435 | −0.281 | −3.54 |
| Glanville et al., PMID:28636589 | 2017 Nature (GLIPH) | NLVPMVATV | 36 | 0.028 | 0.294 | −0.267 | −3.51 |
| Dash et al., PMID:28636592 | 2017 Nature (TCRdist) | GILGFVFTL | 137 | 0.482 | 0.620 | −0.138 | −3.34 |
| Glanville et al., PMID:28636589 | 2017 Nature | GILGFVFTL | 75 | 0.427 | 0.608 | −0.181 | −3.22 |
| Minervina et al., PMID:35383307 | 2022 Nat Immunol | YLQPRTFLL | 171 | 0.339 | 0.461 | −0.122 | −3.19 |
| Rowntree et al., PMID:35750048 | 2022 Immunity | YLQPRTFLL | 47 | 0.234 | 0.430 | −0.196 | −2.71 |
| Trautmann et al., PMID:16237109 | 2005 J Immunol | NLVPMVATV | 14 | **0.000** | 0.275 | −0.275 | −2.31 |

### Reliable

| study | year, venue | epitope | n | validated | expected | excess | z |
|---|---|---|---:|---:|---:|---:|---:|
| Francis et al., PMID:34793243 | 2022 Sci Immunol | YLQPRTFLL | 64 | 0.734 | 0.338 | +0.397 | **+6.71** |
| Francis et al., PMID:34793243 | 2022 Sci Immunol | GILGFVFTL | 75 | 0.853 | 0.522 | +0.332 | +5.75 |
| Minervina et al., PMID:35383307 | 2022 Nat Immunol | TTDPSFLGRY | 213 | 0.568 | 0.380 | +0.188 | +5.66 |
| Wang et al., PMID:29483513 | 2018 Nat Commun | GILGFVFTL | 52 | **0.923** | 0.532 | +0.391 | +5.65 |
| Koo et al., PMID:36711524 | 2023 bioRxiv (nanovial) | NLVPMVATV | 16 | 0.750 | 0.235 | +0.515 | +4.86 |
| Dash et al., PMID:28636592 | 2017 Nature | NLVPMVATV | 53 | 0.472 | 0.218 | +0.254 | +4.48 |
| `github.com` (unattributed deposit) | — | NLVPMVATV | 15 | 0.533 | 0.248 | +0.285 | +2.56 |
| Nguyen et al., PMID:33951417 | 2021 Immunity | YLQPRTFLL | 46 | 0.543 | 0.386 | +0.157 | +2.19 |
| Francis et al., PMID:34793243 | 2022 Sci Immunol | GLCTLVAML | 51 | 0.706 | 0.594 | +0.112 | +1.63 |
| Grant et al., PMID:27645996 | 2016 JBC | GILGFVFTL | 10 | 0.800 | 0.572 | +0.228 | +1.45 |

### The single most important caveat

**This is not a ranking of laboratories.** Three groups appear on both sides:

| study | epitope | excess | epitope | excess |
|---|---|---:|---|---:|
| Dash et al. 2017 | NLVPMVATV | **+0.254** | GILGFVFTL | **−0.138** |
| Minervina et al. 2022 | TTDPSFLGRY | **+0.188** | YLQPRTFLL | **−0.122** |
| Grant et al. 2016 | GILGFVFTL | +0.228 | GLCTLVAML | −0.066 |

The unit that validates well or badly is the **(study, epitope) assay instance**,
not the group. A laboratory's annotations for the pMHC it had good reagents and a
strong population for are reliable; the same laboratory's annotations for a
secondary epitope in the same experiment often are not.

---

## 3. Why — the reasons that survive checking

**1. Independent replication is the strongest single predictor.**

| epitope | 1 reference | rate | >1 reference | rate |
|---|---:|---:|---:|---:|
| GILGFVFTL | 445 | 0.578 | 129 | **0.961** |
| GLCTLVAML | 152 | 0.632 | 13 | **1.000** |
| NLVPMVATV | 301 | 0.262 | 14 | **0.714** |
| YLQPRTFLL | 375 | 0.405 | 6 | **1.000** |

Epitope-stratified Mantel–Haenszel **OR 17.1**. It is not a citation count
standing in for convergence: within-epitope ρ with top-8 density is +0.200,
against a typicality gate of 0.6378. The mechanism is plain — two groups with
different donors, reagents and sorting pipelines rarely produce the same wrong
αβ-to-pMHC assignment.

**2. Paired-chain evidence quality: amplicon-seq beats Sanger.**

| level | n | validated | expected | excess | p | distinct studies | largest share |
|---|---:|---:|---:|---:|---:|---:|---:|
| `amplicon-seq` | 380 | 0.613 | 0.444 | **+0.169** | <0.0005 | 9 | 51.8% |
| `sanger` | 701 | 0.462 | 0.552 | **−0.090** | <0.0005 | 12 | 36.2% |

These are the only method levels spread across enough independent groups to be
read as a method effect rather than one study wearing a method label. The
mechanism is the αβ pairing: Sanger-era work sequences a limited number of sorted
clones and infers the pairing through a cloning step, whereas amplicon-based
single-cell work observes both chains in the same droplet with UMI support and
depth. A wrong pairing produces a receptor that never had the claimed
specificity, which is exactly what a revalidation screen detects.

**3. Studies whose own evidence was functional validate best.** Koo et al. 2023
screened by **nanovial-based affinity *and function***, and its NLVPMVATV
annotations validate at 0.750 against an expected 0.235 (+0.515) — the largest
positive excess in the table. Annotations whose original evidence was multimer
binding alone are being re-tested on a different axis (signalling); annotations
already filtered on function are not.

**4. Secondary epitopes in multiplexed panels are the worst single source.** The
largest and worst contributor, Chour et al. 2020, is a **SARS-CoV-2** study that
nonetheless deposited **143 CMV NLVPMVATV** annotations, of which **8.4%**
survived. The pattern is a multiplexed multimer panel in which the paper's focus
channel receives the scrutiny and the off-focus channels do not. This is also the
most likely explanation for the within-study sign flips in §2.

**5. Rare or impaired populations.** Nguyen et al. 2023 sampled vaccine responses
in B-cell-depleted haematological-malignancy patients: 0.154 against 0.435. When
the true antigen-specific population is small, the ratio of sort noise to real
events rises, and the surviving annotations are correspondingly weaker.

**6. Peer review is not nothing, but it is not the variable either.** The worst
contributor is an unreviewed preprint still unpublished six years later; but
*Nature* (Glanville, Dash), *Immunity* (Rowntree) and *Nat Immunol* (Minervina)
all appear on the negative side, and a *bioRxiv* preprint (Koo) is the single
best contributor. Venue does not order reliability.

---

## 4. Why — the reasons that look real and are not

**VDJdb's own confidence score carries no usable ordering.** Score 3, its
highest, validates at 0.147 against an expected 0.402 — but **91.7% of score-3
rows come from one study**, and scores 0, 1, 2 run +0.106, +0.023, +0.287 with no
monotone trend. The score describes documentation completeness, not whether the
annotation is correct.

**Several alarming "method effects" are one study wearing a method label:**

| level | excess | distinct studies | largest study's share |
|---|---:|---:|---:|
| `rna-seq` | −0.304 | 2 | **96.6%** |
| `tissue = lung` | −0.371 | 1 | **100%** |
| `vdjdb.score = 3` | −0.255 | 5 | **91.7%** |
| `cell.subset` missing | −0.175 | 3 | 94.3% |
| `dextramer-sort` | −0.059 | 3 | 72.6% |

**Era alone buys nothing.** Spearman between PMID-as-recency and the excess z is
**+0.204, permutation p = 0.392** over the 20 cells carrying a PMID. Median PMIDs
do differ by platform (sanger ≈ 29.5M ≈ 2018, amplicon-seq ≈ 34.8M ≈ 2021), so
platform and era cannot be fully separated at 15 studies — but recency on its own
does not predict validation.

---

## 5. Limits

- Everything above is measured on the **five public epitopes only**. The private
  epitopes have a different study composition: CINGVCWTV matches only 79 of its
  243 rows in VDJdb at all, and cross-study replication is **0/243 (CIN),
  1/229 (LLW), 0/130 (LTD)** — and 0/387 for TTDPSFLGRY on the public side.
- **None of this transfers to a per-receptor feature.** The strongest signal,
  replication, is structurally absent from the private split, because each
  private epitope was contributed by only two or three papers. The second,
  sequencing platform, is a property of 15 studies rather than of a TCR. The
  competition entry was not changed on the basis of any of it.
- Attribution is restricted to the 1,615 single-reference rows; the 23 cells are
  those with n ≥ 10.
- `pos = 0` equals "did not signal in the T-RAP Jurkat reporter". A receptor that
  is genuinely specific but depends on a co-receptor, a CD8 contribution or a
  donor-specific context would be scored as a false positive here.
