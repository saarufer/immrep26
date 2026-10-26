# Why a TCR is a singleton, and whether the radius is principled

## Part 1a - predicting singleton status from the TCR itself

| feature               | cause                     |   mean_auc |    min |    max |
|:----------------------|:--------------------------|-----------:|-------:|-------:|
| cdr3_len_total        | chain-feature atypicality |     0.6347 | 0.4339 | 0.7959 |
| ms_iptm_sd            | structurally odd          |     0.3698 | 0.2046 | 0.7518 |
| typicality            | chain-feature atypicality |     0.3641 | 0.2452 | 0.4748 |
| complex_plddt         | structurally odd          |     0.3635 | 0.2813 | 0.5077 |
| consist_tcr_rmsd_mean | structurally odd          |     0.3567 | 0.1934 | 0.5976 |
| log10_pgen_paired     | hard to generate          |     0.3351 | 0.192  | 0.4515 |
| pae_tcr_pep_min       | structurally odd          |     0.3185 | 0.2306 | 0.4062 |
| ms_iptm_mean          | structurally odd          |     0.3054 | 0.1987 | 0.5239 |
| iptm                  | structurally odd          |     0.2986 | 0.2152 | 0.3728 |

All nine together, in-sample within epitope: **0.850** (an upper bound).

## Part 1b - the sampling explanation

| epitope_aa   |   cohort_n |   singleton_rate |
|:-------------|-----------:|-----------------:|
| CINGVCWTV    |        243 |           0.9753 |
| GILGFVFTL    |        670 |           0.5239 |
| GLCTLVAML    |        192 |           0.4479 |
| LLWNGPMAV    |        229 |           0.9563 |
| LTDEMIAQY    |        130 |           0.8538 |
| NLVPMVATV    |        366 |           0.8169 |
| TTDPSFLGRY   |        387 |           0.9535 |
| YLQPRTFLL    |        419 |           0.8162 |

Singleton rate against cohort size: rho **-0.309** (n=8).

## Part 2a - distance distribution and the background

| epitope_aa   |   nn3_q10 |   nn3_median |   nn3_q90 |   bg_nn3_median |   enrich_median |
|:-------------|----------:|-------------:|----------:|----------------:|----------------:|
| CINGVCWTV    |    0.31   |       0.4471 |    0.6102 |          0.5148 |          1.1622 |
| GILGFVFTL    |    0.0392 |       0.1899 |    0.5614 |          0.4792 |          2.4893 |
| GLCTLVAML    |    0.0484 |       0.2135 |    0.5675 |          0.4691 |          2.0204 |
| LLWNGPMAV    |    0.272  |       0.4142 |    0.6445 |          0.5015 |          1.21   |
| LTDEMIAQY    |    0.236  |       0.5946 |    0.6775 |          0.5125 |          0.8773 |
| NLVPMVATV    |    0.124  |       0.5351 |    0.6258 |          0.5027 |          0.975  |
| TTDPSFLGRY   |    0.2697 |       0.5324 |    0.6215 |          0.5106 |          0.9663 |
| YLQPRTFLL    |    0.1112 |       0.452  |    0.62   |          0.5123 |          1.1266 |

| epitope    |   half_background_radius |   frac_within_012 |   bg_frac_within_012 |   enrichment_at_012 |
|:-----------|-------------------------:|------------------:|---------------------:|--------------------:|
| CINGVCWTV  |                     0.21 |            0      |               0.0041 |                 0   |
| GILGFVFTL  |                     0.47 |            0.3925 |               0.003  |               131.5 |
| GLCTLVAML  |                     0.45 |            0.3333 |               0      |               inf   |
| LLWNGPMAV  |                     0.49 |            0.0131 |               0      |               inf   |
| LTDEMIAQY  |                     0.45 |            0.0462 |               0      |               inf   |
| NLVPMVATV  |                     0.45 |            0.0984 |               0.0082 |                12   |
| TTDPSFLGRY |                     0.46 |            0.0233 |               0      |               inf   |
| YLQPRTFLL  |                     0.49 |            0.1217 |               0      |               inf   |

## Part 2b - validation rate by nn3 decile

|   decile |   n |   nn3_mean |   rate_pooled |   rate_epitope_mean |   epitopes_declining |
|---------:|----:|-----------:|--------------:|--------------------:|---------------------:|
|        0 | 181 |     0.0778 |        0.9227 |              0.9056 |                    0 |
|        1 | 178 |     0.1489 |        0.8371 |              0.8208 |                    3 |
|        2 | 177 |     0.2196 |        0.8079 |              0.7874 |                    3 |
|        3 | 178 |     0.2823 |        0.7247 |              0.701  |                    2 |
|        4 | 178 |     0.3505 |        0.5337 |              0.5174 |                    4 |
|        5 | 177 |     0.4105 |        0.4689 |              0.447  |                    5 |
|        6 | 177 |     0.4736 |        0.2429 |              0.2307 |                    5 |
|        7 | 178 |     0.5238 |        0.3146 |              0.3173 |                    1 |
|        8 | 177 |     0.576  |        0.1695 |              0.1974 |                    4 |
|        9 | 181 |     0.6336 |        0.0829 |              0.1167 |                    5 |

## Part 2b2 - label-free enrichment against the published bands

| label_free_index   |   rho_vs_published_auroc |      p |   n |
|:-------------------|-------------------------:|-------:|----:|
| enrich_median      |                   0.3086 | 0.4571 |   8 |
| frac_within_012    |                   0.2315 | 0.5813 |   8 |
| nn3_median         |                  -0.3086 | 0.4571 |   8 |

## Part 2c - binary indicator against the continuous variable

|   radius |   binary_auc |   clustered_frac |
|---------:|-------------:|-----------------:|
|     0.02 |       0.5223 |           0.0303 |
|     0.04 |       0.561  |           0.0952 |
|     0.06 |       0.5911 |           0.1445 |
|     0.08 |       0.6205 |           0.1768 |
|     0.1  |       0.6523 |           0.2136 |
|     0.12 |       0.6725 |           0.236  |
|     0.14 |       0.6865 |           0.2561 |
|     0.16 |       0.6918 |           0.2769 |
|     0.18 |       0.7051 |           0.2929 |
|     0.2  |       0.7127 |           0.3137 |
|     0.22 |       0.7167 |           0.3335 |
|     0.24 |       0.7317 |           0.3528 |
|     0.26 |       0.7408 |           0.3729 |
|     0.28 |       0.7467 |           0.3964 |
|     0.3  |       0.751  |           0.4127 |
|     0.32 |       0.746  |           0.4336 |
|     0.34 |       0.749  |           0.4579 |
|     0.36 |       0.7544 |           0.4784 |
|     0.38 |       0.7464 |           0.4954 |
|     0.4  |       0.7376 |           0.5296 |

- best binary: radius 0.36, AUC **0.7544**
- continuous nn3_dist: **0.8334**
- difference **+0.0790**
