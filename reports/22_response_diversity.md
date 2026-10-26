# Does typicality work only where the response is motif-dominated?

## Part 1 - label-free diversity, rarefied to the smallest cohort

| epitope    |   n | split   |   largest_family_share |   family_simpson |   trbv_simpson |   len_simpson |   published_tcrdist_mid |
|:-----------|----:|:--------|-----------------------:|-----------------:|---------------:|--------------:|------------------------:|
| CINGVCWTV  | 243 | private |                 0.0124 |           0.9922 |         0.9441 |        0.8405 |                   0.815 |
| GILGFVFTL  | 670 | public  |                 0.1064 |           0.9707 |         0.5699 |        0.5136 |                   0.9   |
| GLCTLVAML  | 192 | public  |                 0.1654 |           0.9537 |         0.8945 |        0.7437 |                   0.67  |
| LLWNGPMAV  | 229 | private |                 0.0185 |           0.992  |         0.946  |        0.7489 |                   0.815 |
| LTDEMIAQY  | 130 | private |                 0.0538 |           0.9891 |         0.9464 |        0.8232 |                   0.67  |
| NLVPMVATV  | 366 | public  |                 0.0464 |           0.9896 |         0.9558 |        0.8413 |                   0.815 |
| TTDPSFLGRY | 387 | public  |                 0.0177 |           0.992  |         0.8841 |        0.8031 |                   0.815 |
| YLQPRTFLL  | 419 | public  |                 0.0364 |           0.9899 |         0.9238 |        0.783  |                   0.9   |

## Part 2 - between-epitope correlation (n=5, cannot be significant)

| index                              |   rho_vs_typicality |   p_exact |   p_familywise |   n |   rho_vs_tcrdist_published |
|:-----------------------------------|--------------------:|----------:|---------------:|----:|---------------------------:|
| largest_family_share               |                 0.1 |    0.95   |         1      |   5 |                   nan      |
| family_simpson                     |                -0.1 |    0.95   |         1      |   5 |                   nan      |
| trbv_simpson                       |                -1   |    0.0167 |         0.065  |   5 |                   nan      |
| len_simpson                        |                -0.7 |    0.2333 |         0.6545 |   5 |                   nan      |
| largest_family_share (GIL dropped) |                -0.4 |    0.75   |       nan      |   4 |                   nan      |
| family_simpson (GIL dropped)       |                 0.4 |    0.75   |       nan      |   4 |                   nan      |
| trbv_simpson (GIL dropped)         |                -1   |    0.0833 |       nan      |   4 |                   nan      |
| len_simpson (GIL dropped)          |                -0.4 |    0.75   |       nan      |   4 |                   nan      |
| largest_family_share               |               nan   |  nan      |       nan      |   8 |                    -0.2315 |
| trbv_simpson                       |               nan   |  nan      |       nan      |   8 |                    -0.3858 |

## Part 3 - within-epitope conditional AUC

| epitope    |   auc_none |   pairs_none |   auc_trbv |   pairs_trbv |   auc_family |   pairs_family |   auc_trbv_x_len |   pairs_trbv_x_len |   pair_retention_trbv |   auc_in_largest_family |   auc_is_dominant_trbv |    drop |   drop_lo |   drop_hi |
|:-----------|-----------:|-------------:|-----------:|-------------:|-------------:|---------------:|-----------------:|-------------------:|----------------------:|------------------------:|-----------------------:|--------:|----------:|----------:|
| GILGFVFTL  |     0.779  |        73914 |     0.664  |        19778 |       0.7743 |            257 |           0.7004 |              12493 |                0.2676 |                  0.5501 |                 0.7196 |  0.115  |    0.0556 |    0.1826 |
| GLCTLVAML  |     0.6412 |         6160 |     0.6763 |          519 |       0.5441 |             68 |           0.5163 |                246 |                0.0843 |                  0.6093 |                 0.5604 | -0.0351 |   -0.1609 |    0.1568 |
| NLVPMVATV  |     0.5178 |        20566 |     0.6198 |          643 |       0.6667 |              9 |           0.7137 |                124 |                0.0313 |                  0.5934 |                 0.601  | -0.1019 |   -0.1999 |    0.0452 |
| TTDPSFLGRY |     0.6674 |        29240 |     0.6836 |         2642 |     nan      |              0 |           0.7176 |                779 |                0.0904 |                  0.5118 |                 0.6392 | -0.0162 |   -0.0833 |    0.0592 |
| YLQPRTFLL  |     0.5836 |        35392 |     0.5853 |         2240 |       0.25   |              8 |           0.565  |                869 |                0.0633 |                  0.538  |                 0.5829 | -0.0017 |   -0.0703 |    0.0723 |

## Part 3b - the control: TRBV-conditioning drop, by score

| epitope    |   density alone |   submitted blend |   typicality |
|:-----------|----------------:|------------------:|-------------:|
| GILGFVFTL  |          0.1191 |            0.1409 |       0.115  |
| GLCTLVAML  |         -0.1406 |           -0.1295 |      -0.0351 |
| NLVPMVATV  |         -0.0242 |           -0.0169 |      -0.1019 |
| TTDPSFLGRY |         -0.0067 |            0.0241 |      -0.0162 |
| YLQPRTFLL  |          0.0298 |            0.0386 |      -0.0017 |

Specific to typicality on GILGFVFTL: **False**.

## Part 4 - oracle epitope-conditioned weight

- frozen entry **0.8498**
- oracle per-epitope weight **0.8523** (+0.0025)
- oracle single global weight 0.8503 (w=0.10)
- conditioning buys an oracle +0.0020
- private variance gate True (CINGV:243 LLWNG:229 LTDEM:130)
- oracle below the +0.02 bar, so no label-free estimate was fitted
