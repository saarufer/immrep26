# H5: antigen-conditioned alpha-beta pairing compatibility

## Gate A — cross-chain dependence, real against shuffled

| epitope    |   n |   pmi_real |   pmi_shuffled |   pmi_stratified |   cont_real |   cont_shuffled |   pmi_excess |   pmi_excess_vs_strata |
|:-----------|----:|-----------:|---------------:|-----------------:|------------:|----------------:|-------------:|-----------------------:|
| CINGVCWTV  | 243 |     3.4099 |         2.8078 |           3.4114 |     -2.4221 |         -2.4173 |       0.6021 |                -0.0015 |
| GILGFVFTL  | 670 |     2.1271 |         1.0128 |           2.0968 |     -2.1256 |         -2.159  |       1.1143 |                 0.0303 |
| GLCTLVAML  | 192 |     5.076  |         2.4584 |           5.0169 |     -2.4182 |         -2.4689 |       2.6176 |                 0.0591 |
| LLWNGPMAV  | 229 |     2.9114 |         2.3965 |           2.875  |     -2.3259 |         -2.3304 |       0.5148 |                 0.0364 |
| LTDEMIAQY  | 130 |     5.3307 |         4.557  |           5.3307 |     -2.4375 |         -2.4692 |       0.7737 |                 0      |
| NLVPMVATV  | 366 |     3.7827 |         2.5985 |           3.7762 |     -2.4276 |         -2.4361 |       1.1842 |                 0.0065 |
| TTDPSFLGRY | 387 |     2.3946 |         1.9708 |           2.3772 |     -2.3481 |         -2.3509 |       0.4238 |                 0.0173 |
| YLQPRTFLL  | 419 |     2.6445 |         1.774  |           2.6218 |     -2.3086 |         -2.3432 |       0.8705 |                 0.0227 |

## Gate B — separation of functional labels

| score   |   mean_auc |   GILGFVFTL |   GLCTLVAML |   NLVPMVATV |   TTDPSFLGRY |   YLQPRTFLL |
|:--------|-----------:|------------:|------------:|------------:|-------------:|------------:|
| pmi     |     0.44   |      0.3582 |      0.3494 |      0.648  |       0.4656 |      0.3789 |
| cont    |     0.6336 |      0.7845 |      0.5846 |      0.5605 |       0.6371 |      0.6013 |

## Gate C — viability and nested gain

| score   | viable   | private_levels                |   nested_gain |
|:--------|:---------|:------------------------------|--------------:|
| pmi     | True     | CINGV:230 LLWNG:189 LTDEM:113 |       -0.0013 |
| cont    | True     | CINGV:243 LLWNG:229 LTDEM:130 |       -0.0007 |

Best nested gain -0.0007 against +0.02: does not clear the bar; the entry stands.
