# Classifying singletons only

Labelled singletons: 1281 across five epitopes. The splice reproduces the frozen score exactly when fed itself.

## Single features within labelled singletons

| feature                       |   mean_auc |   worst_epitope | neighbour_free   |
|:------------------------------|-----------:|----------------:|:-----------------|
| frozen                        |     0.7941 |          0.6549 | False            |
| nn3_dist                      |     0.7793 |          0.5994 | False            |
| nn1_dist                      |     0.7711 |          0.5953 | False            |
| density                       |     0.7702 |          0.6043 | False            |
| cons_ms_iptm_mean             |     0.7292 |          0.644  | True             |
| cons_ms_tcr_pep_iptm_mean     |     0.7276 |          0.6373 | True             |
| complex_ipde                  |     0.7142 |          0.6259 | True             |
| iptm_tcr_pep_max              |     0.7125 |          0.6176 | True             |
| iptm                          |     0.7084 |          0.6167 | True             |
| iptm_tcr_pep_mean             |     0.7076 |          0.6113 | True             |
| ptm                           |     0.7061 |          0.6167 | True             |
| iptm_a_pep                    |     0.7057 |          0.5962 | True             |
| iptm_b_pep                    |     0.7055 |          0.6244 | True             |
| iptm_a_mhc                    |     0.7054 |          0.5881 | True             |
| cons_ms_confidence_score_sd   |     0.7053 |          0.6078 | True             |
| iptm_b_mhc                    |     0.7051 |          0.6188 | True             |
| cons_ms_iptm_sd               |     0.7047 |          0.6088 | True             |
| cons_ms_confidence_score_mean |     0.7029 |          0.6302 | True             |
| plddt_peptide                 |     0.701  |          0.5997 | True             |
| confidence_score              |     0.6998 |          0.6345 | True             |

## Singleton-only classifiers, nested LOEO

| features            | model             | top_k   |   singleton_auc |   worst_epitope |   spliced_macro |
|:--------------------|:------------------|:--------|----------------:|----------------:|----------------:|
| all features        | logistic L2       | 5       |          0.7949 |          0.6396 |          0.851  |
| all features        | logistic L2 (C=1) | 5       |          0.7943 |          0.6431 |          0.8492 |
| all features        | logistic L2 (C=1) | 15      |          0.7921 |          0.6473 |          0.8474 |
| all features        | logistic L2       | 15      |          0.7888 |          0.6403 |          0.8455 |
| all features        | boosting depth 2  | 5       |          0.7863 |          0.6604 |          0.843  |
| all features        | logistic L2       | all     |          0.7927 |          0.693  |          0.8413 |
| all features        | boosting depth 2  | 15      |          0.782  |          0.6653 |          0.8405 |
| all features        | logistic L2 (C=1) | all     |          0.7864 |          0.6951 |          0.8379 |
| all features        | boosting depth 2  | all     |          0.7811 |          0.6888 |          0.8379 |
| neighbour-free only | logistic L2       | 15      |          0.7277 |          0.6348 |          0.7995 |
| neighbour-free only | logistic L2 (C=1) | 5       |          0.7246 |          0.6407 |          0.7992 |
| neighbour-free only | logistic L2       | 5       |          0.725  |          0.6387 |          0.799  |
| neighbour-free only | logistic L2 (C=1) | 15      |          0.7234 |          0.6215 |          0.7962 |
| neighbour-free only | logistic L2       | all     |          0.7189 |          0.5986 |          0.7912 |
| neighbour-free only | boosting depth 2  | 15      |          0.7034 |          0.6146 |          0.7869 |
| neighbour-free only | boosting depth 2  | 5       |          0.7021 |          0.6228 |          0.7859 |
| neighbour-free only | logistic L2 (C=1) | all     |          0.704  |          0.589  |          0.7819 |
| neighbour-free only | boosting depth 2  | all     |          0.6771 |          0.593  |          0.7708 |

- frozen entry within singletons **0.7941**
- frozen entry macro **0.8498**
- best spliced macro **0.8510** (+0.0012)
