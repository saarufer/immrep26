# Pairwise AUC error decomposition

Mean per-epitope AUC is the probability a validated TCR outranks a
non-validating one from the same cohort, so 1 - AUC is a set of inverted
pairs. Grouping them by cluster state gives a mutually exclusive, complete
partition whose losses sum to 1 - AUC (asserted in the script).

## Share of pairs

| epitope    |   pos_clu_neg_sing |   diff_cluster |   same_cluster |   both_singleton |   pos_sing_neg_clu |
|:-----------|-------------------:|---------------:|---------------:|-----------------:|-------------------:|
| GILGFVFTL  |               50.2 |            8.5 |            0.3 |             34.8 |                6.1 |
| GLCTLVAML  |               43.1 |           17.7 |            1.1 |             26.6 |               11.6 |
| NLVPMVATV  |               54.1 |            3   |            0   |             40.6 |                2.3 |
| TTDPSFLGRY |                5.9 |            0   |            0   |             94.1 |                0   |
| YLQPRTFLL  |               39.5 |            1.6 |            0   |             56.5 |                2.4 |
| MEAN       |               38.6 |            6.2 |            0.3 |             50.5 |                4.5 |

## AUC loss contributed

| epitope    |   pos_clu_neg_sing |   diff_cluster |   same_cluster |   both_singleton |   pos_sing_neg_clu |   1-AUC |
|:-----------|-------------------:|---------------:|---------------:|-----------------:|-------------------:|--------:|
| GILGFVFTL  |             0.0169 |         0.0193 |         0.0013 |           0.0602 |             0.033  |  0.1307 |
| GLCTLVAML  |             0.0214 |         0.0446 |         0.0024 |           0.0951 |             0.095  |  0.2586 |
| NLVPMVATV  |             0.0062 |         0.0079 |         0      |           0.1022 |             0.0188 |  0.1351 |
| TTDPSFLGRY |             0.0001 |         0      |         0      |           0.1597 |             0      |  0.1598 |
| YLQPRTFLL  |             0.0034 |         0.003  |         0.0001 |           0.048  |             0.012  |  0.0666 |
| MEAN       |             0.0096 |         0.015  |         0.0008 |           0.0931 |             0.0318 |  0.1502 |

## Error rate within class

| epitope    |   pos_clu_neg_sing |   diff_cluster |   same_cluster |   both_singleton |   pos_sing_neg_clu |
|:-----------|-------------------:|---------------:|---------------:|-----------------:|-------------------:|
| GILGFVFTL  |              0.034 |          0.227 |          0.405 |            0.173 |              0.54  |
| GLCTLVAML  |              0.05  |          0.253 |          0.221 |            0.358 |              0.819 |
| NLVPMVATV  |              0.011 |          0.263 |          0.111 |            0.252 |              0.825 |
| TTDPSFLGRY |              0.002 |        nan     |        nan     |            0.17  |            nan     |
| YLQPRTFLL  |              0.009 |          0.182 |          0.5   |            0.085 |              0.509 |
| MEAN       |              0.021 |          0.231 |          0.309 |            0.207 |              0.673 |
