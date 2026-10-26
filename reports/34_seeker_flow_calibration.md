# F1 - SEEKER screen-to-flow calibration

## The conditional sampling

| tcr   |   candidates |   flow_measured |   measured_pct |   median_r1_r0_measured |   median_r1_r0_all | has_high_copy   |
|:------|-------------:|----------------:|---------------:|------------------------:|-------------------:|:----------------|
| TCR59 |        19322 |             192 |         0.0099 |                 125.807 |            14.8425 | True            |
| TCR63 |        16369 |             192 |         0.0117 |                 170.226 |            14.418  | False           |

Flow was run on top-ranked candidates, so every number below describes ranked candidates, not the library.

## Screening statistic against flow activation

| tcr   | predictor        |   n |   spearman_raw |   cv_isotonic |   cv_linear_log |
|:------|:-----------------|----:|---------------:|--------------:|----------------:|
| TCR59 | r1_r0_ratio      | 192 |        -0.4505 |        0.4483 |          0.448  |
| TCR59 | r2_r1_ratio      | 192 |         0.6865 |        0.6348 |          0.6824 |
| TCR59 | count_pct_round1 | 192 |        -0.4804 |        0.4616 |          0.4812 |
| TCR59 | count_pct_round2 | 192 |         0.3781 |        0.3451 |          0.3718 |
| TCR59 | count_round0     | 192 |         0.1243 |        0.0686 |          0.1058 |
| TCR63 | r1_r0_ratio      | 192 |         0.1427 |        0.066  |         -0.1866 |
| TCR63 | r2_r1_ratio      | 192 |         0.324  |        0.4271 |          0.3199 |
| TCR63 | count_pct_round1 | 192 |        -0.0455 |       -0.0993 |         -0.0237 |
| TCR63 | count_pct_round2 | 192 |         0.5372 |        0.4877 |          0.5205 |
| TCR63 | count_round0     | 192 |        -0.2202 |        0.1593 |          0.1816 |

## Single-copy against high-copy (TCR59)

| condition   |   n |   median_activation |   spearman_with_r1_r0 |
|:------------|----:|--------------------:|----------------------:|
| single      | 192 |               22.05 |               -0.4505 |
| high        | 192 |               17.35 |               -0.4041 |

Rank agreement between conditions rho **+0.923**; **0 of 192** peptides are inactive at single copy (<5%) and active at high copy (>25%).
