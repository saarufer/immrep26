# Cohort thinning: does density survive a sparser repertoire?

Public cohorts are randomly thinned to simulate the private split's lower
candidate coverage, with real labels retained for scoring. `thin_reference`
drops TCRs from the reference set only and still scores every labelled row,
isolating coverage from the noise of scoring fewer rows.

## AUC vs retained fraction (thin_reference, density only)

| epitope    |   0.12 |    0.2 |    0.3 |   0.45 |    0.6 |    0.8 |    1.0 |
|:-----------|-------:|-------:|-------:|-------:|-------:|-------:|-------:|
| GILGFVFTL  | 0.8648 | 0.8679 | 0.8701 | 0.8697 | 0.8683 | 0.8683 | 0.8678 |
| GLCTLVAML  | 0.6917 | 0.6883 | 0.699  | 0.6897 | 0.7001 | 0.7038 | 0.713  |
| NLVPMVATV  | 0.6997 | 0.7406 | 0.7669 | 0.7874 | 0.8153 | 0.8248 | 0.8436 |
| TTDPSFLGRY | 0.73   | 0.7518 | 0.7861 | 0.7973 | 0.8165 | 0.8349 | 0.8449 |
| YLQPRTFLL  | 0.922  | 0.9234 | 0.9261 | 0.9272 | 0.928  | 0.9284 | 0.9285 |

## Mean over epitopes

|   fraction |   ('auc_blend', 'thin_both') |   ('auc_blend', 'thin_reference') |   ('auc_density', 'thin_both') |   ('auc_density', 'thin_reference') |
|-----------:|-----------------------------:|----------------------------------:|-------------------------------:|------------------------------------:|
|       0.12 |                       0.8024 |                            0.7975 |                         0.7903 |                              0.7816 |
|       0.2  |                       0.8135 |                            0.8083 |                         0.8011 |                              0.7944 |
|       0.3  |                       0.8308 |                            0.8216 |                         0.8193 |                              0.8097 |
|       0.45 |                       0.8337 |                            0.8256 |                         0.8234 |                              0.8143 |
|       0.6  |                       0.8381 |                            0.8359 |                         0.8288 |                              0.8257 |
|       0.8  |                       0.8445 |                            0.8418 |                         0.8352 |                              0.832  |
|       1    |                       0.8487 |                            0.8487 |                         0.8396 |                              0.8396 |

## Calibration by achieved close-neighbour fraction

| bin          |   n |   mean_close_frac |   auc_density |   auc_blend |   auc_density_sd |
|:-------------|----:|------------------:|--------------:|------------:|-----------------:|
| (0.0, 0.02]  | 159 |            0.0109 |        0.804  |      0.8138 |           0.0715 |
| (0.02, 0.04] | 147 |            0.0291 |        0.8338 |      0.8417 |           0.0746 |
| (0.04, 0.07] |  80 |            0.0539 |        0.8191 |      0.835  |           0.0924 |
| (0.07, 0.12] | 155 |            0.0971 |        0.833  |      0.8463 |           0.0868 |
| (0.12, 0.2]  | 227 |            0.1545 |        0.8382 |      0.8494 |           0.0838 |
| (0.2, 0.35]  | 182 |            0.2763 |        0.7954 |      0.8072 |           0.09   |
| (0.35, 1.01] | 165 |            0.4158 |        0.779  |      0.7939 |           0.0848 |
