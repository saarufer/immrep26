# E54 - is family purity more than the continuous density relation?

## The graph (full cohort) against the purity set (labelled only)

| epitope    |   cohort_rows |   labelled |   unlabelled_in_graph |   families_in_graph |   families_ge3_labelled |   families_ge2_labelled |
|:-----------|--------------:|-----------:|----------------------:|--------------------:|------------------------:|------------------------:|
| CINGVCWTV  |           243 |          0 |                   243 |                 240 |                       0 |                       0 |
| GILGFVFTL  |           670 |        575 |                    95 |                 390 |                      21 |                      35 |
| GLCTLVAML  |           192 |        166 |                    26 |                 108 |                       8 |                      17 |
| LLWNGPMAV  |           229 |          0 |                   229 |                 223 |                       0 |                       0 |
| LTDEMIAQY  |           130 |          0 |                   130 |                 118 |                       0 |                       0 |
| NLVPMVATV  |           366 |        317 |                    49 |                 320 |                       4 |                      20 |
| TTDPSFLGRY |           387 |        342 |                    45 |                 375 |                       2 |                       5 |
| YLQPRTFLL  |           419 |        382 |                    37 |                 361 |                      10 |                      19 |

Families with >=3 labelled members: **45**, holding 390 TCRs. Observed: pure-positive **26**, pure-negative 3, mixed 16, pair concordance **0.9092**.

## The density-to-label relation the null preserves

|   bin10 |   n |   pos_rate |   nn3_mean |
|--------:|----:|-----------:|-----------:|
|       0 |  41 |     0.9268 |     0.026  |
|       1 |  40 |     0.925  |     0.0393 |
|       2 |  38 |     0.8947 |     0.0476 |
|       3 |  39 |     0.9487 |     0.0572 |
|       4 |  39 |     0.8462 |     0.0665 |
|       5 |  39 |     0.9487 |     0.0768 |
|       6 |  38 |     0.8421 |     0.0878 |
|       7 |  38 |     0.9474 |     0.0995 |
|       8 |  38 |     0.9737 |     0.1132 |
|       9 |  40 |     0.825  |     0.1429 |

## Permutation within (epitope x nn3 bin)

|   n_bins |   n_strata |   movable_rows |   movable_frac |   pure_pos_obs |   pure_pos_null |   pure_pos_p |   pure_pos_z |   pure_neg_obs |   pure_neg_null |   pure_neg_p |   pure_neg_z |   mixed_obs |   mixed_null |   mixed_p |   mixed_z |   concordance_obs |   concordance_null |   concordance_p |   concordance_z |
|---------:|-----------:|---------------:|---------------:|---------------:|----------------:|-------------:|-------------:|---------------:|----------------:|-------------:|-------------:|------------:|-------------:|----------:|----------:|------------------:|-------------------:|----------------:|----------------:|
|        1 |          5 |            383 |         0.9821 |             26 |         26.4785 |       0.7066 |      -0.2436 |              3 |          0.026  |       0.0005 |      18.3301 |          16 |      18.4955 |    0.1539 |   -1.2635 |            0.9092 |             0.829  |          0.0005 |          4.6791 |
|        3 |         15 |            353 |         0.9051 |             26 |         26.6575 |       0.7416 |      -0.3613 |              3 |          0.012  |       0.0005 |      27.4417 |          16 |      18.3305 |    0.1659 |   -1.2753 |            0.9092 |             0.8272 |          0.0005 |          5.0288 |
|        5 |         25 |            318 |         0.8154 |             26 |         25.884  |       0.5767 |       0.0613 |              3 |          0.0135 |       0.0005 |      25.879  |          16 |      19.1025 |    0.0885 |   -1.6293 |            0.9092 |             0.8361 |          0.0005 |          4.3498 |
|       10 |         47 |            252 |         0.6462 |             26 |         24.465  |       0.2699 |       0.8791 |              3 |          0.021  |       0.0005 |      20.7764 |          16 |      20.514  |    0.0135 |   -2.5569 |            0.9092 |             0.8541 |          0.0005 |          3.5482 |
|       20 |         87 |            176 |         0.4513 |             26 |         24.581  |       0.2889 |       0.8646 |              3 |          0.0525 |       0.0005 |      13.0847 |          16 |      20.3665 |    0.009  |   -2.5928 |            0.9092 |             0.8609 |          0.002  |          3.1118 |

## Where the extra agreement sits, by family size

| family_size   |   families |   tcrs |   pairs |   concordance |   null |   excess |      p |
|:--------------|-----------:|-------:|--------:|--------------:|-------:|---------:|-------:|
| 3-4           |         23 |     75 |      87 |        0.8621 | 0.8103 |   0.0518 | 0.1617 |
| 5-9           |         12 |     73 |     190 |        0.7947 | 0.8087 |  -0.0139 | 0.6487 |
| 10-19         |          4 |     57 |     393 |        0.9237 | 0.8757 |   0.048  | 0.1956 |
| 20-+          |          6 |    185 |    2996 |        0.9159 | 0.8566 |   0.0593 | 0.006  |

Wider permutation pool (all labelled rows in the stratum): concordance 0.9092 against 0.8312, p 0.0005.


## The families carrying the excess

| family       |   labelled |   validated |   nn3_median |   distinct_references | shared_by_all   | dominant_TRBV    |
|:-------------|-----------:|------------:|-------------:|----------------------:|:----------------|:-----------------|
| GILGFVFTL#5  |         50 |       0.92  |       0.0636 |                    13 | none            | TRBV19 (50/50)   |
| GILGFVFTL#1  |         31 |       1     |       0.0664 |                    12 | none            | TRBV19 (31/31)   |
| GILGFVFTL#36 |         27 |       1     |       0.0353 |                    13 | none            | TRBV19 (27/27)   |
| GILGFVFTL#27 |         27 |       1     |       0.0482 |                    10 | none            | TRBV19 (27/27)   |
| GLCTLVAML#2  |         27 |       0.963 |       0.0606 |                    11 | none            | TRBV20-1 (27/27) |
| GILGFVFTL#30 |         23 |       0.913 |       0.033  |                    19 | none            | TRBV19 (23/23)   |

## Provenance

| family         |   n |   purity | all_share_a_reference   |   n_distinct_references |
|:---------------|----:|---------:|:------------------------|------------------------:|
| GILGFVFTL#1    |  31 |   1      | False                   |                      12 |
| GILGFVFTL#127  |   3 |   0      | True                    |                       1 |
| GILGFVFTL#139  |   3 |   0.6667 | False                   |                       7 |
| GILGFVFTL#16   |   3 |   1      | False                   |                       2 |
| GILGFVFTL#174  |   3 |   1      | True                    |                       3 |
| GILGFVFTL#2    |   7 |   1      | False                   |                       6 |
| GILGFVFTL#220  |   6 |   0.6667 | False                   |                       7 |
| GILGFVFTL#228  |   3 |   1      | False                   |                       2 |
| GILGFVFTL#27   |  27 |   1      | False                   |                      10 |
| GILGFVFTL#30   |  23 |   0.913  | False                   |                      19 |
| GILGFVFTL#31   |   4 |   1      | False                   |                       2 |
| GILGFVFTL#34   |   3 |   1      | False                   |                       2 |
| GILGFVFTL#36   |  27 |   1      | False                   |                      13 |
| GILGFVFTL#4    |   8 |   0.75   | False                   |                       9 |
| GILGFVFTL#5    |  50 |   0.92   | False                   |                      13 |
| GILGFVFTL#52   |   4 |   1      | False                   |                       4 |
| GILGFVFTL#54   |   5 |   1      | False                   |                       3 |
| GILGFVFTL#56   |   5 |   1      | False                   |                       7 |
| GILGFVFTL#62   |   4 |   0      | True                    |                       1 |
| GILGFVFTL#77   |   3 |   1      | False                   |                       3 |
| GILGFVFTL#8    |   6 |   0.8333 | False                   |                       5 |
| GLCTLVAML#2    |  27 |   0.963  | False                   |                      11 |
| GLCTLVAML#20   |   5 |   0.8    | False                   |                       2 |
| GLCTLVAML#3    |  17 |   0.8824 | False                   |                       9 |
| GLCTLVAML#32   |   3 |   0.6667 | True                    |                       1 |
| GLCTLVAML#33   |   3 |   1      | False                   |                       4 |
| GLCTLVAML#36   |   6 |   0      | True                    |                       1 |
| GLCTLVAML#49   |   3 |   0.6667 | True                    |                       1 |
| GLCTLVAML#54   |   3 |   0.6667 | False                   |                       2 |
| NLVPMVATV#16   |  17 |   1      | False                   |                       9 |
| NLVPMVATV#17   |   6 |   0.8333 | False                   |                       8 |
| NLVPMVATV#20   |   6 |   1      | False                   |                       8 |
| NLVPMVATV#98   |   3 |   0.6667 | False                   |                       3 |
| TTDPSFLGRY#156 |   4 |   1      | False                   |                       2 |
| TTDPSFLGRY#62  |   3 |   1      | True                    |                       1 |
| YLQPRTFLL#11   |  11 |   1      | False                   |                       6 |
| YLQPRTFLL#126  |   4 |   1      | False                   |                       2 |
| YLQPRTFLL#195  |   3 |   1      | False                   |                       2 |
| YLQPRTFLL#22   |   3 |   1      | False                   |                       3 |
| YLQPRTFLL#27   |   7 |   1      | False                   |                       2 |
| YLQPRTFLL#35   |   3 |   0.6667 | False                   |                       2 |
| YLQPRTFLL#50   |   6 |   0.8333 | False                   |                       4 |
| YLQPRTFLL#67   |   4 |   1      | False                   |                       3 |
| YLQPRTFLL#7    |  12 |   1      | False                   |                       7 |
| YLQPRTFLL#72   |   3 |   1      | False                   |                       2 |
