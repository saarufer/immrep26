# Pool expansion: self-hit and duplication audit

Plan card Q02. Every test TCR came from VDJdb, so expanding the cohort with
VDJdb risks counting the query itself, or one rearrangement curated several
times, as independent neighbours.

## Self-hit rate and pool composition

| epitope    |   n_query |   pool_records |   pool_distinct_clonotypes |   query_with_exact_selfhit |   mean_exact_selfhits |   mean_clonotype_copies |   pool_refs |   pool_top_ref_share |   new_records_per_query |   new_clonotypes_per_query |
|:-----------|----------:|---------------:|---------------------------:|---------------------------:|----------------------:|------------------------:|------------:|---------------------:|------------------------:|---------------------------:|
| CINGVCWTV  |       243 |            254 |                        254 |                      0.942 |                 0.942 |                   1     |           3 |                0.625 |                   1.045 |                      1.045 |
| GILGFVFTL  |       670 |           6062 |                       3674 |                      0.99  |                 2.193 |                   2.258 |          23 |                0.349 |                   9.048 |                      5.484 |
| GLCTLVAML  |       192 |           1030 |                        650 |                      0.958 |                 1.151 |                   1.552 |          16 |                0.508 |                   5.365 |                      3.385 |
| LLWNGPMAV  |       229 |            421 |                        420 |                      0.878 |                 0.882 |                   1.004 |           4 |                0.468 |                   1.838 |                      1.834 |
| LTDEMIAQY  |       130 |            135 |                        135 |                      1     |                 1     |                   1     |           2 |                0.978 |                   1.038 |                      1.038 |
| NLVPMVATV  |       366 |          12793 |                      10245 |                      0.945 |                 1.383 |                   1.508 |          30 |                0.948 |                  34.954 |                     27.992 |
| TTDPSFLGRY |       387 |            424 |                        412 |                      1     |                 1.021 |                   1.028 |           3 |                0.599 |                   1.096 |                      1.065 |
| YLQPRTFLL  |       419 |            492 |                        472 |                      0.998 |                 1.038 |                   1.045 |          12 |                0.413 |                   1.174 |                      1.126 |

## Records vs independent clonotypes vs studies

| epitope    |   records/query |   clonotypes/query |   record:clonotype |   distinct_refs |   top_ref_share |
|:-----------|----------------:|-------------------:|-------------------:|----------------:|----------------:|
| CINGVCWTV  |           1.045 |              1.045 |              1     |               3 |           0.625 |
| GILGFVFTL  |           9.048 |              5.484 |              1.65  |              23 |           0.349 |
| GLCTLVAML  |           5.365 |              3.385 |              1.585 |              16 |           0.508 |
| LLWNGPMAV  |           1.838 |              1.834 |              1.002 |               4 |           0.468 |
| LTDEMIAQY  |           1.038 |              1.038 |              1     |               2 |           0.978 |
| NLVPMVATV  |          34.954 |             27.992 |              1.249 |              30 |           0.948 |
| TTDPSFLGRY |           1.096 |              1.065 |              1.029 |               3 |           0.599 |
| YLQPRTFLL  |           1.174 |              1.126 |              1.042 |              12 |           0.413 |
