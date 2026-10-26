# E58 - generation-adjusted motif enrichment

## Part 1 - generation probability by group (GILGFVFTL)

| group                    |   n |   median_log10_pgen_beta |   median_cdr3b_length |
|:-------------------------|----:|-------------------------:|----------------------:|
| in a large TRBV19 family | 205 |                  -7.0226 |                    13 |
| TRBV19, no large family  | 318 |                  -7.5774 |                    13 |
| other V, no large family | 147 |                  -9.1222 |                    15 |

## Part 2 - motif enrichment over an OLGA background

| family       |   members |   L | consensus     |   radius |   observed_in_cohort |   olga_background |   enrichment |
|:-------------|----------:|----:|:--------------|---------:|---------------------:|------------------:|-------------:|
| GILGFVFTL#5  |        71 |  13 | CASSIRSTGEQYF |        3 |             0.408955 |           0.00055 |      743.555 |
| GILGFVFTL#1  |        31 |  12 | CASSIGSYGYTF  |        1 |             0.041791 |           4e-05   |     1044.78  |
| GILGFVFTL#27 |        29 |  13 | CASSIRSTDTQYF |        1 |             0.080597 |           2e-05   |     4029.85  |
| GLCTLVAML#2  |        32 |  13 | CSARDGTGNGYTF |        3 |             0.125    |           7.5e-05 |     1666.67  |
| GILGFVFTL#36 |        37 |  13 | CASSIRSAYEQYF |        1 |             0.159701 |           3.5e-05 |     4562.9   |
| GILGFVFTL#30 |        37 |  13 | CASSIRSSYEQYF |        1 |             0.20597  |           5.5e-05 |     3744.91  |

## Part 3 - against the matched control

| family       |   enrichment |   control_median |   control_p95 |   p_vs_control |
|:-------------|-------------:|-----------------:|--------------:|---------------:|
| GILGFVFTL#5  |      743.555 |          34.0441 |       78.6295 |          0.005 |
| GILGFVFTL#1  |     1044.78  |          38.2325 |      369.244  |          0.005 |
| GILGFVFTL#27 |     4029.85  |          37.1734 |      313.968  |          0.005 |
| GLCTLVAML#2  |     1666.67  |         nan      |      nan      |        nan     |
| GILGFVFTL#36 |     4562.9   |          34.0441 |      129.742  |          0.005 |
| GILGFVFTL#30 |     3744.91  |          34.0441 |      129.742  |          0.005 |

Motifs are defined on CDR3beta only, because OLGA generates no V gene; the families' TRBV19 restriction is reported separately and is not folded into the enrichment.
