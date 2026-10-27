# E75 - structural comparison of GLCTLVAML against the other epitopes

## Docking geometry (A*02:01 9-mers)

| epitope   |   n |   crossing_angle_deg_median |   crossing_angle_deg_sd |   tcr_height_over_peptide_median |   tcr_height_over_peptide_sd |   contacts_a_pep_median |   contacts_a_pep_sd |   contacts_b_pep_median |   contacts_b_pep_sd |   cdr3_share_a_median |   cdr3_share_a_sd |   cdr3_share_b_median |   cdr3_share_b_sd |
|:----------|----:|----------------------------:|------------------------:|---------------------------------:|-----------------------------:|------------------------:|--------------------:|------------------------:|--------------------:|----------------------:|------------------:|----------------------:|------------------:|
| CINGVCWTV | 243 |                     75.0989 |                 12.375  |                          21.7587 |                       0.8819 |                    37   |             20.1201 |                    55   |             20.7671 |                0.6508 |            0.3103 |                1      |            0.1151 |
| GILGFVFTL | 670 |                     81.0828 |                 11.1033 |                          21.8112 |                       0.6676 |                    32.5 |             13.8188 |                    51   |             13.2988 |                1      |            0.2048 |                0.4651 |            0.2425 |
| GLCTLVAML | 192 |                     69.4844 |                 15.1335 |                          20.8107 |                       0.9628 |                    32   |             15.8741 |                    68.5 |             21.004  |                0.6088 |            0.2418 |                0.8948 |            0.1463 |
| LLWNGPMAV | 229 |                     64.9576 |                  9.7926 |                          21.0111 |                       0.8462 |                    51   |             13.6535 |                    51   |             17.1509 |                0.3103 |            0.2696 |                0.9846 |            0.0902 |
| NLVPMVATV | 366 |                     65.5757 |                 14.9377 |                          21.6113 |                       1.0092 |                    36   |             18.7081 |                    56   |             24.1704 |                0.7878 |            0.258  |                0.8794 |            0.2164 |
| YLQPRTFLL | 419 |                     73.7842 |                 15.7457 |                          21.2613 |                       1.2537 |                    71   |             21.9819 |                    57   |             20.2041 |                0.6471 |            0.2024 |                0.9487 |            0.1922 |

## Contact footprint

| epitope   | hottest   |   top3_share |   effective_positions |
|:----------|:----------|-------------:|----------------------:|
| CINGVCWTV | P7        |       0.8233 |                4.6204 |
| GILGFVFTL | P4        |       0.8096 |                4.4855 |
| GLCTLVAML | P4        |       0.8275 |                4.4777 |
| LLWNGPMAV | P4        |       0.7766 |                4.9552 |
| NLVPMVATV | P4        |       0.7835 |                4.9106 |
| YLQPRTFLL | P4        |       0.7551 |                5.3356 |

## Ensemble tightness against Boltz confidence

| epitope    |   n |   tightness |   iptm |   residual |
|:-----------|----:|------------:|-------:|-----------:|
| CINGVCWTV  | 243 |      0.1416 | 0.8834 |    -0.0214 |
| GILGFVFTL  | 670 |      0.0856 | 0.9326 |    -0.0105 |
| GLCTLVAML  | 192 |      0.1076 | 0.9241 |     0      |
| LLWNGPMAV  | 229 |      0.0783 | 0.9472 |     0.0021 |
| LTDEMIAQY  | 130 |      0.1451 | 0.9096 |     0.0177 |
| NLVPMVATV  | 366 |      0.1497 | 0.8913 |    -0.0024 |
| TTDPSFLGRY | 387 |      0.1799 | 0.8851 |     0.0193 |
| YLQPRTFLL  | 419 |      0.1319 | 0.9028 |    -0.0047 |

Spearman(ipTM, tightness) = **-0.810**; footprint width against sequence convergence rho -0.486, exact p 0.356.
