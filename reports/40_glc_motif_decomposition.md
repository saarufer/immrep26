# E69 - GLCTLVAML motif decomposition

Groups fixed before any error was examined. Cohort 192, labelled 166.

## Pre-fixed groups

| group                     |   n |   validated |   mean_D_rank |   AUC_within_group |   misranked_share |
|:--------------------------|----:|------------:|--------------:|-------------------:|------------------:|
| public motif family       |  27 |      0.963  |        0.8343 |             0.9231 |            0.03   |
| TRAV5 only                |  46 |      0.8261 |        0.6421 |             0.6513 |            0.233  |
| neither                   |  80 |      0.5125 |        0.3147 |             0.6266 |            0.6121 |
| TRAV5+TRBV20-1, no family |   8 |      0.5    |        0.5692 |           nan      |            0.0871 |
| TRBV20-1 only             |   5 |      0.2    |        0.3083 |           nan      |            0.0377 |

## Residual-error test, conditioning on D's own score

| term          |   coefficient |      p |   pseudo_r2 |
|:--------------|--------------:|-------:|------------:|
| TRAV5         |        0.1866 | 0.3759 |      0.1708 |
| TRBV20_1      |       -0.0436 | 0.8961 |      0.1677 |
| in_big_family |       -1.3689 | 0      |      0.2288 |

Surviving terms: ['in_big_family'].

## Does the surviving term exist on the private split?

| epitope    | split   |   in_big_family |
|:-----------|:--------|----------------:|
| CINGVCWTV  | private |          0      |
| GILGFVFTL  | public  |          0.306  |
| GLCTLVAML  | public  |          0.1667 |
| LLWNGPMAV  | private |          0      |
| LTDEMIAQY  | private |          0      |
| NLVPMVATV  | public  |          0      |
| TTDPSFLGRY | public  |          0      |
| YLQPRTFLL  | public  |          0      |

## Direction, by label

| side                        |   n |   in_family |   coefficient |   p |
|:----------------------------|----:|------------:|--------------:|----:|
| positives (scored too LOW)  | 110 |          26 |       -1.2879 |   0 |
| negatives (scored too HIGH) |  56 |           1 |      nan      | nan |
