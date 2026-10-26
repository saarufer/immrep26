# A private estimate anchored on the published baseline

## Does the paper's tcrdist3 match ours?

| epitope    |   published_lo |   published_hi |   measured | inside   |    gap |
|:-----------|---------------:|---------------:|-----------:|:---------|-------:|
| GILGFVFTL  |           0.87 |           0.93 |     0.8707 | True     | 0      |
| GLCTLVAML  |           0.66 |           0.68 |     0.6901 | False    | 0.0101 |
| NLVPMVATV  |           0.8  |           0.83 |     0.8089 | True     | 0      |
| TTDPSFLGRY |           0.8  |           0.83 |     0.8262 | True     | 0      |
| YLQPRTFLL  |           0.87 |           0.93 |     0.9278 | True     | 0      |

## Shift and estimate

- tcrdist3 public macro, measured: **0.8247**
- tcrdist3 private macro, published: **0.753-0.780**
- shift: **-0.071 to -0.045**
- our public: 0.8498, margin +0.0251
- **private estimate 0.778-0.805**, midpoint 0.792

Source: Immunity 2026, doi:10.1016/j.immuni.2026.09.007, Figure 4 results text. Published literature, independent of any other team.
