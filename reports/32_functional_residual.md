# E59 - functional residual of the public motifs (GILGFVFTL only)

GLCTLVAML is **not estimable** (E56: 13 matched controls).

Raw: motif 0.962 (n=158) against matched TRBV19 non-motif 0.745 (n=278), difference **+0.217**, odds ratio 8.7, Fisher p 6.43e-10.

## Matching ladder

| matching                                      |   difference |   motif_tcrs_covered |   strata_used |
|:----------------------------------------------|-------------:|---------------------:|--------------:|
| CDR3b length                                  |       0.2098 |                  158 |             2 |
| length x Pgen tercile                         |       0.2296 |                  158 |             6 |
| length x Pgen x nn3 tercile (density as well) |       0.0824 |                  137 |             9 |

The density-matched row controls a **mediator** and is not the causal estimate.

## Leave-one-out sensitivity

| level   | removed            |   difference |   n_motif_left |
|:--------|:-------------------|-------------:|---------------:|
| family  | GILGFVFTL#1        |       0.2082 |            127 |
| family  | GILGFVFTL#27       |       0.2096 |            131 |
| family  | GILGFVFTL#30       |       0.2258 |            135 |
| family  | GILGFVFTL#36       |       0.2096 |            131 |
| family  | GILGFVFTL#5        |       0.2369 |            108 |
| motif   | other:CASSIGSYGYTF |       0.2082 |            127 |
| motif   | CASSIRS-13         |       0.2554 |             31 |

Motif groups: other:CASSIGSYGYTF = ['GILGFVFTL#1']; CASSIRS-13 = ['GILGFVFTL#27', 'GILGFVFTL#30', 'GILGFVFTL#36', 'GILGFVFTL#5']
