# E67/E68 - configuration D's error anatomy

Components rebuilt and verified against the submitted file at r=0.999783.

## E67 - loss by cluster state

| class          |   pairs |    loss |   mean_error_rate |   pair_share |   loss_share |
|:---------------|--------:|--------:|------------------:|-------------:|-------------:|
| clu_clu        |    9251 |  1951.5 |            0.2153 |       0.056  |       0.1015 |
| posclu_negsing |   67871 |  1262.5 |            0.0258 |       0.4107 |       0.0656 |
| possing_negclu |    6493 |  3709   |            0.6342 |       0.0393 |       0.1929 |
| sing_sing      |   81657 | 12308.5 |            0.1916 |       0.4941 |       0.64   |

## E68 gate A - misranked pairs by disagreement quartile

| epitope    |   U_q0 |   U_q1 |   U_q2 |   U_q3 |
|:-----------|-------:|-------:|-------:|-------:|
| GILGFVFTL  |  17.67 |  37.02 |  32.09 |  30.24 |
| GLCTLVAML  |  12.48 |  16.39 |  19.22 |  24.1  |
| NLVPMVATV  |  13.44 |  11.05 |  21.37 |  20.05 |
| TTDPSFLGRY |  19.83 |  29.49 |  32.78 |  21.02 |
| YLQPRTFLL  |  11.26 |  13.56 |  11.82 |  11.56 |

## E68 gate B

U coefficient +1.1484, cluster-robust p 0.0047; gate passes.

## E68 gate C - is U directional?

| side                        |   n |   U_coefficient |   p |
|:----------------------------|----:|----------------:|----:|
| positives (scored too LOW)  | 910 |          1.122  |   0 |
| negatives (scored too HIGH) | 872 |          1.5732 |   0 |

Same sign on both sides: **True**.

## E68 gate D - oracle price of a U correction

| epitope    |   oracle_c |    auc |
|:-----------|-----------:|-------:|
| GILGFVFTL  |       0.1  | 0.8874 |
| GLCTLVAML  |      -0.05 | 0.7578 |
| NLVPMVATV  |       0    | 0.8729 |
| TTDPSFLGRY |       0.15 | 0.8537 |
| YLQPRTFLL  |      -0.15 | 0.9354 |

D as submitted 0.8600; oracle per-epitope 0.8614 (**+0.0014**), against a +0.02 bar.
