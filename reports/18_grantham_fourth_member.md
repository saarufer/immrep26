# Grantham as a fourth density member

## Members, standalone

|                |   mean AUC |
|:---------------|-----------:|
| SCEPTR default |     0.8396 |
| SCEPTR small   |     0.8354 |
| tcrdist3       |     0.8329 |
| Grantham CDR3  |     0.8022 |
| V-domain align |     0.8354 |

## Within-epitope Spearman

|                |   SCEPTR default |   SCEPTR small |   tcrdist3 |   Grantham CDR3 |   V-domain align |
|:---------------|-----------------:|---------------:|-----------:|----------------:|-----------------:|
| SCEPTR default |            1     |          0.907 |      0.878 |           0.779 |            0.831 |
| SCEPTR small   |            0.907 |          1     |      0.829 |           0.737 |            0.838 |
| tcrdist3       |            0.878 |          0.829 |      1     |           0.836 |            0.818 |
| Grantham CDR3  |            0.779 |          0.737 |      0.836 |           1     |            0.664 |
| V-domain align |            0.831 |          0.838 |      0.818 |           0.664 |            1     |

## Configurations

| configuration              |   density only |   with structure |   GILGFVFTL |   GLCTLVAML |   NLVPMVATV |   TTDPSFLGRY |   YLQPRTFLL |   vs submitted |
|:---------------------------|---------------:|-----------------:|------------:|------------:|------------:|-------------:|------------:|---------------:|
| 3 members (submitted)      |         0.843  |           0.8498 |      0.8693 |      0.7414 |      0.8649 |       0.8402 |      0.9334 |         0      |
| 4 members + Grantham       |         0.8405 |           0.8475 |      0.8722 |      0.736  |      0.8608 |       0.8355 |      0.9328 |        -0.0024 |
| 4 members + V-domain (E26) |         0.846  |           0.8524 |      0.869  |      0.7464 |      0.8711 |       0.8406 |      0.9347 |         0.0025 |
| 5 members                  |         0.8449 |           0.8505 |      0.8714 |      0.7388 |      0.8669 |       0.8406 |      0.935  |         0.0007 |
| Grantham replaces tcrdist3 |         0.8387 |           0.8469 |      0.8727 |      0.7371 |      0.856  |       0.8364 |      0.9324 |        -0.0029 |

## Nested LOEO over member choice

0.8524

A priori 4-member+Grantham: 0.8475, delta -0.0024 against the submitted 0.8498. does not clear +0.02; the entry stands.
