# E31 under family-level resampling

## Independence

| group      |   sites |   pairs |   tcrs |   families |
|:-----------|--------:|--------:|-------:|-----------:|
| control    |     767 |     522 |    302 |         79 |
| discordant |     103 |      71 |     87 |         33 |

Families contributing to both groups: 15.

## Effect with the family as the sampling unit

| quantity                        |   discordant |   control |   effect |   ci_lo |   ci_hi |   p_family_perm | crosses_zero   |
|:--------------------------------|-------------:|----------:|---------:|--------:|--------:|----------------:|:---------------|
| peptide contacts at the site    |       6.534  |    2.9048 |   3.6292 |  0.4257 |  7.8319 |          0.043  | False          |
| site faces the peptide (binary) |       0.3398 |    0.2308 |   0.109  | -0.029  |  0.2745 |          0.2357 | True           |

## Leave-one-epitope-out

| dropped    |   sites_lost |   contacts_effect |   facing_effect |
|:-----------|-------------:|------------------:|----------------:|
| GILGFVFTL  |          641 |             4.6   |           0.16  |
| GLCTLVAML  |           96 |             3.248 |           0.12  |
| NLVPMVATV  |           73 |             2.352 |           0.026 |
| TTDPSFLGRY |            3 |             3.618 |           0.108 |
| YLQPRTFLL  |           57 |             4.062 |           0.123 |
| (none)     |            0 |             3.629 |           0.109 |

## Per epitope

|                             |   sites |   families |   contacts |   facing |
|:----------------------------|--------:|-----------:|-----------:|---------:|
| ('GILGFVFTL', 'control')    |     576 |         42 |      1.545 |    0.177 |
| ('GILGFVFTL', 'discordant') |      65 |         17 |      3.569 |    0.215 |
| ('GLCTLVAML', 'control')    |      86 |          7 |      7.012 |    0.453 |
| ('GLCTLVAML', 'discordant') |      10 |          6 |     14.9   |    0.5   |
| ('NLVPMVATV', 'control')    |      47 |         13 |      7.638 |    0.383 |
| ('NLVPMVATV', 'discordant') |      26 |          9 |     11.231 |    0.615 |
| ('TTDPSFLGRY', 'control')   |       3 |          2 |      0     |    0     |
| ('YLQPRTFLL', 'control')    |      55 |         15 |      6.836 |    0.327 |
| ('YLQPRTFLL', 'discordant') |       2 |          1 |      0     |    0     |
