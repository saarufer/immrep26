# E31 under family-level resampling

## Independence

| group      |   sites |   pairs |   tcrs |   families |
|:-----------|--------:|--------:|-------:|-----------:|
| control    |     767 |     522 |    302 |         79 |
| discordant |     103 |      71 |     87 |         33 |

Families contributing to both groups: 15.

## Effect with the family as the sampling unit

| quantity                        |   site_effect |   ci_lo |   ci_hi |   family_disc |   family_ctl |   p_family_cell |   paired_mean |   paired_median |   p_paired |   paired_positive |   paired_n |
|:--------------------------------|--------------:|--------:|--------:|--------------:|-------------:|----------------:|--------------:|----------------:|-----------:|------------------:|-----------:|
| peptide contacts at the site    |        3.6292 |  0.4257 |  7.8319 |        9.0432 |       6.1917 |          0.4635 |        0.0874 |               0 |     0.6744 |                 5 |         15 |
| site faces the peptide (binary) |        0.109  | -0.0302 |  0.2704 |        0.4569 |       0.3443 |          0.2486 |        0.123  |               0 |     0.063  |                 6 |         15 |

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
