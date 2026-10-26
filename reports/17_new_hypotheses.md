# Four new hypotheses

Submitted blend baseline 0.8498; inclusion threshold +0.02 under strict nested LOEO.

H3 is closed by its own stop condition: VDJdb has no nucleotide column, and subject identifiers are absent on two of three private epitopes.

| feature                                    |      auc | private_levels                | viable   |   rho_with_blend |   nested_gain |
|:-------------------------------------------|---------:|:------------------------------|:---------|-----------------:|--------------:|
| alpha-led partner diversity r=0.1          |   0.6465 | CINGV:20 LLWNG:43 LTDEM:16    | True     |           0.5001 |        0      |
| beta-led (control) partner diversity r=0.1 |   0.5981 | CINGV:5 LLWNG:10 LTDEM:5      | True     |           0.3036 |        0      |
| beta-led (control) partner diversity r=0.2 |   0.6637 | CINGV:17 LLWNG:21 LTDEM:10    | True     |           0.4255 |        0      |
| H4 len_deviation (negated)                 |   0.5128 | CINGV:12 LLWNG:9 LTDEM:9      | True     |           0.0123 |        0      |
| H4 hydrophobic_frac (negated)              |   0.4958 | CINGV:75 LLWNG:59 LTDEM:54    | True     |          -0.0213 |        0      |
| alpha-led partner diversity r=0.2          |   0.7038 | CINGV:45 LLWNG:102 LTDEM:21   | True     |           0.5767 |       -0.0003 |
| H4 nglyc_domain (negated)                  |   0.5121 | CINGV:4 LLWNG:4 LTDEM:5       | True     |           0.0002 |       -0.0007 |
| alpha neighbour count r=0.2                |   0.7493 | CINGV:9 LLWNG:21 LTDEM:6      | True     |           0.6938 |       -0.0009 |
| alpha neighbour count r=0.1                |   0.6952 | CINGV:6 LLWNG:8 LTDEM:6       | True     |           0.6117 |       -0.0009 |
| identity (D04 reference)                   |   0.7989 | CINGV:243 LLWNG:229 LTDEM:130 | True     |           0.777  |       -0.0015 |
| BLOSUM62 (D04 reference)                   |   0.7984 | CINGV:243 LLWNG:229 LTDEM:130 | True     |           0.7793 |       -0.0016 |
| Grantham (physicochemical)                 |   0.8022 | CINGV:243 LLWNG:229 LTDEM:130 | True     |           0.7734 |       -0.0022 |
| H4 liability_sum (negated)                 |   0.5155 | CINGV:3 LLWNG:3 LTDEM:2       | True     |           0.0196 |       -0.0025 |
| H4 extra_cys (negated)                     |   0.5106 | CINGV:2 LLWNG:2 LTDEM:2       | True     |           0.0124 |       -0.0047 |
| H4 nglyc_cdr3 (negated)                    |   0.5013 | CINGV:2 LLWNG:1 LTDEM:2       | False    |          -0.0056 |      nan      |
| H4 noncanonical_ends (negated)             | nan      | CINGV:1 LLWNG:1 LTDEM:1       | False    |         nan      |      nan      |
