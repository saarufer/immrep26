# E64 - E61's calibration under family holdout

## Is sequence-family leakage real here?

| tcr   |   measured |   modal_length |   at_modal_length |   hamming1_pairs |   distinct_genes |   distinct_origins |
|:------|-----------:|---------------:|------------------:|-----------------:|-----------------:|-------------------:|
| TCR59 |        192 |              9 |                68 |                1 |              160 |                  2 |
| TCR63 |        192 |              9 |                80 |                4 |              157 |                  2 |

Essentially no near-neighbour structure, so E61's plain KFold was not leaking. Gene and origin are the groupings that exist.

## Held-out Spearman by validation scheme and model

|                                                  |   gene holdout |   origin holdout (human<->virus) |   plain KFold (E61) |
|:-------------------------------------------------|---------------:|---------------------------------:|--------------------:|
| ('TCR59', 'isotonic, direction from data')       |         0.6393 |                           0.665  |              0.6348 |
| ('TCR59', 'linear (log x)')                      |         0.6835 |                           0.6798 |              0.6824 |
| ('TCR59', 'saturating sigmoid')                  |         0.5881 |                           0.6568 |              0.5776 |
| ('TCR59', 'two-piece linear, knot bounded in x') |         0.658  |                           0.6406 |              0.6485 |
| ('TCR63', 'isotonic, direction from data')       |         0.4237 |                           0.5135 |              0.4271 |
| ('TCR63', 'linear (log x)')                      |         0.3091 |                           0.3185 |              0.3199 |
| ('TCR63', 'saturating sigmoid')                  |         0.4048 |                           0.3676 |              0.3863 |
| ('TCR63', 'two-piece linear, knot bounded in x') |         0.5812 |                           0.6323 |              0.5671 |
