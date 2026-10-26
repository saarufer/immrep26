# Pair-conditioned weights: the pairwise formulation never tried

Anchor: the pair decomposition reproduces the frozen entry to 0.000001.

## Label-free classes

| class          |   pair_share |   contribution_at_frozen_weights |   auc_within_class |   loss |   loss_share |
|:---------------|-------------:|---------------------------------:|-------------------:|-------:|-------------:|
| both_singleton |       0.492  |                           0.4023 |             0.8176 | 0.0897 |       0.5976 |
| one_clustered  |       0.4363 |                           0.3942 |             0.9035 | 0.0421 |       0.2805 |
| same_cluster   |       0.003  |                           0.0022 |             0.7336 | 0.0008 |       0.0054 |
| diff_cluster   |       0.0686 |                           0.0511 |             0.7449 | 0.0175 |       0.1165 |

## Oracle, one weight pair per class

| class          |   oracle_w_structure |   oracle_w_typicality |   oracle_w_density |   gain |   within_class_auc |
|:---------------|---------------------:|----------------------:|-------------------:|-------:|-------------------:|
| both_singleton |                 0.2  |                  0    |               0.8  | 0.001  |             0.8196 |
| one_clustered  |                 0.2  |                  0.15 |               0.65 | 0.0016 |             0.9071 |
| same_cluster   |                 0.05 |                  0.1  |               0.85 | 0.0002 |             0.7887 |
| diff_cluster   |                 0.4  |                  0.35 |               0.25 | 0.0012 |             0.7624 |

- frozen entry **0.8498**
- oracle class-conditioned **0.8537** (+0.0039)
- oracle single global weight 0.8515 (structure 0.20, typicality 0.10)
- pair-conditioning buys an oracle +0.0022
- oracle below +0.02, so nothing was fitted

## Mutual distance within both_singleton

| distance_tercile   |   oracle_w_structure |   oracle_w_typicality |   mean_within_band_auc |
|:-------------------|---------------------:|----------------------:|-----------------------:|
| near               |                  0.2 |                     0 |                 0.7634 |
| middle             |                  0.3 |                     0 |                 0.8142 |
| far                |                  0.3 |                     0 |                 0.809  |
