# A public-maximising second entry, with the bet priced

| configuration                                |   public |   realisable |   pool_dependent |   vs_frozen_public |   vs_frozen_realisable |   realisable_iedb |   vs_frozen_iedb |
|:---------------------------------------------|---------:|-------------:|-----------------:|-------------------:|-----------------------:|------------------:|-----------------:|
| D  density3 + pooled m30 + m50 + structure   |   0.8601 |       0.8561 |           0.0039 |             0.0102 |                 0.0063 |            0.8565 |           0.0067 |
| C  HeidAI rebuild + structure 0.3            |   0.8588 |       0.8534 |           0.0055 |             0.009  |                 0.0035 |            0.854  |           0.0042 |
| A  frozen, density3 -> pooled m30            |   0.8583 |       0.8499 |           0.0084 |             0.0085 |                 0.0001 |            0.8509 |           0.0011 |
| E  (density3 + pooled m30 + m50)/3 + str 0.3 |   0.8582 |       0.8548 |           0.0034 |             0.0083 |                 0.005  |            0.8551 |           0.0052 |
| B  pooled m30 + m50 + structure              |   0.8572 |       0.8521 |           0.0052 |             0.0074 |                 0.0023 |            0.8526 |           0.0028 |
| F  D with no pool at all                     |   0.8557 |     nan      |         nan      |             0.0059 |               nan      |            0.8557 |           0.0059 |
| HeidAI rebuild (m30 pooled + m50)            |   0.8549 |       0.8476 |           0.0073 |             0.0051 |                -0.0022 |            0.8484 |          -0.0015 |
| frozen entry (density3 + structure)          |   0.8498 |       0.8498 |          -0      |             0      |                 0      |            0.8498 |           0      |

## Within both-singleton pairs, which decide the private split

| configuration                                |   both_singleton_auc |   vs_frozen |   worst_epitope_delta |
|:---------------------------------------------|---------------------:|------------:|----------------------:|
| D  density3 + pooled m30 + m50 + structure   |               0.8084 |      0.0143 |                0.0021 |
| C  HeidAI rebuild + structure 0.3            |               0.8083 |      0.0143 |               -0.0004 |
| B  pooled m30 + m50 + structure              |               0.8075 |      0.0134 |               -0.0056 |
| E  (density3 + pooled m30 + m50)/3 + str 0.3 |               0.8056 |      0.0115 |                0.0026 |
| F  D with no pool at all                     |               0.8027 |      0.0087 |                0.0007 |
| A  frozen, density3 -> pooled m30            |               0.8005 |      0.0064 |               -0.0038 |
| HeidAI rebuild (m30 pooled + m50)            |               0.7984 |      0.0043 |               -0.0076 |
| frozen entry (density3 + structure)          |               0.7941 |      0      |                0      |

`pool_dependent` is the share of a candidate's public score the
private epitopes cannot receive at the multipliers E27 counted
(1.038 / 1.045 / 1.834). A second entry wagers that count is low.

Public best: **D  density3 + pooled m30 + m50 + structure** at **0.8601**, against the
frozen entry's 0.8498 and
the rebuilt rival configuration's 0.8549.

Gates on it: private variance True (CINGV:236 LLWNG:220 LTDEM:126); typicality carries something typicality does not.
