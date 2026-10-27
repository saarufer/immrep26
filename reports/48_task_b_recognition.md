# E82 - task B, recognition-profile distance

| representation   | pretrained_on                       |    rho |      p |
|:-----------------|:------------------------------------|-------:|-------:|
| default          | real repertoire                     | 0.0153 | 0.8498 |
| blosum           | real, BLOSUM embedding              | 0.0087 | 0.9144 |
| shuffled_data    | real, alpha-beta pairing randomised | 0.0284 | 0.725  |
| synthetic_data   | OLGA synthetic only                 | 0.0391 | 0.6278 |
| cdr3_only        | real, CDR3 only                     | 0.0345 | 0.669  |
| b_sceptr         | real, beta chain only               | 0.0425 | 0.5982 |
| Levenshtein      | E63's baseline                      | 0.0149 | 0.8537 |

156 pairs, beta-only, non-independent (21 receptors, 16 distinct sequences).
