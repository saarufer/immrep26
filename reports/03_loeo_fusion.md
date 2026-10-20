# LOEO: density vs cluster shape vs fusion

Mean per-epitope AUC, leave-one-epitope-out. Columns are the held-out epitope.

|                                   |   GILGFVFTL |   GLCTLVAML |   NLVPMVATV |   TTDPSFLGRY |   YLQPRTFLL |   mean |
|:----------------------------------|------------:|------------:|------------:|-------------:|------------:|-------:|
| [single] default__kern0.1         |      0.8708 |      0.7018 |      0.8384 |       0.842  |      0.9279 | 0.8362 |
| [logreg C=0.1] density only       |      0.8698 |      0.7081 |      0.834  |       0.8394 |      0.9283 | 0.8359 |
| [logreg C=0.1] density + cluster  |      0.8673 |      0.7102 |      0.8143 |       0.8308 |      0.9272 | 0.8299 |
| [hgb depth=2] density + cluster   |      0.8381 |      0.6821 |      0.8425 |       0.8288 |      0.9023 | 0.8188 |
| [single] a_sceptr__ntopk15        |      0.8249 |      0.7248 |      0.8077 |       0.7591 |      0.9167 | 0.8067 |
| [logreg C=0.1] cluster shape only |      0.7405 |      0.7002 |      0.7613 |       0.5294 |      0.6914 | 0.6846 |
| [single] clu_entropy_b            |      0.5581 |      0.6573 |      0.5216 |       0.5059 |      0.5335 | 0.5553 |
