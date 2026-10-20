# Master cross-epitope comparison

Mean per-epitope ROC AUC over the five public epitopes, the competition's
own metric. No method is tuned on the epitope it is scored on; kernel
bandwidth is fixed a priori at 0.1 (SCEPTR cosine) / 20 (tcrdist3), which
the AUC-vs-bandwidth curve shows is flat between 0.05 and 0.15.

Fusions are zero-parameter rank averages with weights fixed in advance,
because fitting a combiner on five epitopes loses more than it gains.

| method                                        |   GILGFVFTL |   GLCTLVAML |   NLVPMVATV |   TTDPSFLGRY |   YLQPRTFLL |   mean |
|:----------------------------------------------|------------:|------------:|------------:|-------------:|------------:|-------:|
| FUSION: 3-space density + structure @0.2      |      0.8732 |      0.7336 |      0.8643 |       0.8504 |      0.9325 | 0.8508 |
| FUSION: 2-space density + structure @0.2      |      0.8733 |      0.739  |      0.8579 |       0.8451 |      0.9326 | 0.8496 |
| FUSION: 3-space density                       |      0.8684 |      0.7147 |      0.8529 |       0.8513 |      0.9276 | 0.843  |
| FUSION: 2-space density                       |      0.8704 |      0.7167 |      0.8458 |       0.8463 |      0.9281 | 0.8415 |
| density: SCEPTR default                       |      0.8678 |      0.713  |      0.8436 |       0.8449 |      0.9285 | 0.8396 |
| density: SCEPTR small                         |      0.8679 |      0.7221 |      0.8387 |       0.824  |      0.9246 | 0.8354 |
| FUSION: 3-space density + specificity @0.2    |      0.8564 |      0.7235 |      0.8476 |       0.8324 |      0.9143 | 0.8348 |
| density: tcrdist3 kernel                      |      0.8593 |      0.6987 |      0.8497 |       0.8367 |      0.9202 | 0.8329 |
| tcrdist3 nndist (organisers' baseline)        |      0.8707 |      0.6901 |      0.8089 |       0.8262 |      0.9278 | 0.8247 |
| density: SCEPTR alpha-chain only              |      0.7973 |      0.7024 |      0.8142 |       0.771  |      0.9079 | 0.7986 |
| structure: TCR-peptide min-PAE (AF3 analogue) |      0.7792 |      0.7362 |      0.7628 |       0.6294 |      0.8746 | 0.7565 |
| density: SCEPTR beta-chain only               |      0.8029 |      0.5852 |      0.7672 |       0.7623 |      0.8352 | 0.7505 |
| structure: whole-complex confidence           |      0.7719 |      0.7349 |      0.7185 |       0.6503 |      0.8578 | 0.7467 |
| metadata: vdjdb_score alone                   |      0.5017 |      0.5811 |      0.399  |       0.5    |      0.4427 | 0.4849 |
| specificity model (target epitope held out)   |      0.5331 |      0.5307 |      0.5528 |       0.5068 |      0.2825 | 0.4812 |
| generation probability (-log10 Pgen paired)   |      0.325  |      0.4851 |      0.4777 |       0.3104 |      0.4847 | 0.4165 |

## Reading it

- Every sequence-neighbourhood density lands in 0.75-0.84 regardless of
  representation. The signal is one-dimensional and saturated.
- The alpha chain beats the beta chain, inverting the usual expectation.
- `vdjdb_score` and the specificity model sit below chance: both are
  quantities whose sign flips between epitopes.
- Structure is a genuine but non-orthogonal signal; it adds ~0.007 when
  blended at low weight and nothing when used alone.
- Selecting any of these choices by cross-validation on five folds costs
  about 0.005, so the configuration is fixed a priori instead.
