# E56.1 - provenance audit, the gate on E56-E59

2472 of 2636 competition TCRs (93.8%) match a VDJdb record on (cdr3 alpha, cdr3 beta, epitope), giving 3638 records.

## Field coverage

| field                 |   populated |   distinct_values |
|:----------------------|------------:|------------------:|
| reference.id          |      1      |                57 |
| meta.study.id         |      0.2644 |                23 |
| meta.subject.id       |      0.6611 |               252 |
| meta.replica.id       |      0.0368 |                23 |
| meta.clone.id         |      0.215  |               599 |
| method.identification |      1      |                18 |

## Per epitope

| epitope    | split   |   tcrs_traced |   records |   distinct_references |   largest_reference_share |   subject_populated |   distinct_subjects |   study_populated | subject_resolvable   |
|:-----------|:--------|--------------:|----------:|----------------------:|--------------------------:|--------------------:|--------------------:|------------------:|:---------------------|
| CINGVCWTV  | private |            79 |        79 |                     2 |                    0.6962 |              0.3038 |                   5 |            0      | False                |
| GILGFVFTL  | public  |           670 |      1513 |                    21 |                    0.2022 |              0.9128 |                 158 |            0.3999 | True                 |
| GLCTLVAML  | public  |           192 |       298 |                    14 |                    0.302  |              0.6644 |                  20 |            0.0268 | True                 |
| LLWNGPMAV  | private |           229 |       230 |                     3 |                    0.7783 |              0      |                   0 |            0      | False                |
| LTDEMIAQY  | private |           130 |       130 |                     2 |                    0.9769 |              0      |                   0 |            0      | False                |
| NLVPMVATV  | public  |           366 |       552 |                    24 |                    0.5435 |              0.8931 |                  43 |            0.6322 | True                 |
| TTDPSFLGRY | public  |           387 |       398 |                     2 |                    0.6181 |              0.3543 |                   8 |            0      | False                |
| YLQPRTFLL  | public  |           419 |       438 |                    10 |                    0.4292 |              0.3836 |                  39 |            0      | False                |

## E54's six families, recounted

| family       |   labelled |   records |   distinct_references |   subject_populated |   distinct_subjects |   study_populated |   distinct_studies | subject_resolvable   |
|:-------------|-----------:|----------:|----------------------:|--------------------:|--------------------:|------------------:|-------------------:|:---------------------|
| GILGFVFTL#5  |         50 |       201 |                    13 |              0.9204 |                  67 |            0.4726 |                 15 | True                 |
| GILGFVFTL#1  |         31 |       113 |                    12 |              0.885  |                  50 |            0.4248 |                 10 | True                 |
| GILGFVFTL#27 |         27 |       100 |                    10 |              0.96   |                  52 |            0.58   |                 11 | True                 |
| GILGFVFTL#36 |         27 |       157 |                    14 |              0.9745 |                  68 |            0.4841 |                 11 | True                 |
| GLCTLVAML#2  |         27 |        56 |                    11 |              0.6071 |                  12 |            0      |                  0 | True                 |
| GILGFVFTL#30 |         23 |       284 |                    17 |              0.9613 |                 110 |            0.5246 |                 12 | True                 |

## Gate

- [PASS] most records traceable to a reference
- [PASS] donors resolvable on the public epitopes that hold the families
- [FAIL] literal study-block resampling (E56.2 as written)
- [FAIL] anything study-aware transfers to the private split

The unit that answers the independence question is the donor, not the study identifier; a first version of this gate tested the wrong field.

- E57 on 6 of 6 families: **estimable** with donor as the unit
- E56.2 as literally written (study blocks): **not estimable**, study id 26.4% populated; substitute donor blocks
- any provenance-aware method on the private split: **not estimable**
