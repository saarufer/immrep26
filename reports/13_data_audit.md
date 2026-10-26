# Data-layer audit

57 checks across the artifacts behind all 39 experiments: 2 NOTE, 55 PASS.

| group       | check                                                                | status   | detail                                                   |
|:------------|:---------------------------------------------------------------------|:---------|:---------------------------------------------------------|
| core        | test.csv row count                                                   | PASS     | 2636 vs expected 2636                                    |
| core        | tcr_id unique                                                        | PASS     | 2636 distinct                                            |
| core        | per-epitope counts match the data description                        | PASS     | all 8 match                                              |
| core        | 8 epitopes, 5 public and 3 private                                   | PASS     | 8 epitopes                                               |
| core        | no missing sequence or HLA fields                                    | PASS     | complete                                                 |
| core        | split sizes                                                          | PASS     | {'private': 602, 'public': 2034}                         |
| labels      | scored rows = public split minus ambiguous                           | PASS     | 1782 + 252 ambiguous = 2034                              |
| labels      | labelled ids are a subset of the public split                        | PASS     | 1782 ids                                                 |
| labels      | no private row carries a label                                       | PASS     | private split is unlabelled, as the competition requires |
| labels      | labels are binary                                                    | PASS     | positives 910, rate 0.511                                |
| labels      | every public epitope has both classes                                | PASS     | all 5 mixed                                              |
| metric      | tcrdist3 baseline reproduces the organisers' number                  | PASS     | 0.8247 vs published 0.824                                |
| metric      | density3_struct03.csv reproduces its leaderboard score               | PASS     | 0.8498 truncates to 0.849, leaderboard showed 0.849      |
| metric      | density_sceptr_default_sigma0.1.csv reproduces its leaderboard score | PASS     | 0.8396 truncates to 0.839, leaderboard showed 0.839      |
| sequence    | pgen.csv: ids all present in test.csv                                | PASS     | 2636 distinct                                            |
| sequence    | pgen.csv: covers every test row                                      | PASS     | 2636/2636                                                |
| sequence    | pgen.csv: all numeric values finite                                  | PASS     | finite                                                   |
| sequence    | specificity_scores.csv: ids all present in test.csv                  | PASS     | 2636 distinct                                            |
| sequence    | specificity_scores.csv: all numeric values finite                    | PASS     | finite                                                   |
| sequence    | sceptr/ distance or embedding cache                                  | PASS     | 7 files                                                  |
| sequence    | tcrdist/ distance or embedding cache                                 | PASS     | 8 files                                                  |
| structure   | boltz_features.csv: ids all present in test.csv                      | PASS     | 2636 distinct                                            |
| structure   | boltz_features.csv covers every test row                             | PASS     | 2636/2636                                                |
| structure   | boltz_features.csv has no duplicate ids                              | PASS     | 2636 rows                                                |
| structure   | min-PAE sign convention is uniform                                   | PASS     | all non-positive, i.e. negated so higher is better       |
| structure   | min-PAE magnitude in a physical range                                | PASS     | 0.32-14.80 A                                             |
| structure   | ipTM within [0, 1]                                                   | PASS     | 0.544-0.981                                              |
| structure   | folded complexes match the test set                                  | PASS     | 2636/2636                                                |
| structure   | boltz_consistency.csv: ids all present in test.csv                   | PASS     | 2636 distinct                                            |
| structure   | structure_contacts.csv: ids all present in test.csv                  | PASS     | 2636 distinct                                            |
| structure   | minimisation_robustness.csv: ids all present in test.csv             | PASS     | 69 distinct                                              |
| structure   | contact_fingerprints.npz loads                                       | PASS     | 2 arrays                                                 |
| structure   | decoy manifest: ids all present in test.csv                          | PASS     | 902 distinct                                             |
| structure   | decoy folds complete                                                 | PASS     | 2706/2706 target-decoy complexes                         |
| structure   | both decoy phases ran                                                | PASS     | {1: 300, 2: 602}                                         |
| structure   | relaxed structures present                                           | PASS     | 69 PDBs                                                  |
| structure   | every relaxed structure has 4 chains                                 | PASS     | all 4-chain                                              |
| structure   | every relaxed structure has an 8-11 residue peptide                  | PASS     | found                                                    |
| structure   | relaxed ids key into test.csv                                        | PASS     | all present                                              |
| md          | trajectory count                                                     | PASS     | 207 (69 complexes x 3 replicas)                          |
| md          | every complex has all 3 replicas                                     | PASS     | 69 complexes                                             |
| md          | every trajectory has 100 frames                                      | PASS     | frames [100]                                             |
| md          | md_residue_contacts.csv: ids all present in test.csv                 | PASS     | 69 distinct                                              |
| md          | distances finite and physical                                        | PASS     | 2.74-21.22 A                                             |
| md          | contact count and contact fraction agree                             | PASS     | count>0 iff some frame within 4.5 A                      |
| md          | no duplicate (tcr_id, replica, chain, resid)                         | PASS     | 5391 rows                                                |
| md          | CDR3 offsets contiguous from 0 in every chain                        | PASS     |                                                          |
| md          | interface trace matches the trajectory count                         | PASS     | 207 traces                                               |
| md          | pair table has 20 discordant and 20 control                          | PASS     | 40 pairs                                                 |
| md          | site residues stored per pair member                                 | PASS     | no shared union                                          |
| md          | pairs whose members number the site differently                      | NOTE     | 5/40 -- each would be mis-probed by a union              |
| md          | every pair member was simulated                                      | PASS     | 69 members, 0 missing                                    |
| submissions | density3_struct03.csv: validates                                     | PASS     | 2636 rows, public 0.8498                                 |
| submissions | density4_struct03.csv: validates                                     | PASS     | 2636 rows, public 0.8525                                 |
| submissions | density_sceptr_default_sigma0.1.csv: validates                       | PASS     | 2636 rows, public 0.8396                                 |
| reports     | every report a script declares exists                                | PASS     | all present                                              |
| reports     | report files on disk                                                 | NOTE     | 19 files, 15 declared by scripts                         |
