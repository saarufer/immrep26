# Minimisation robustness of the interface descriptors

69 complexes relaxed under one fixed protocol (PDBFixer terminal atoms and hydrogens at pH 7, Amber14 with GBn2 implicit solvent, L-BFGS), measured at [0, 200, 1000] cumulative iterations.

Stability differs sharply between descriptors, which decides which of them a conclusion may rest on.

| descriptor   |   median |change| 0->1000 |   spearman(raw, relaxed) |   moves >20% |   n |
|:-------------|--------------------------:|-------------------------:|-------------:|----:|
| contacts     |                    0.0921 |                   0.7927 |       0.1739 |  69 |
| hbonds       |                    0.3333 |                   0.7877 |       0.6377 |  69 |
| salt_bridges |                    0      |                   1      |       0      |  69 |
| buried_area  |                    0.0163 |                   0.9605 |       0      |  69 |

The position-specific test of E31 is not here: it needs contact at the *mutated CDR3 position*, which this script enabled by saving the relaxed coordinates to `data/interim/relaxed/`. See the E37 entry and its revision in `reports/00_experiment_log.md`, and `scripts/35_md_analysis.py --static`.
