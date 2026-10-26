# Rebuilding the rival configuration in our harness

## Public scores

| configuration             |   public |
|:--------------------------|---------:|
| m30 pool-expanded top-8   |   0.8526 |
| m50 V-domain top-8        |   0.8319 |
| m30 + m50 (their S04)     |   0.8549 |
| m30 only, no pool (floor) |   0.8414 |
| our frozen entry          |   0.8498 |
| tcrdist3 baseline         |   0.8247 |

## After private-realisable downsampling

| private epitope   |   multiplier |   their config |
|:------------------|-------------:|---------------:|
| LTDEMIAQY         |        1.038 |         0.8469 |
| CINGVCWTV         |        1.045 |         0.8469 |
| LLWNGPMAV         |        1.834 |         0.8491 |

## Under E47's published anchor

| configuration                |   public |   margin used |   private lo |   private hi |   private mid |
|:-----------------------------|---------:|--------------:|-------------:|-------------:|--------------:|
| our frozen entry             |   0.8498 |        0.0251 |       0.7784 |       0.8051 |        0.7918 |
| their S04, public margin     |   0.8549 |        0.0302 |       0.7835 |       0.8102 |        0.7968 |
| their S04, realisable margin |   0.8549 |        0.0229 |       0.7762 |       0.8029 |        0.7896 |
