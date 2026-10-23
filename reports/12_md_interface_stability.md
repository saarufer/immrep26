# C11: is the interface difference dynamically stable?

69 complexes, 207 trajectories.

## Group comparison at the mutated site

| descriptor                           |   discordant |   control |   n_disc |   n_ctl |      p |
|:-------------------------------------|-------------:|----------:|---------:|--------:|-------:|
| contacts, time-averaged              |       2.4562 |    1.332  |       62 |      64 | 0.1588 |
| distance, time-averaged (A)          |       6.817  |    6.8034 |       62 |      64 | 0.5451 |
| fraction of frames in contact        |       0.2878 |    0.1536 |       62 |      64 | 0.1486 |
| touches at any point in MD (binary)  |       0.5968 |    0.5938 |       62 |      64 | 1      |
| frames in contact | faced statically |       0.7582 |    0.7081 |       20 |       9 | 0.9435 |
| frames in contact | did not face     |       0.0638 |    0.0628 |       42 |      55 | 0.5787 |

## Effect against replica-to-replica spread

| quantity         |   effect |   replica_sd_of_effect | per_replica           | sign_consistent   | verdict                                                 |
|:-----------------|---------:|-----------------------:|:----------------------|:------------------|:--------------------------------------------------------|
| contacts         |   1.1242 |                 0.3068 | [1.267, 0.772, 1.334] | True              | exceeds spread                                          |
| distance         |   0.0136 |                 0.1701 | [-0.025, 0.2, -0.134] | False             | WITHIN SPREAD -- no result; sign flips between replicas |
| contact fraction |   0.1342 |                 0.025  | [0.117, 0.122, 0.163] | True              | exceeds spread                                          |

## By epitope

| epitope   |   ('count_mean', 'control') |   ('count_mean', 'discordant') |   ('dist_mean', 'control') |   ('dist_mean', 'discordant') |
|:----------|----------------------------:|-------------------------------:|---------------------------:|------------------------------:|
| GILGFVFTL |                       0.641 |                          1.161 |                      7.134 |                         7.589 |
| GLCTLVAML |                       9.357 |                          2.835 |                      3.537 |                         5.716 |
| NLVPMVATV |                       0.187 |                          5.322 |                      6.769 |                         4.881 |
| YLQPRTFLL |                       4.442 |                          0     |                      5.753 |                        10.616 |
