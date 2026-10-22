# C11: is the interface difference dynamically stable?

28 complexes, 83 trajectories.

## Group comparison at the mutated site

| descriptor                    |   discordant |   control |   n_disc |   n_ctl |      p |
|:------------------------------|-------------:|----------:|---------:|--------:|-------:|
| contacts, time-averaged       |       1.2687 |    0.8548 |       28 |      27 | 0.8681 |
| distance, time-averaged (A)   |       7.3674 |    7.2063 |       28 |      27 | 0.8795 |
| fraction of frames in contact |       0.1731 |    0.114  |       28 |      27 | 0.9721 |

## Effect against replica-to-replica spread

| quantity         |   effect |   replica_sd_of_effect | per_replica            | sign_consistent   | verdict                                                 |
|:-----------------|---------:|-----------------------:|:-----------------------|:------------------|:--------------------------------------------------------|
| contacts         |   0.3992 |                 0.2029 | [0.411, 0.191, 0.596]  | True              | exceeds spread                                          |
| distance         |   0.1506 |                 0.2934 | [0.113, 0.461, -0.122] | False             | WITHIN SPREAD -- no result; sign flips between replicas |
| contact fraction |   0.0566 |                 0.0191 | [0.049, 0.043, 0.078]  | True              | exceeds spread                                          |

## By epitope

| epitope   |   ('count_mean', 'control') |   ('count_mean', 'discordant') |   ('dist_mean', 'control') |   ('dist_mean', 'discordant') |
|:----------|----------------------------:|-------------------------------:|---------------------------:|------------------------------:|
| GILGFVFTL |                       0.855 |                          1.269 |                      7.206 |                         7.367 |
