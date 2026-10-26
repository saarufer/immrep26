# E55 - how much of "singleton" is sample size?

## Observed, at each cohort's own size

| epitope    |   n | split   |   singleton_all |   singleton_pos |   singleton_neg |
|:-----------|----:|:--------|----------------:|----------------:|----------------:|
| CINGVCWTV  | 243 | private |          0.9753 |        nan      |        nan      |
| GILGFVFTL  | 670 | public  |          0.5239 |          0.399  |          0.8505 |
| GLCTLVAML  | 192 | public  |          0.4479 |          0.3545 |          0.6607 |
| LLWNGPMAV  | 229 | private |          0.9563 |        nan      |        nan      |
| LTDEMIAQY  | 130 | private |          0.8538 |        nan      |        nan      |
| NLVPMVATV  | 366 | public  |          0.8169 |          0.4286 |          0.9425 |
| TTDPSFLGRY | 387 | public  |          0.9535 |          0.9176 |          1      |
| YLQPRTFLL  | 419 | public  |          0.8162 |          0.5886 |          0.9598 |

## Public cohorts thinned to the private sizes

| private epitope   |   n |   actual singleton rate |   public thinned, mean |   public thinned, max |   public epitopes available |
|:------------------|----:|------------------------:|-----------------------:|----------------------:|----------------------------:|
| LTDEMIAQY         | 130 |                  0.8538 |                 0.7953 |                0.9758 |                           5 |
| LLWNGPMAV         | 229 |                  0.9563 |                 0.8274 |                0.9643 |                           4 |
| CINGVCWTV         | 243 |                  0.9753 |                 0.8247 |                0.9631 |                           4 |

## Full sweep

| epitope    |   n |   singleton_all |   singleton_pos |   singleton_neg |   lost_neighbour_pos |   lost_neighbour_neg |
|:-----------|----:|----------------:|----------------:|----------------:|---------------------:|---------------------:|
| GILGFVFTL  |  80 |          0.7422 |          0.6708 |          0.9553 |               0.4522 |               0.7195 |
| GILGFVFTL  | 100 |          0.7153 |          0.6403 |          0.9474 |               0.4075 |               0.681  |
| GILGFVFTL  | 130 |          0.697  |          0.6159 |          0.9395 |               0.357  |               0.5973 |
| GILGFVFTL  | 160 |          0.6649 |          0.5708 |          0.9323 |               0.292  |               0.5455 |
| GILGFVFTL  | 192 |          0.6463 |          0.5483 |          0.9213 |               0.2577 |               0.5002 |
| GILGFVFTL  | 229 |          0.6321 |          0.5311 |          0.9103 |               0.217  |               0.409  |
| GILGFVFTL  | 243 |          0.6249 |          0.5201 |          0.9105 |               0.2029 |               0.4054 |
| GILGFVFTL  | 300 |          0.6036 |          0.4962 |          0.8951 |               0.164  |               0.3159 |
| GILGFVFTL  | 366 |          0.5812 |          0.4661 |          0.8824 |               0.1149 |               0.2248 |
| GILGFVFTL  | 387 |          0.5743 |          0.4569 |          0.8807 |               0.1016 |               0.2138 |
| GILGFVFTL  | 419 |          0.5707 |          0.4549 |          0.8784 |               0.0916 |               0.1862 |
| GILGFVFTL  | 500 |          0.5516 |          0.4307 |          0.8665 |               0.0527 |               0.1106 |
| GILGFVFTL  | 670 |          0.5239 |          0.399  |          0.8505 |               0      |               0      |
| GLCTLVAML  |  80 |          0.5901 |          0.5028 |          0.7607 |               0.2347 |               0.3095 |
| GLCTLVAML  | 100 |          0.5588 |          0.4723 |          0.7371 |               0.1878 |               0.2264 |
| GLCTLVAML  | 130 |          0.5054 |          0.4134 |          0.7005 |               0.0931 |               0.1186 |
| GLCTLVAML  | 160 |          0.4757 |          0.3833 |          0.6826 |               0.0445 |               0.0606 |
| GLCTLVAML  | 192 |          0.4479 |          0.3545 |          0.6607 |               0      |               0      |
| NLVPMVATV  |  80 |          0.9202 |          0.7319 |          0.9814 |               0.544  |               0.7204 |
| NLVPMVATV  | 100 |          0.9113 |          0.6958 |          0.9836 |               0.4793 |               0.7257 |
| NLVPMVATV  | 130 |          0.8943 |          0.6445 |          0.977  |               0.3817 |               0.6195 |
| NLVPMVATV  | 160 |          0.8797 |          0.6001 |          0.9716 |               0.3045 |               0.5383 |
| NLVPMVATV  | 192 |          0.8685 |          0.5675 |          0.9681 |               0.2466 |               0.4607 |
| NLVPMVATV  | 229 |          0.8544 |          0.5298 |          0.9606 |               0.1843 |               0.3348 |
| NLVPMVATV  | 243 |          0.8526 |          0.5232 |          0.9598 |               0.1638 |               0.3161 |
| NLVPMVATV  | 300 |          0.835  |          0.4714 |          0.9532 |               0.078  |               0.1787 |
| NLVPMVATV  | 366 |          0.8169 |          0.4286 |          0.9425 |               0      |               0      |
| TTDPSFLGRY |  80 |          0.9836 |          0.9729 |          1      |               0.7348 |             nan      |
| TTDPSFLGRY | 100 |          0.9791 |          0.9642 |          1      |               0.6444 |             nan      |
| TTDPSFLGRY | 130 |          0.9758 |          0.9581 |          1      |               0.5286 |             nan      |
| TTDPSFLGRY | 160 |          0.9741 |          0.956  |          1      |               0.4916 |             nan      |
| TTDPSFLGRY | 192 |          0.969  |          0.9462 |          1      |               0.3762 |             nan      |
| TTDPSFLGRY | 229 |          0.9643 |          0.9381 |          1      |               0.2693 |             nan      |
| TTDPSFLGRY | 243 |          0.9631 |          0.937  |          1      |               0.2572 |             nan      |
| TTDPSFLGRY | 300 |          0.9587 |          0.928  |          1      |               0.1356 |             nan      |
| TTDPSFLGRY | 366 |          0.9549 |          0.9207 |          1      |               0.0321 |             nan      |
| TTDPSFLGRY | 387 |          0.9535 |          0.9176 |          1      |               0      |             nan      |
| YLQPRTFLL  |  80 |          0.9298 |          0.836  |          0.9923 |               0.6095 |               0.8408 |
| YLQPRTFLL  | 100 |          0.9224 |          0.8191 |          0.9883 |               0.5697 |               0.7675 |
| YLQPRTFLL  | 130 |          0.9041 |          0.774  |          0.9859 |               0.4562 |               0.7015 |
| YLQPRTFLL  | 160 |          0.8884 |          0.7418 |          0.9821 |               0.3739 |               0.5902 |
| YLQPRTFLL  | 192 |          0.8752 |          0.7065 |          0.9814 |               0.2921 |               0.557  |
| YLQPRTFLL  | 229 |          0.8587 |          0.6748 |          0.9753 |               0.2159 |               0.4389 |
| YLQPRTFLL  | 243 |          0.8582 |          0.6737 |          0.9761 |               0.1995 |               0.4221 |
| YLQPRTFLL  | 300 |          0.8389 |          0.6296 |          0.9707 |               0.1091 |               0.278  |
| YLQPRTFLL  | 366 |          0.8257 |          0.6073 |          0.9633 |               0.0455 |               0.1027 |
| YLQPRTFLL  | 387 |          0.8207 |          0.5974 |          0.9616 |               0.0236 |               0.0616 |
| YLQPRTFLL  | 419 |          0.8162 |          0.5886 |          0.9598 |               0      |               0      |

## Positives against negatives, full size against n=130

| epitope    |   n_full |   singleton_pos_full |   singleton_pos_130 |   pos_increase |   singleton_neg_full |   singleton_neg_130 |   neg_increase |
|:-----------|---------:|---------------------:|--------------------:|---------------:|---------------------:|--------------------:|---------------:|
| GILGFVFTL  |      670 |               0.399  |              0.6159 |         0.217  |               0.8505 |              0.9395 |         0.0889 |
| GLCTLVAML  |      192 |               0.3545 |              0.4134 |         0.0588 |               0.6607 |              0.7005 |         0.0398 |
| NLVPMVATV  |      366 |               0.4286 |              0.6445 |         0.2159 |               0.9425 |              0.977  |         0.0345 |
| TTDPSFLGRY |      387 |               0.9176 |              0.9581 |         0.0404 |               1      |              1      |         0      |
| YLQPRTFLL  |      419 |               0.5886 |              0.774  |         0.1854 |               0.9598 |              0.9859 |         0.0261 |

Caveat: this subsamples a finite candidate library and cannot reconstruct a real antigen-specific repertoire, so it bounds the sampling explanation from below.
