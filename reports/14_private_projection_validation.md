# Does E30's private projection transfer?

## Observable and unobservable quantities

| epitope    |   q_cohort |    pi |     a |     b |     lam |   s_same |
|:-----------|-----------:|------:|------:|------:|--------:|---------:|
| GILGFVFTL  |      0.476 | 0.663 | 0.601 | 0.149 |   4.021 |    0.039 |
| GLCTLVAML  |      0.552 | 0.663 | 0.645 | 0.339 |   1.902 |    0.05  |
| NLVPMVATV  |      0.183 | 0.287 | 0.571 | 0.058 |   9.934 |    0.013 |
| TTDPSFLGRY |      0.047 | 0.497 | 0.082 | 0     | nan     |    0     |
| YLQPRTFLL  |      0.184 | 0.414 | 0.411 | 0.04  |  10.239 |    0.014 |

`q` is label-free. `pi`, `a`, `b` are not, and the class composition needs `a` and `b` separately.

## Leave-one-epitope-out transfer test

| epitope    |   actual |   full |   oracle_composition |   oracle_errors |   d_full |   d_oracle_composition |   d_oracle_errors |
|:-----------|---------:|-------:|---------------------:|----------------:|---------:|-----------------------:|------------------:|
| GILGFVFTL  |   0.8693 | 0.9537 |               0.8763 |          0.9375 |   0.0844 |                 0.007  |            0.0682 |
| GLCTLVAML  |   0.7414 | 0.9543 |               0.833  |          0.9293 |   0.2129 |                 0.0916 |            0.1879 |
| NLVPMVATV  |   0.8649 | 0.8559 |               0.9015 |          0.7823 |  -0.009  |                 0.0366 |           -0.0826 |
| TTDPSFLGRY |   0.8402 | 0.847  |               0.8523 |        nan      |   0.0068 |                 0.0122 |          nan      |
| YLQPRTFLL  |   0.9334 | 0.8299 |               0.8637 |          0.9115 |  -0.1035 |                -0.0697 |           -0.0219 |

|                    |   mean_error |    mae |   max_abs |
|:-------------------|-------------:|-------:|----------:|
| full               |       0.0383 | 0.0833 |    0.2129 |
| oracle_composition |       0.0156 | 0.0434 |    0.0916 |
| oracle_errors      |       0.0379 | 0.0902 |    0.1879 |

## Private sensitivity surface

|   lambda |   0.3 |   0.4 |   0.5 |   0.6 |   0.7 |
|---------:|------:|------:|------:|------:|------:|
|  1.90239 | 0.829 | 0.83  | 0.83  | 0.83  | 0.83  |
|  6.97744 | 0.852 | 0.849 | 0.847 | 0.845 | 0.844 |
| 10.2391  | 0.857 | 0.852 | 0.85  | 0.848 | 0.846 |
