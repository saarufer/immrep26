# E79 - reference annotation quality weighting

| epitope    | split   |   n_refs |   n_meth |   score |
|:-----------|:--------|---------:|---------:|--------:|
| CINGVCWTV  | private |        2 |        1 |       2 |
| GILGFVFTL  | public  |       11 |        7 |       3 |
| GLCTLVAML  | public  |        7 |        4 |       4 |
| LLWNGPMAV  | private |        2 |        2 |       1 |
| LTDEMIAQY  | private |        1 |        1 |       1 |
| NLVPMVATV  | public  |        5 |        4 |       4 |
| TTDPSFLGRY | public  |        1 |        1 |       1 |
| YLQPRTFLL  | public  |        3 |        2 |       3 |

| weighting   |   macro |   delta |
|:------------|--------:|--------:|
| unweighted  |  0.8414 |  0      |
| n_refs      |  0.8421 |  0.0007 |
| score       |  0.8388 | -0.0026 |
| n_meth      |  0.8409 | -0.0005 |

Into D: oracle 0.8600 at w=0.00.
