# Predicted against crystallographic contacts

TCR, peptide and MHC sequences taken from solved complexes and folded
through the same Boltz pipeline used for the competition, then compared to
the interface measured in the crystal.

| complex        | epitope   |   rho_counts |   rho_mindist |   jaccard_touched_positions |   pred_contacts |   real_contacts |
|:---------------|:----------|-------------:|--------------:|----------------------------:|----------------:|----------------:|
| 1oga_GILGFVFTL | GILGFVFTL |        0.992 |         0.991 |                       1     |              75 |              72 |
| 2vlj_GILGFVFTL | GILGFVFTL |        0.999 |         0.994 |                       1     |              78 |              74 |
| 2vlk_GILGFVFTL | GILGFVFTL |        0.996 |         0.983 |                       1     |              77 |              68 |
| 2vlr_GILGFVFTL | GILGFVFTL |        0.992 |         0.997 |                       1     |              71 |              66 |
| 3gsn_NLVPMVATV | NLVPMVATV |        0.964 |         0.983 |                       0.9   |             130 |             121 |
| 3o4l_GLCTLVAML | GLCTLVAML |        0.955 |         0.989 |                       0.875 |              87 |              84 |
| 5e6i_GILGFVFTL | GILGFVFTL |        0.811 |         0.982 |                       0.667 |              88 |              77 |
| 5euo_GILGFVFTL | GILGFVFTL |        0.997 |         0.983 |                       1     |              99 |              88 |
| 5isz_GILGFVFTL | GILGFVFTL |        0.986 |         0.983 |                       1     |             102 |             100 |
| 5jhd_GILGFVFTL | GILGFVFTL |        0.997 |         0.98  |                       1     |              81 |              82 |
| 5tez_GILGFVFTL | GILGFVFTL |        0.948 |         0.971 |                       0.875 |              79 |              73 |
| 7n1f_YLQPRTFLL | YLQPRTFLL |        0.992 |         0.973 |                       1     |             131 |             151 |
| 7n6e_YLQPRTFLL | YLQPRTFLL |        0.996 |         0.988 |                       1     |             144 |             155 |
| 7pbe_YLQPRTFLL | YLQPRTFLL |        0.976 |         0.98  |                       0.9   |             147 |             165 |
| 9x7r_AMDLGIHKV | AMDLGIHKV |        0.84  |         0.848 |                       0.7   |             146 |             104 |
