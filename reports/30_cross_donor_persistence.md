# E57 - cross-donor persistence of the six public families

Donor = reference.id + meta.subject.id. Identical alpha-beta clonotypes (6 in the cohort) never count as each other's neighbours.

## Cross-donor support by group

| group                         |   n |   has_cross_donor_cons |   median_donors_cons |   has_cross_donor_perm |   median_donors_perm |   median_neighbours |
|:------------------------------|----:|-----------------------:|---------------------:|-----------------------:|---------------------:|--------------------:|
| GIL large TRBV19 families     | 205 |                 0.9707 |                 19   |                 0.9951 |                 19   |                13   |
| GIL TRBV19, no large family   | 318 |                 0.2421 |                  0   |                 0.2673 |                  0   |                 0   |
| GIL other V, no large family  | 147 |                 0.0136 |                  0   |                 0.0272 |                  0   |                 0   |
| GLCTLVAML#2 (12 donors, thin) |  32 |                 0.9062 |                  3.5 |                 1      |                  3.5 |                 7.5 |

## Leave-one-donor-out

| family       |   labelled |   donors |   donors_tested |   survival |   worst_donor |
|:-------------|-----------:|---------:|----------------:|-----------:|--------------:|
| GILGFVFTL#5  |         50 |       68 |              68 |     1      |        1      |
| GILGFVFTL#1  |         31 |       51 |              51 |     0.9982 |        0.9091 |
| GLCTLVAML#2  |         27 |       12 |              12 |     1      |        1      |
| GILGFVFTL#27 |         27 |       53 |              53 |     1      |        1      |
| GILGFVFTL#36 |         27 |       69 |              69 |     1      |        1      |
| GILGFVFTL#30 |         23 |      112 |             112 |     1      |        1      |

## Cross-donor support against validation (CONFOUNDED with having any neighbour)

| epitope   | group                   |   n |   validated |
|:----------|:------------------------|----:|------------:|
| GILGFVFTL | has cross-donor support | 223 |      0.9462 |
| GILGFVFTL | no cross-donor support  | 352 |      0.483  |
| GLCTLVAML | has cross-donor support |  50 |      0.94   |
| GLCTLVAML | no cross-donor support  | 116 |      0.5431 |

## The same, conditioned on already having a neighbour

| epitope   | group                      |   n |   validated |   fisher_p |   odds_ratio |
|:----------|:---------------------------|----:|------------:|-----------:|-------------:|
| GILGFVFTL | cross-donor neighbour      | 223 |      0.9462 |     0      |      16.6065 |
| GILGFVFTL | only same-donor neighbours |  35 |      0.5143 |     0      |      16.6065 |
| GLCTLVAML | cross-donor neighbour      |  50 |      0.94   |     0.0001 |      10.4444 |
| GLCTLVAML | only same-donor neighbours |  40 |      0.6    |     0.0001 |      10.4444 |
