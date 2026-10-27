# E70 - validating TCR63's hinge

## Synthetic recovery (must find a real hinge, must not invent one)

| scenario                 |   true_change |   recovered_change_median |   held_out_gain_median |
|:-------------------------|--------------:|--------------------------:|-----------------------:|
| true hinge, change = -60 |           -60 |                  -61.6044 |                 0.1484 |
| true hinge, change = -25 |           -25 |                  -30.495  |                -0.0209 |
| pure line, change = 0    |             0 |                    0.0437 |                -0.0059 |

## Real data

| tcr   | validation     |   linear |   hinge |    gain |   knot_ci_lo |   knot_ci_hi |
|:------|:---------------|---------:|--------:|--------:|-------------:|-------------:|
| TCR59 | gene holdout   |   0.6835 |  0.6634 | -0.0201 |     nan      |     nan      |
| TCR59 | origin holdout |   0.6798 |  0.6415 | -0.0382 |       0.2299 |       0.5774 |
| TCR63 | gene holdout   |   0.3091 |  0.5757 |  0.2666 |     nan      |     nan      |
| TCR63 | origin holdout |   0.3185 |  0.6285 |  0.31   |       0.2699 |       0.4469 |

## Saturation or reversal

| tcr   |   knot |   slope_after |     ci_lo |    ci_hi |   p_negative |   rho_above_knot |   n_above |
|:------|-------:|--------------:|----------:|---------:|-------------:|-----------------:|----------:|
| TCR59 | 0.3272 |       12.8595 |  -45.7509 |  52.9005 |        0.355 |           0.0504 |        94 |
| TCR63 | 0.3206 |      -65.4587 | -134.022  | -30.2143 |        1     |          -0.36   |       119 |

## Per-fold fit status

| tcr   | validation     |   fold | converged   |   knot |   slope_before |   slope_after |   slope_change |   n_left |   n_right | moved_off_start   | knot_at_bound   |
|:------|:---------------|-------:|:------------|-------:|---------------:|--------------:|---------------:|---------:|----------:|:------------------|:----------------|
| TCR59 | gene holdout   |      0 | True        | 0.3216 |        130.53  |       17.7587 |      -112.772  |       80 |        73 | True              | False           |
| TCR59 | gene holdout   |      1 | True        | 0.3272 |        139.168 |        4.2064 |      -134.961  |       82 |        71 | True              | False           |
| TCR59 | gene holdout   |      2 | True        | 0.5156 |        100.932 |      -18.2227 |      -119.155  |       96 |        58 | True              | False           |
| TCR59 | gene holdout   |      3 | True        | 0.3216 |        113.668 |       35.8151 |       -77.8529 |       78 |        76 | True              | False           |
| TCR59 | gene holdout   |      4 | True        | 0.3216 |        133.817 |        9.7302 |      -124.087  |       74 |        80 | True              | False           |
| TCR59 | origin holdout |      0 | True        | 0.4618 |        120.09  |      -59.5344 |      -179.624  |       41 |        28 | True              | False           |
| TCR59 | origin holdout |      1 | True        | 0.3272 |        129.012 |       18.0838 |      -110.929  |       60 |        63 | True              | False           |
| TCR63 | gene holdout   |      0 | True        | 0.3116 |        249.768 |      -60.4237 |      -310.192  |       57 |        96 | True              | False           |
| TCR63 | gene holdout   |      1 | True        | 0.3223 |        242.191 |      -52.2004 |      -294.391  |       59 |        94 | True              | False           |
| TCR63 | gene holdout   |      2 | True        | 0.388  |        213.119 |     -105.694  |      -318.812  |       79 |        75 | True              | False           |
| TCR63 | gene holdout   |      3 | True        | 0.3066 |        262.558 |      -60.6172 |      -323.175  |       53 |       101 | True              | False           |
| TCR63 | gene holdout   |      4 | True        | 0.3737 |        217.14  |      -80.9562 |      -298.097  |       73 |        81 | True              | False           |
| TCR63 | origin holdout |      0 | True        | 0.2789 |        294.461 |      -86.6023 |      -381.063  |       17 |        25 | True              | False           |
| TCR63 | origin holdout |      1 | True        | 0.367  |        216.21  |      -67.24   |      -283.45   |       62 |        88 | True              | False           |
