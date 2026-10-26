# Gate 0 - external data audit (PRP + SEEKER)

Both datasets are CC0-1.0. Manifest: `data/derived/external_manifest.tsv`, 76 files.

## PRP selection files

| tcr       | anchor   |   peptides |   in_own_input_library |   input_coverage |   reads_R1 |   reads_R4 | engineered   |
|:----------|:---------|-----------:|-----------------------:|-----------------:|-----------:|-----------:|:-------------|
| TCR019.1  | P2_P8    |     250802 |                  28751 |           0.1146 |     257620 |     708814 | False        |
| TCR135.1  | P2_P8    |      82846 |                  11830 |           0.1428 |      67708 |     175015 | False        |
| TCR135.3  | P2_P8    |     174812 |                  21638 |           0.1238 |     165754 |     325353 | False        |
| TCR135.4  | P2_P8    |     553682 |                  46326 |           0.0837 |     613008 |     949325 | False        |
| TCR135.7  | P2_P8    |     474951 |                  42191 |           0.0888 |     626326 |    1254066 | False        |
| TCR135.8  | P2_P8    |     341438 |                  35557 |           0.1041 |     355345 |     312551 | False        |
| TCR19.2C1 | P2_P8    |     478596 |                  39283 |           0.0821 |     486247 |     592045 | True         |
| TCR19.2C2 | P2_P8    |     358832 |                  32473 |           0.0905 |     281830 |     515814 | True         |
| TCR19.2C3 | P2_P8    |     315397 |                  30994 |           0.0983 |     209681 |     825954 | True         |
| TCR19.2C4 | P2_P8    |     330492 |                  31502 |           0.0953 |     281503 |     536763 | True         |
| TCR19.2C5 | P2_P8    |     370349 |                  36350 |           0.0982 |     437819 |     298079 | True         |
| TCR19.2   | P2_P8    |     776210 |                  53476 |           0.0689 |    1011651 |     938256 | False        |
| TCR26.1   | P2_P8    |     362011 |                  35932 |           0.0993 |     550853 |     680257 | False        |
| TCR26.2   | P2_P8    |     522892 |                  46050 |           0.0881 |     881373 |     522689 | False        |
| TCR28.1   | P2_P8    |     328650 |                  34013 |           0.1035 |     420215 |     353731 | False        |
| TCR4.1    | P2_P8    |     540627 |                  46476 |           0.086  |     676704 |     537136 | False        |
| TCR4.2    | P2_P9    |     685323 |                   1104 |           0.0016 |     540273 |     389149 | False        |
| TCR4.3    | P2_P9    |     667191 |                   1059 |           0.0016 |     461143 |     364765 | False        |
| TCR4.4    | P2_P9    |     627776 |                    988 |           0.0016 |     381159 |     411706 | False        |
| TCR8.4    | P2_P8    |     438629 |                  38964 |           0.0888 |     602211 |     492423 | False        |
| TCR9.1    | P2_P8    |     456867 |                  39754 |           0.087  |     569180 |     422526 | False        |

Anchor libraries: P2_P8 272,888 peptides, P2_P9 151,914, sharing 256. Kept apart.


**Data quirk:** TCR4.2, TCR4.3, TCR4.4 ship without a header row; their first line is already data, and they are exactly the three P2/P9 files. A naive directory read loses one peptide from each and mislabels their columns.

## SEEKER B27 screens

| tcr   |   candidates |   has_round0 |   round0_nonzero |   add_flow_single_copy |   add_flow_single_copy_pct | add_flow_single_copy_absent   |   add_flow_high_copy |   add_flow_high_copy_pct | add_flow_high_copy_absent   |
|:------|-------------:|-------------:|-----------------:|-----------------------:|---------------------------:|:------------------------------|---------------------:|-------------------------:|:----------------------------|
| TCR59 |        19322 |        19322 |            19322 |                    192 |                     0.0099 | False                         |                  192 |                   0.0099 | False                       |
| TCR63 |        16369 |        16369 |            16369 |                    192 |                     0.0117 | False                         |                    0 |                   0      | True                        |

6 top-10,000 tables are truncated candidate lists; absence is not a measured negative.

## Overlap

- shared peptide sequences: **19**
- shared TCRs: **0**

## Gate 0

- [FAIL] PRP input library available for the enrichment denominator
- [PASS] SEEKER R0 exists
- [PASS] SEEKER flow measurements exist
- [FAIL] same TCR x same peptide x same HLA across the two sources

**Verdict: profiles only** - build recognition profiles and the within-SEEKER screen-to-flow calibration; do not fit a cross-source binding to activation map.
