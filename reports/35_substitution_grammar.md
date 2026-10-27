# E62 - substitution grammar on the 4TCR NNK block

Response is `r3_r2_ratio` (later-round gain, per E61). Holdout is by connected component of the Hamming-1 graph.

## Held-out Spearman by model

| model                |     1G4 |    DMF5 |     EBV |     Flu |
|:---------------------|--------:|--------:|--------:|--------:|
| BLOSUM62             | -0.0789 | -0.0061 | -0.0242 |  0.0135 |
| Grantham             | -0.0028 |  0.0063 | -0.0073 | -0.006  |
| physicochemical      |  0.0747 | -0.1799 | -0.0151 |  0.0405 |
| position + Grantham  |  0.0944 |  0.2042 | -0.0008 |  0.0209 |
| position only        |  0.0889 |  0.1908 | -0.0295 |  0.009  |
| position x chemistry |  0.1264 |  0.177  | -0.0161 |  0.0623 |

## Negative control (chemistry shuffled within position)

| screen   |    real |   chemistry_shuffled |
|:---------|--------:|---------------------:|
| 1G4      |  0.0944 |               0.0799 |
| DMF5     |  0.2042 |               0.2088 |
| EBV      | -0.0008 |              -0.0315 |
| Flu      |  0.0209 |               0.0121 |

## Cross-receptor agreement of substitution effects

|      |    1G4 |   DMF5 |    EBV |    Flu |
|:-----|-------:|-------:|-------:|-------:|
| 1G4  | 1      | 0.1861 | 0.002  | 0.313  |
| DMF5 | 0.1861 | 1      | 0.212  | 0.1917 |
| EBV  | 0.002  | 0.212  | 1      | 0.1102 |
| Flu  | 0.313  | 0.1917 | 0.1102 | 1      |
