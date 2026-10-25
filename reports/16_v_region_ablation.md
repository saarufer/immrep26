# Where in the V domain does the alignment signal sit?

Framework contacts nothing in the interface; CDR1, CDR2 and CDR2.5 contact peptide and MHC. Splitting E26's CDR3-masked region between them tests whether its signal is germline confounding or compatible with recognition.

## Mean per-epitope AUC

|           |   sigma=0.02 |   sigma=0.05 |   sigma=0.1 |   sigma=0.2 |   best |
|:----------|-------------:|-------------:|------------:|------------:|-------:|
| full      |       0.8243 |       0.8304 |      0.8354 |      0.8343 | 0.8354 |
| no_cdr3   |       0.8199 |       0.826  |      0.8276 |      0.8264 | 0.8276 |
| loops     |       0.8074 |       0.8092 |      0.8122 |      0.8173 | 0.8173 |
| framework |       0.8154 |       0.8246 |      0.827  |      0.8223 | 0.827  |

## Within-epitope Spearman between regions (sigma=0.1)

|           |   full |   no_cdr3 |   loops |   framework |
|:----------|-------:|----------:|--------:|------------:|
| full      |  1     |     0.985 |   0.949 |       0.984 |
| no_cdr3   |  0.985 |     1     |   0.97  |       0.994 |
| loops     |  0.949 |     0.97  |   1     |       0.961 |
| framework |  0.984 |     0.994 |   0.961 |       1     |

## Verdict

Framework carries it. Framework contacts nothing in the interface, so this is germline confounding and E26's original reading stands.
