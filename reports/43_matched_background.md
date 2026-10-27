# E78 - matched-background neighbourhood enrichment

Background matched per query on TRBV stem x CDR3-beta length bin x Pgen tercile, relaxed when fewer than 30 candidates.

## Match levels

|                               |    n |
|:------------------------------|-----:|
| len+pgen                      | 2039 |
| trbv+len                      |  380 |
| trbv+len+pgen                 |  150 |
| len only                      |   52 |
| any foreign (E8's background) |   15 |

## Feature performance

| feature                       |   macro |   GIL |    GLC |    NLV |    TTD |    YLQ |
|:------------------------------|--------:|------:|-------:|-------:|-------:|-------:|
| own-cohort top-8 (raw)        |  0.8414 | 0.88  | 0.714  | 0.833  | 0.8505 | 0.9296 |
| matched background alone      |  0.5881 | 0.596 | 0.6942 | 0.5736 | 0.6011 | 0.4755 |
| enrichment = own - background |  0.8154 | 0.883 | 0.6753 | 0.8122 | 0.7928 | 0.9139 |

- D as submitted 0.8600
- enrichment's within-epitope rho with D: +0.826
- oracle global weight: 0.8600 (w=0.00)
- nested LOEO: 0.8587
