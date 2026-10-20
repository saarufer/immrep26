# Metadata viability inside private epitopes

VDJdb `2026-06-03`, 75,855 paired human clonotypes.
Cells are the number of distinct non-null levels within that epitope's cohort.
A field counts as usable only if it has >=2 levels in all three private epitopes and its dominant level holds <95% of rows.

| feature | CINGVCWTV | LLWNGPMAV | LTDEMIAQY | GILGFVFTL | YLQPRTFLL | TTDPSFLGRY | NLVPMVATV | GLCTLVAML | max skew (private) | usable |
|---|---|---|---|---|---|---|---|---|---|---|
| `vdjdb.score` | 2 | 1 | 1 | 3 | 4 | 1 | 4 | 4 | 1.00 | **no** |
| `reference.id` | 3 | 4 | 2 | 20 | 9 | 3 | 25 | 14 | 0.98 | **no** |
| `method.identification` | 3 | 2 | 2 | 8 | 3 | 2 | 12 | 9 | 0.98 | **no** |
| `method.verification` | 0 | 0 | 0 | 2 | 2 | 0 | 5 | 2 | nan | **no** |
| `method.singlecell` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 | **no** |
| `method.sequencing` | 3 | 1 | 1 | 3 | 2 | 2 | 3 | 3 | 1.00 | **no** |
| `method.frequency` | 16 | 10 | 0 | 287 | 42 | 0 | 81 | 91 | 0.84 | **no** |
| `meta.cell.subset` | 4 | 1 | 1 | 20 | 2 | 4 | 1 | 6 | 1.00 | **no** |
| `meta.subject.cohort` | 4 | 1 | 0 | 7 | 3 | 0 | 10 | 4 | 1.00 | **no** |
| `meta.subject.id` | 10 | 0 | 0 | 167 | 46 | 8 | 47 | 27 | 0.83 | **no** |
| `meta.tissue` | 1 | 1 | 1 | 5 | 1 | 1 | 3 | 3 | 1.00 | **no** |
| `meta.donor.MHC` | 5 | 1 | 0 | 43 | 1 | 0 | 2 | 4 | 1.00 | **no** |
| `mhc.a` | 1 | 1 | 1 | 3 | 1 | 1 | 3 | 2 | 1.00 | **no** |

**Usable on all three private epitopes:** none

**Dead on at least one private epitope:** `meta.cell.subset`, `meta.donor.MHC`, `meta.subject.cohort`, `meta.subject.id`, `meta.tissue`, `method.frequency`, `method.identification`, `method.sequencing`, `method.singlecell`, `method.verification`, `mhc.a`, `reference.id`, `vdjdb.score`

## Consequence

Annotation-provenance features cannot rank the private epitopes, so leave-one-epitope-out scores on the public epitopes overstate what they will deliver on the final ranking. Treat them as an ablation / diagnostic, not as a submission component, and keep the ranking model on paired-abTCR sequence and structure, which transfer.
