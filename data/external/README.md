# External data (NOT committed)

Reproduce by running `scripts/00_fetch_external.sh`, except TCRvdb which must be
requested. Sizes and deletion guidance are in `../../ARCHIVE.md`.

| Dataset | Source | Licence / access | Local path |
|---|---|---|---|
| VDJdb `2026-06-03` | https://github.com/antigenomics/vdjdb-db/releases | CC-BY-NC-SA (see bundled LICENSE) | `vdjdb/vdjdb-2026-06-03/` |
| VDJdb `2025-02-21`, `2024-06-13`, `2023-06-01` | same | same | `vdjdb/` — pre-revalidation snapshot (E92) only |
| TCRvdb `2026-10-01` | https://tcrvdb.nkiai.nl (**request required**) | **Redistribution prohibited.** Cite doi:10.1016/j.immuni.2026.09.007 | `tcrvdb/` |
| IMMREP26 competition files | Kaggle `immrep26` | Kaggle competition rules | `../raw/` |
| PDB sequence index + 33 mmCIF entries | https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz; https://files.rcsb.org/download/<ID>.cif | public domain (CC0) | `pdb/` — crystallography comparison (E34) |
| PRP yeast-display peptide counts | Wang N. et al. 2026, Nat Biotechnol doi:10.1038/s41587-026-03128-x; data doi:10.5061/dryad.h70rxwdzk | Dryad CC0 | `prp/` |
| SEEKER library-on-library screens | Li S. et al. 2026, bioRxiv doi:10.64898/2026.09.25.753432; data doi:10.5061/dryad.79cnp5jct | Dryad CC0 | `seeker/` |
| IMGT/HLA class-I α chains; UniProt P01892, P30443 | https://www.ebi.ac.uk/ipd/imgt/hla/; https://www.uniprot.org/ | IMGT terms; UniProt CC-BY | `hla/` |

Neither PRP nor SEEKER contributed to any submitted model; they were evaluated
as independent information sources in E63 and E61 respectively
(`reports/00_experiment_log.md`).
