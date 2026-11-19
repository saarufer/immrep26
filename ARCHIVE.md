# Archiving this project — what can be deleted, and how to get it back

Written 2026-11-19, when the project was closed. The competition ended
2026-11-02; the final results are in `reports/107_final_results.md`.

The working tree is **38 GB**, plus **7.6 GB** of Boltz model weights outside
it. Almost all of that is regenerable or re-downloadable. Everything needed to
read the project's conclusions — `reports/`, `scripts/`, `src/`, and the one
committed feature cache — is under 10 MB and is already on GitHub.

Short version: **deleting `data/` and `.venv/` recovers 37.6 GB and loses
nothing that cannot be rebuilt.** The only expensive item is `data/interim/boltz/`,
and the features extracted from it are committed, so even that is safe.

## A. Downloaded — delete freely, every source is recorded below

| Path | Size | Source | Re-fetch |
|---|---:|---|---|
| `data/external/vdjdb/vdjdb-2026-06-03/` | 406 MB | VDJdb release 2026-06-03, [antigenomics/vdjdb-db](https://github.com/antigenomics/vdjdb-db/releases) | `scripts/00_fetch_external.sh` |
| `data/external/vdjdb/vdjdb-2025-02-21/`, `-2024-06-13/`, `-2023-06-01/`, `tmp_2025-02-21/` | 501 MB | older VDJdb releases, same source; used only for the pre-revalidation snapshot comparison in E92 | same URL pattern, substitute the tag |
| `data/external/vdjdb-2026-06-03.zip`, `vdjdb-2023-06-01.zip` | 51 MB | the downloaded archives themselves | re-download |
| `data/external/pdb/pdb_seqres.txt(.gz)` | 459 MB | RCSB derived data, `https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz` | `curl` the same URL |
| `data/external/pdb/cif/` (33 entries) | 11 MB | RCSB mmCIF for the solved TCR–pMHC structures used in E34's crystallography comparison | `https://files.rcsb.org/download/<ID>.cif`; IDs in `epitope_entries.json` |
| `data/external/prp/` | 351 MB | Wang N. et al. 2026, *Deep peptide recognition profiling…*, Nat Biotechnol [doi:10.1038/s41587-026-03128-x](https://doi.org/10.1038/s41587-026-03128-x); data [doi:10.5061/dryad.h70rxwdzk](https://doi.org/10.5061/dryad.h70rxwdzk); code [akds/PRP-TCR-Specificity](https://github.com/akds/PRP-TCR-Specificity) | Dryad |
| `data/external/seeker/` | 139 MB | Li S. et al. 2026, *SEEKER*, bioRxiv [doi:10.64898/2026.09.25.753432](https://doi.org/10.64898/2026.09.25.753432); data [doi:10.5061/dryad.79cnp5jct](https://doi.org/10.5061/dryad.79cnp5jct); code [ldi2009/SEEKER-2026](https://github.com/ldi2009/SEEKER-2026) @ `3f4d822`, Zenodo [10.5281/zenodo.23216513](https://doi.org/10.5281/zenodo.23216513) | Dryad |
| `data/external/hla/` | 3.2 MB | IMGT/HLA `A_prot.fasta`; UniProt `P01892` (HLA-A\*02:01), `P30443` (HLA-A\*01:01) | IMGT/HLA FTP, UniProt |
| `data/raw/` | 1.6 MB | Kaggle competition `immrep26` | `kaggle competitions download -c immrep26` |

**`data/external/tcrvdb/` and `data/external/TCRvdb/` (1.1 MB) are different.**
TCRvdb is request-only from [tcrvdb.nkiai.nl](https://tcrvdb.nkiai.nl) and its
terms **prohibit redistribution**. It is gitignored and was never committed.
Delete it with the rest, or keep it locally — but do not copy it anywhere, and
cite [doi:10.1016/j.immuni.2026.09.007](https://doi.org/10.1016/j.immuni.2026.09.007)
in any work that used it. Re-obtaining it means making the request again.

## B. Generated — deleting costs compute, not provenance

| Path | Size | Produced by | Cost to rebuild |
|---|---:|---|---|
| `data/interim/boltz/` | **23 GB** | `scripts/08`–`10`, Boltz-2 complexes for all 2,636 pairs (`inputs/`, `out/`, `out5/`) | ~10 h on an RTX 5090, 14 s per complex |
| `data/interim/decoy/` | 3.9 GB | anchor-matched decoy peptides, the specificity control (E36) | a few GPU-hours |
| `data/interim/realval/`, `relaxed/` | 82 MB | crystallography comparison (E34) and OpenMM relaxation of the E31 interface result (E37) | minutes to hours |
| `data/interim/sceptr/`, `vregion/`, `fullseq/`, `tcrdist/`, `grantham/` | 55 MB | embeddings and distance matrices, `scripts/05` and friends | minutes on GPU |
| `data/interim/olga_beta_200k.npy` | 4.9 MB | OLGA synthetic repertoire | minutes |
| the ~20 small CSVs in `data/interim/` | ~6 MB | one per feature family; each names its producing script in the log | cheap except the structural ones |

**Do not delete `data/interim/boltz_features.csv` (1.8 MB).** It is the one
derived cache that *is* committed, carrying `tcr_id` plus 34 numeric structural
descriptors — no sequences, no labels. It exists precisely so that
`make submission` reproduces D from a fresh clone without a GPU, which is what
makes deleting the 23 GB safe.

Also keep `submissions/` (1.3 MB, gitignored): these are the files actually
submitted, and their SHA-256 values are quoted in the method documents sent to
the organisers.

## C. Environments and model weights

| Path | Size | Note |
|---|---:|---|
| `.venv/` | 8.6 GB | mostly `nvidia` 4.2 GB + `torch` 1.6 GB + `triton` 641 MB. Rebuild per the README. |
| `.mmenv-openmm/` | 1.1 GB | conda-forge env for the OpenMM CUDA plugin that the PyPI wheel lacks; only needed to redo the relaxation in E37. |
| `~/.boltz/` | 7.6 GB | Boltz-2 weights: `boltz2_conf.ckpt` 2.2 GB, `boltz2_aff.ckpt` 2.0 GB, `mols.tar` + `mols/` 3.6 GB. **Auto-redownloaded on first run**, so safe to delete. |
| `~/.cache/hermes/` | 2.4 MB | HERMES Wigner matrices from E148. |

Two caches are **not** this project's and should be left alone:
`~/.cache/uv` (17 GB, shared package cache across all your projects) and
`~/.cache/huggingface` (1.2 GB — inspected 2026-11-19, it holds MedCPT,
DeBERTa, arcinstitute and Xaira artefacts from other work, nothing from
IMMREP26).

## D. Doing it

```bash
cd ~/project_Immrep26

# 1. keep the one committed cache and the submitted files
cp data/interim/boltz_features.csv /tmp/keep_boltz_features.csv
tar czf ~/immrep26_submissions.tgz submissions/

# 2. the 27 GB of generated structures
rm -rf data/interim/boltz data/interim/decoy data/interim/realval data/interim/relaxed

# 3. the 1.9 GB of re-downloadable external data
rm -rf data/external/vdjdb data/external/pdb data/external/prp data/external/seeker \
       data/external/*.zip

# 4. environments (8.6 GB + 1.1 GB + 7.6 GB)
rm -rf .venv .mmenv-openmm ~/.boltz

# 5. restore the committed cache
mkdir -p data/interim && cp /tmp/keep_boltz_features.csv data/interim/boltz_features.csv
```

That leaves roughly 15 MB on disk and the full repository on GitHub. `git status`
should stay clean throughout: everything removed above is gitignored.

TCRvdb (step 3 deliberately omits it) is a separate decision — see §A.
