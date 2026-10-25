# IMMREP26 — reproduction entry points.
#
# The predictor is frozen (2026-10-25). `make submission` rebuilds it from the
# cached features and checks the result is byte-identical to what was entered;
# that check is the project's reproducibility claim, so it is the default
# target. Everything expensive is opt-in and named for what it costs.
#
# See docs/REPRODUCING.md for what each tier needs and why the GPU stages are
# separated.

PY      := .venv/bin/python
MMENV   := ./.mmenv-openmm
FROZEN  := submissions/density3_struct03.csv
SHA     := 523a4692df02467b3d9b4fe0e51136a4aa4df3af7ddae2262cf8ed282ac5608c

.DEFAULT_GOAL := check
.PHONY: check test audit submission verify qc analysis fold-main fold-multi \
        fold-decoys relax md clean-interim help

help:
	@sed -n 's/^## //p' $(MAKEFILE_LIST)

## check        guards + data audit + rebuild and verify the frozen submission
check: test audit submission

## test         fast guards on the metric, the data invariants and the frozen file
test:
	$(PY) -m pytest tests/ -q

## audit        57 data-layer checks across every experiment
audit:
	$(PY) scripts/36_data_audit.py

## submission   rebuild the frozen submission and assert it is byte-identical
submission:
	$(PY) scripts/18_final_submission.py
	@echo "$(SHA)  $(FROZEN)" | sha256sum --check - \
	  || { echo "FROZEN SUBMISSION CHANGED — do not submit this"; exit 1; }

## qc           id alignment and per-complex provenance on the submission
qc:
	$(PY) scripts/21_submission_qc.py

## verify       re-run the analyses that the write-up quotes (CPU, ~20 min)
verify:
	$(PY) scripts/02_confirm_metric.py
	$(PY) scripts/20_auc_error_decomposition.py
	$(PY) scripts/37_private_projection_validation.py
	$(PY) scripts/38_e31_family_robustness.py
	$(PY) scripts/39_v_region_ablation.py

## analysis     every CPU analysis stage, in dependency order (hours)
analysis:
	$(PY) scripts/01_audit_metadata_viability.py
	$(PY) scripts/03_tcrdist_matrices.py
	$(PY) scripts/04_pgen.py
	$(PY) scripts/05_sceptr_embeddings.py
	$(PY) scripts/13_strict_loeo.py
	$(PY) scripts/15_master_comparison.py
	$(PY) scripts/19_cohort_thinning.py
	$(PY) scripts/22_fullseq_alignment_distance.py
	$(PY) scripts/23_pool_selfhit_audit.py
	$(PY) scripts/25_family_holdout_and_chain_discordance.py
	$(PY) scripts/27_remaining_sequence_cards.py
	$(PY) scripts/33_pool_negative_control.py

# ---- GPU stages. Each is hours; none is needed to rebuild the submission ----

## fold-main    Boltz-2, 2,636 complexes, ~7.4 GPU-hours
fold-main:
	$(PY) scripts/08_boltz_inputs.py && bash scripts/09_boltz_run.sh
	$(PY) scripts/10_boltz_features.py

## fold-multi   5-sample re-fold for docking consistency, ~13.2 GPU-hours
fold-multi:
	bash scripts/16_boltz_multisample.sh && $(PY) scripts/17_boltz_consistency.py

## fold-decoys  anchor-matched decoy peptides, 2,706 folds, ~7 GPU-hours
fold-decoys:
	$(PY) scripts/30_control_peptide_inputs.py && bash scripts/31_fold_decoys.sh

## relax        PDBFixer + Amber14/GBn2 minimisation of 69 complexes (CPU, ~1 h)
relax:
	$(PY) scripts/32_minimisation_robustness.py
	$(PY) scripts/32_minimisation_robustness.py --report

## md           414 ns of implicit-solvent MD, ~10 h. Needs the conda OpenMM env.
md:
	MAMBA_ROOT_PREFIX=$$HOME/micromamba PYTHONNOUSERSITE=1 \
	  micromamba run -p $(MMENV) python scripts/34_md_interface_stability.py
	$(PY) scripts/35_md_analysis.py --static
	$(PY) scripts/35_md_analysis.py --md

## clean-interim  drop derived caches, keeping external data and submissions
clean-interim:
	rm -rf data/interim/vregion data/interim/relaxed
	rm -f data/interim/md_*.csv data/interim/*.log
