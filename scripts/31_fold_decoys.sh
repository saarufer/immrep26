#!/usr/bin/env bash
# Fold each sampled TCR against its anchor-matched decoy peptides. Single
# process, one diffusion sample: the main run used 5.5 GB of 32 GB, so this
# leaves ample headroom and there is no reason to risk a second process.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate
exec boltz predict data/interim/decoy/inputs \
  --out_dir data/interim/decoy/out \
  --output_format pdb --write_full_pae \
  --diffusion_samples 1 --recycling_steps 3 --num_workers 4 --seed 0
