#!/usr/bin/env bash
# Re-fold every complex with 5 diffusion samples, to measure docking *consistency*
# rather than confidence. Boltz runs the trunk once and samples the diffusion head,
# so 5 samples cost only ~1.36x one sample (19 s vs 14 s per complex): ~10-14 h total.
#
# Kept in a separate out_dir so the single-sample results stay intact.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate
exec boltz predict data/interim/boltz/inputs \
  --out_dir data/interim/boltz/out5 \
  --output_format pdb \
  --write_full_pae \
  --diffusion_samples 5 \
  --recycling_steps 3 \
  --num_workers 4 \
  --seed 0
