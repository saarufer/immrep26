#!/usr/bin/env bash
# Fold all 2,636 TCR-peptide-MHC complexes. ~14 s each on an RTX 5090, so about
# 10 h end to end. Boltz skips inputs whose predictions already exist unless
# --override is passed, so this script is safe to re-run after an interruption.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate
exec boltz predict data/interim/boltz/inputs \
  --out_dir data/interim/boltz/out \
  --output_format pdb \
  --write_full_pae \
  --diffusion_samples 1 \
  --recycling_steps 3 \
  --num_workers 4 \
  --seed 0
