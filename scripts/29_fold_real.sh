#!/usr/bin/env bash
# Fold the crystallographic TCR-pMHC complexes through the same pipeline used
# for the competition, so the comparison in scripts/28 --compare is apples to
# apples: same model, same single-sequence mode, same chain order.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate
exec boltz predict data/interim/realval/inputs \
  --out_dir data/interim/realval/out \
  --output_format pdb --write_full_pae \
  --diffusion_samples 1 --recycling_steps 3 --num_workers 2 --seed 0
