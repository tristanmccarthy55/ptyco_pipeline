#!/usr/bin/env bash
# Launch the test campaign: abTEM sim -> synthetic recon (chained dependency).
# The recon runs only if the sim COMPLETES OK (afterok); if the sim fails it is
# auto-cancelled (DependencyNeverSatisfied).
#
# Usage (from the repo root, on a Blythe login node):
#   bash run_campaign.sh
#   SLICE_THICKNESS=2 SCAN_STEP=0.2 bash run_campaign.sh   # override sim knobs
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"            # so SLURM_SUBMIT_DIR = repo root for both jobs
mkdir -p logs

# Proper-physics defaults (override via env before the command).
SLICE_THICKNESS="${SLICE_THICKNESS:-0.5}"   # Å  (fine slicing for real diffraction)
SCAN_STEP="${SCAN_STEP:-0.4}"               # Å  (~90% overlap test scan)

echo "Submitting sim   (slice=${SLICE_THICKNESS} Å, step=${SCAN_STEP} Å) ..."
SIM_ID=$(sbatch --parsable \
    --export=ALL,SLICE_THICKNESS="${SLICE_THICKNESS}",SCAN_STEP="${SCAN_STEP}" \
    sim/run_sim.slurm)
echo "  sim job   : ${SIM_ID}"

echo "Submitting recon (afterok:${SIM_ID}) ..."
REC_ID=$(sbatch --parsable \
    --dependency="afterok:${SIM_ID}" \
    run_recon_synthetic.slurm)
echo "  recon job : ${REC_ID}  (waits for sim ${SIM_ID} to finish OK)"

echo
echo "Monitor : squeue -u \$USER"
echo "Sim log : tail -f logs/abtem_sim_${SIM_ID}.out"
echo "Rec log : tail -f logs/ptycho_synth_${REC_ID}.out"
