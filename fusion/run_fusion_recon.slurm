#!/usr/bin/env bash
#SBATCH --job-name=fusion_recon
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# No -o/-e here: SLURM will not create a missing log directory and the job would die at
# 00:00:00 with no output. run_fusion_recon.sh passes --output/--error explicitly.

# Multislice HOLLOW ptychography of the fusion toy membrane. Submitted by
# fusion/run_fusion_recon.sh, which sets SIM_BASE / HSA / NLAYERS (and passes through
# NITER, GROUPING, REGLAYER, PROBE_MODES, BETA_LSQ, PROBE_START).

set -uo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -d "${REPO_DIR}/ptycho" ] || REPO_DIR="$(cd "${REPO_DIR}/.." && pwd)"
PTYCHO_DIR="${REPO_DIR}/ptycho"

module purge
module load MATLAB/2024b_Update_3

echo "[$(date)] job ${SLURM_JOB_ID:-local} on $(hostname)"
echo "  SIM_BASE : ${SIM_BASE:-<unset>}"
echo "  HSA      : ${HSA:-0.75}   (hollow semi-angle / alpha; 0 = full detector)"
echo "  NLAYERS  : ${NLAYERS:-<driver default>}"
echo "  GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"

matlab -nodisplay -nosplash -nodesktop \
       -sd "${PTYCHO_DIR}" \
       -batch "run_fusion_hollow" || {
    echo "[$(date)] HOLLOW RECON FAILED (exit $?)"; exit 1; }

echo "[$(date)] done. Read it out with:"
echo "  ~/hyperspy-bundle/bin/python fusion/make_figure.py --recon <this run>/01/*/Niter*.mat \\"
echo "      --budget fusion/runs/fusion/hollow_budget.json --out fusion/fusion_headline.png"
