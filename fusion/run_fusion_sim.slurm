#!/usr/bin/env bash
#SBATCH --job-name=fusion_sim
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
# NOTE: no -o/-e here on purpose. SLURM will NOT create a missing log directory and the
# job dies at 00:00:00 with no output if you point it at one. Pass --output/--error on
# the sbatch command line (run_fusion.sh does), or let it default to the submit dir.

# 4D-STEM simulation of the fusion toy membrane, then the hollow masks + dose split.
# One simulation serves EVERY hollow semi-angle: the hole is a reconstruction-time mask.
#
#   sbatch fusion/run_fusion_sim.slurm
#   SCAN_STEP=0.2 SCAN_WINDOW=32 PHONONS=8 sbatch fusion/run_fusion_sim.slurm

set -uo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -d "${REPO_DIR}/fusion" ] || REPO_DIR="$(cd "${REPO_DIR}/.." && pwd)"
cd "${REPO_DIR}" || exit 1

module purge
module load CUDA/12.6.0

SHARE="${SHARE:-/springbrook/share/physics}"
PY="${CONDA_ENV:-${SHARE}/phucrh/envs/abtem}/bin/python"

OUT_DIR="${OUT_DIR:-${REPO_DIR}/fusion/runs/fusion}"
SCAN_WINDOW="${SCAN_WINDOW:-24}"
SCAN_STEP="${SCAN_STEP:-0.3}"
CONVERGENCE="${CONVERGENCE:-100}"
SLICE_THICKNESS="${SLICE_THICKNESS:-1.0}"
PHONONS="${PHONONS:-0}"
HSA="${HSA:-0 0.50 0.75 0.95}"

echo "[$(date)] job ${SLURM_JOB_ID:-local} on $(hostname)"
echo "  repo    : ${REPO_DIR}"
echo "  python  : ${PY}"
echo "  out     : ${OUT_DIR}"
echo "  scan    : ${SCAN_WINDOW} A window, ${SCAN_STEP} A step"
echo "  optics  : alpha ${CONVERGENCE} mrad, slice ${SLICE_THICKNESS} A, phonons ${PHONONS}"
echo "  HSA     : ${HSA}"
"${PY}" --version

# the specimen is cheap and deterministic -- rebuild it so the sim and the ground truth
# can never drift apart
"${PY}" fusion/toy_sample.py --out-dir fusion/sample || {
    echo "[$(date)] TOY BUILD FAILED"; exit 1; }

"${PY}" fusion/simulate_fusion.py \
    --out-dir "${OUT_DIR}" \
    --scan-window "${SCAN_WINDOW}" --scan-step "${SCAN_STEP}" \
    --convergence "${CONVERGENCE}" --slice-thickness "${SLICE_THICKNESS}" \
    --phonons "${PHONONS}" --hsa ${HSA} --device gpu || {
    echo "[$(date)] FUSION SIM FAILED (exit $?)"; exit 1; }

echo "[$(date)] done -> ${OUT_DIR}"
echo "next:  bash fusion/run_fusion_recon.sh 0.75          # hollow reconstruction"
