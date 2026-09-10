#!/usr/bin/env bash
#SBATCH --job-name=sign_enc
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=8:00:00
# No -o/-e here: SLURM will not create a missing log directory and the job would die at
# 00:00:00 with no output. Pass --output/--error on the sbatch line (run_sign_encoding.sh does).

# Is the SIGN of P_z encoded in the diffraction data, and how does it grow with thickness?
# Two uniform membranes differing only in sign(delta_z), scanned at identical positions.
# NOISELESS: exact patterns are compared, nothing is ever sampled.
#
#   bash fusion/run_sign_encoding.sh                       # 3 5 10 20 40 cells
#   THICKNESS="5 10 20 40 60" N_SCAN=12 bash fusion/run_sign_encoding.sh

set -uo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -d "${REPO_DIR}/fusion" ] || REPO_DIR="$(cd "${REPO_DIR}/.." && pwd)"
cd "${REPO_DIR}" || exit 1

module purge
module load CUDA/12.6.0

SHARE="${SHARE:-/springbrook/share/physics}"
PY="${CONDA_ENV:-${SHARE}/phucrh/envs/abtem}/bin/python"

THICKNESS="${THICKNESS:-3 5 10 20 40}"
N_LAT="${N_LAT:-4}"
N_SCAN="${N_SCAN:-8}"
CONVERGENCE="${CONVERGENCE:-100}"
HSA="${HSA:-0.75}"
THETA="${THETA:-20}"
SLICE_THICKNESS="${SLICE_THICKNESS:-1.0}"
OUT="${OUT:-${REPO_DIR}/fusion/runs/sign_encoding.json}"

echo "[$(date)] job ${SLURM_JOB_ID:-local} on $(hostname)"
echo "  thickness : ${THICKNESS} cells"
echo "  specimen  : ${N_LAT}x${N_LAT} cells in-plane, delta at ${THETA} deg from the beam"
echo "  optics    : alpha ${CONVERGENCE} mrad, hole ${HSA} alpha, slice ${SLICE_THICKNESS} A"
echo "  scan      : ${N_SCAN}x${N_SCAN} positions over one unit cell"
echo "  GPU       : ${CUDA_VISIBLE_DEVICES:-<unset>}"
mkdir -p "$(dirname "${OUT}")"

"${PY}" fusion/sign_encoding.py \
    --thickness ${THICKNESS} \
    --n-lat "${N_LAT}" --n-scan "${N_SCAN}" \
    --convergence "${CONVERGENCE}" --hsa "${HSA}" --theta "${THETA}" \
    --slice-thickness "${SLICE_THICKNESS}" \
    --device gpu --out "${OUT}" || {
    echo "[$(date)] SIGN ENCODING FAILED (exit $?)"; exit 1; }

echo "[$(date)] done -> ${OUT}"
