#!/usr/bin/env bash
#SBATCH --job-name=fusion_gpu
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=8:00:00
# No -o/-e here: SLURM will not create a missing log directory and the job would die at
# 00:00:00 with no output. run_gpu.sh passes them, pointing at fusion/runs/.

# Generic launcher for a fusion GPU experiment. run_gpu.sh sets FUSION_MOD and FUSION_ARGS.
set -uo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -d "${REPO_DIR}/fusion" ] || REPO_DIR="$(cd "${REPO_DIR}/.." && pwd)"
cd "${REPO_DIR}" || exit 1

module purge
module load CUDA/12.6.0
SHARE="${SHARE:-/springbrook/share/physics}"
PY="${CONDA_ENV:-${SHARE}/phucrh/envs/abtem}/bin/python"

echo "[$(date)] job ${SLURM_JOB_ID:-local} on $(hostname)"
echo "  experiment : fusion/${FUSION_MOD}.py ${FUSION_ARGS:-}"
echo "  python     : ${PY}"
echo "  GPU        : ${CUDA_VISIBLE_DEVICES:-<unset>}"

# shellcheck disable=SC2086
"${PY}" "fusion/${FUSION_MOD}.py" ${FUSION_ARGS:-} --device gpu || {
    echo "[$(date)] ${FUSION_MOD} FAILED (exit $?)"; exit 1; }
echo "[$(date)] done"
