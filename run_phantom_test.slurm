#!/usr/bin/env bash
#SBATCH --job-name=orient_test
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/orient_test_%j.out
#SBATCH --error=logs/orient_test_%j.err

# Orientation (8-DOF) test: generate the thin 'F' phantom (abTEM), then reconstruct
# it 8 times — once per detector orientation (custom_data_flip) — in one job.
# Writes to <repo>/sim_phantom/ so it does NOT touch the real sim_out/.

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
mkdir -p "${REPO_DIR}/logs" "${REPO_DIR}/sim_phantom/01"

echo "Job ${SLURM_JOB_ID} on $(hostname)"

# ---- 1. phantom sim (conda python + CUDA) ----------------------------------
module purge
module load CUDA/12.6.0
CONDA_ENV="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}"
PYBIN="${CONDA_ENV}/bin/python"
if [ ! -x "${PYBIN}" ]; then echo "ERROR: env python not found at ${PYBIN}" >&2; exit 1; fi

echo "=== phantom sim ==="
"${PYBIN}" "${REPO_DIR}/sim/simulate_4dstem.py" \
    --phantom --device gpu \
    --out-dir "${REPO_DIR}/sim_phantom/01" \
    --slice-thickness 2 --scan-step 0.4

# ---- 2. 8-orientation sweep (MATLAB) ---------------------------------------
module purge
module load MATLAB/2024b_Update_3
export SIM_BASE="${REPO_DIR}/sim_phantom/"

echo "=== orientation sweep ==="
matlab -nodisplay -nosplash -nodesktop \
       -sd "${REPO_DIR}/ptycho" \
       -batch "run_orientation_sweep"
