#!/usr/bin/env bash
#SBATCH --job-name=ptycho_synthML
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=logs/ptycho_synthML_%j.out
#SBATCH --error=logs/ptycho_synthML_%j.err

# Multislice (depth-resolved) reconstruction of the abTEM-synthetic dataset.
# Uses the SAME sim_out/01/ as the single-slice check (no re-simulation needed).

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PTYCHO_DIR="${REPO_DIR}/ptycho"
mkdir -p "${REPO_DIR}/logs"

module purge
module load MATLAB/2024b_Update_3

export SIM_BASE="${SIM_BASE:-${REPO_DIR}/sim_out/}"

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "REPO_DIR : ${REPO_DIR}"
echo "SIM_BASE : ${SIM_BASE}"
echo "GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"

matlab -nodisplay -nosplash -nodesktop \
       -sd "${PTYCHO_DIR}" \
       -batch "run_synthetic_recon_ML"
