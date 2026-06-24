#!/usr/bin/env bash
#SBATCH --job-name=ptycho_baseline
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ptycho_baseline_%j.out
#SBATCH --error=logs/ptycho_baseline_%j.err

# Baseline PrScO3 reconstruction on Blythe (Warwick HPC).
# Cluster constraints (MATLAB module, L40 GPU, partition/account) follow the
# infrastructure reference in ../config.sh.

# --- resolve repo root (directory containing this script) -------------------
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PTYCHO_DIR="${REPO_DIR}/ptycho"
mkdir -p "${REPO_DIR}/logs"

# --- environment (MATLAB module confirmed available on Blythe) --------------
module purge
module load MATLAB/2024b_Update_3

# Raw PARADIM data lives on the shared filesystem (NOT in git, NOT in $HOME,
# whose quota is tiny). Defaults to one level above the repo; override by
# exporting RAW_DATA before sbatch if it lives elsewhere on the project FS.
export RAW_DATA="${RAW_DATA:-${REPO_DIR}/../sample_data_PrScO3.mat}"

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "REPO_DIR : ${REPO_DIR}"
echo "RAW_DATA : ${RAW_DATA}"
echo "GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"

# --- run MATLAB fully headless ----------------------------------------------
# -sd sets MATLAB's start directory so addpath(pwd)/find_base_package resolve.
# -batch runs non-interactively (no desktop/JVM splash) and returns a non-zero
#  exit code if the script errors, so SLURM correctly marks the job FAILED.
matlab -nodisplay -nosplash -nodesktop \
       -sd "${PTYCHO_DIR}" \
       -batch "run_baseline_PrScO3"
