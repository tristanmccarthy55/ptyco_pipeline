#!/usr/bin/env bash
#SBATCH --job-name=ptycho_synth
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ptycho_synth_%j.out
#SBATCH --error=logs/ptycho_synth_%j.err

# Reconstruct the abTEM-synthetic dataset (hollow = 0 baseline) on Blythe.
# Reads <repo>/sim_out/01/ produced by sim/run_sim.slurm. Geometry comes from
# sim_meta.mat, so this script needs no calibration args.

# --- resolve repo root (sbatch launch dir; spool-dir safe) ------------------
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PTYCHO_DIR="${REPO_DIR}/ptycho"
mkdir -p "${REPO_DIR}/logs"

module purge
module load MATLAB/2024b_Update_3

# Where the sim wrote its per-scan folder(s). Default <repo>/sim_out/ ; the
# driver appends the scan number (01). Override with SIM_BASE if elsewhere.
export SIM_BASE="${SIM_BASE:-${REPO_DIR}/sim_out/}"

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "REPO_DIR : ${REPO_DIR}"
echo "SIM_BASE : ${SIM_BASE}"
echo "GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"

matlab -nodisplay -nosplash -nodesktop \
       -sd "${PTYCHO_DIR}" \
       -batch "run_synthetic_recon"
