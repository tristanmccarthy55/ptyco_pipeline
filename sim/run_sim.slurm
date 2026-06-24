#!/usr/bin/env bash
#SBATCH --job-name=abtem_sim
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/abtem_sim_%j.out
#SBATCH --error=logs/abtem_sim_%j.err

# abTEM 4D-STEM simulation on Blythe (Warwick HPC).
# Cluster constraints (L40 GPU, partition/account, CUDA + conda modules) follow
# the infrastructure reference in ../config.sh.
#
# TEST CAMPAIGN defaults: thick slices + coarse scan for a fast geometry check.
# Override per submission, e.g. for production:
#   sbatch --export=ALL,SLICE_THICKNESS=2,SCAN_STEP=0.1 sim/run_sim.slurm

# --- resolve repo root (sbatch launch dir) ----------------------------------
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SIM_DIR="${REPO_DIR}/sim"
OUT_DIR="${REPO_DIR}/sim_out/01"
mkdir -p "${REPO_DIR}/logs" "${OUT_DIR}"

# --- environment: CUDA + conda env --------------------------------------------
# CUDA toolkit (for cupy's runtime NVRTC/JIT). cupy-cuda12x bundles its own CUDA
# libs but JIT still needs a toolkit on the path.
module purge
module load CUDA/12.6.0

# Conda env with abTEM + GPU (cupy) + dask + h5py + ase + scipy. Lives on the
# big shared FS, NOT $HOME (quota). Create once (interactively):
#   module load Miniconda3/24.7.1-0
#   source "$(conda info --base)/etc/profile.d/conda.sh"
#   export CONDA_PKGS_DIRS="$SHARE/phucrh/conda/pkgs" PIP_CACHE_DIR="$SHARE/phucrh/conda/pip-cache"
#   conda create -p "$SHARE/phucrh/envs/abtem" python=3.11 -y
#   conda activate "$SHARE/phucrh/envs/abtem"
#   pip install abtem cupy-cuda12x dask h5py ase scipy numpy
#
# In the batch shell `conda activate` is unreliable (the shell hook may not load,
# and it then silently falls back to a python without our packages). Call the
# env's python by absolute path instead — deterministic, no hook needed.
CONDA_ENV="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}"
PYBIN="${CONDA_ENV}/bin/python"
if [ ! -x "${PYBIN}" ]; then
    echo "ERROR: env python not found at ${PYBIN}; create the env first." >&2
    exit 1
fi

# Test-campaign defaults (override via --export=ALL,SLICE_THICKNESS=..,SCAN_STEP=..)
SLICE_THICKNESS="${SLICE_THICKNESS:-10}"   # Å  (~7 slices through ~70 Å -> fast)
SCAN_STEP="${SCAN_STEP:-0.4}"              # Å  (50x50 = 2500 positions over 20 Å)

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "REPO_DIR : ${REPO_DIR}"
echo "OUT_DIR  : ${OUT_DIR}"
echo "PYBIN    : ${PYBIN}"
echo "GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "slice=${SLICE_THICKNESS} Å  step=${SCAN_STEP} Å"
"${PYBIN}" -c "import sys,numpy,abtem; print('python:',sys.executable); print('abtem:',abtem.__version__,'numpy:',numpy.__version__)"

"${PYBIN}" "${SIM_DIR}/simulate_4dstem.py" \
    --device gpu \
    --out-dir "${OUT_DIR}" \
    --slice-thickness "${SLICE_THICKNESS}" \
    --scan-step "${SCAN_STEP}"
