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

# abTEM 4D-STEM simulation WORKER on Blythe (Warwick HPC). Normally launched via
# sim/run_sim.sh, which builds the descriptive output dir (sim_out_<tag>/ encoding
# step + phonons) and routes this job's SLURM log INTO that dir. Cluster constraints
# (L40 GPU, partition/account, CUDA) follow the infrastructure reference ../config.sh.
#
#   bash sim/run_sim.sh                              # coherent baseline params
#   SCAN_STEP=0.05 PHONONS=8 bash sim/run_sim.sh     # finer scan + TDS
#
# JOB_DIR is passed by the launcher. If this script is sbatch'd directly without it,
# output falls back to sim_out/ and refuses to clobber a finished sim (OVERWRITE=1
# forces).

# --- resolve repo root + output dir -----------------------------------------
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SIM_DIR="${REPO_DIR}/sim"
JOB_DIR="${JOB_DIR:-${REPO_DIR}/sim_out}"      # launcher passes the descriptive dir

# Tiling: when launched as an array (run_sim_tiled.sh) with N_TILES + TILES_BASE set,
# this task simulates scan x-band SLURM_ARRAY_TASK_ID of N_TILES into its own tile dir.
TILE_ARG=()
if [ -n "${N_TILES:-}" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    TILE_ARG=(--scan-tile "${SLURM_ARRAY_TASK_ID}/${N_TILES}")
    JOB_DIR="${TILES_BASE}/tile_${SLURM_ARRAY_TASK_ID}of${N_TILES}"
fi
OUT_DIR="${JOB_DIR}/01"
mkdir -p "${REPO_DIR}/logs"
if [ -e "${OUT_DIR}/data_dp.hdf5" ] && [ "${OVERWRITE:-0}" != "1" ]; then
    echo "ERROR: ${OUT_DIR}/data_dp.hdf5 already exists — refusing to overwrite a finished sim." >&2
    echo "       Use a new SIM_TAG (via run_sim.sh) or OVERWRITE=1 to force." >&2
    exit 1
fi
mkdir -p "${OUT_DIR}"

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

# Params (defaults = fast test-campaign values; launcher / --export overrides them)
SLICE_THICKNESS="${SLICE_THICKNESS:-10}"   # Å  (~7 slices through ~70 Å -> fast)
SCAN_STEP="${SCAN_STEP:-0.4}"              # Å  (50x50 = 2500 positions over 20 Å)
PHONONS="${PHONONS:-0}"                    # 0 = coherent; 8-16 for a TDS sim
PHONON_SIGMA="${PHONON_SIGMA:-0.08}"       # Å rms thermal displacement

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "REPO_DIR : ${REPO_DIR}"
echo "OUT_DIR  : ${OUT_DIR}"
echo "PYBIN    : ${PYBIN}"
echo "GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "slice=${SLICE_THICKNESS} Å  step=${SCAN_STEP} Å  phonons=${PHONONS} (sigma ${PHONON_SIGMA} Å)"
"${PYBIN}" -c "import sys,numpy,abtem; print('python:',sys.executable); print('abtem:',abtem.__version__,'numpy:',numpy.__version__)"

"${PYBIN}" "${SIM_DIR}/simulate_4dstem.py" \
    --device gpu \
    --out-dir "${OUT_DIR}" \
    --slice-thickness "${SLICE_THICKNESS}" \
    --scan-step "${SCAN_STEP}" \
    --phonons "${PHONONS}" \
    --phonon-sigma "${PHONON_SIGMA}" \
    "${TILE_ARG[@]}"
