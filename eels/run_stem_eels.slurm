#!/usr/bin/env bash
#SBATCH --job-name=stem_eels
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=stem_eels_%j.out
#SBATCH --error=stem_eels_%j.err
#
# Membrane STEM-EELS (simple PbTiO3 unit-cell membrane, up/down P_z domains): pixelated 4D-STEM
# DP + a full O-K EELS spectrum-image (injected CASTEP ELNES + background), into its OWN folder
# under eels/runs/stem_eels_<tag>/. REUSES sim/simulate_4dstem.py for the scattering (leaves it
# untouched). GPU L40 + CUDA/12.6.0 + the abtem env (now with gpaw). Mirrors sim/run_sim.slurm.
#
#   sbatch --export=ALL,SCAN=3x3,NTHICK=4 eels/run_stem_eels.slurm      # small, quick
#   sbatch --export=ALL,SCAN=8x8,DOMAINS=updown eels/run_stem_eels.slurm
# Params (env vars): SCAN=NxN (default 5x5), NLAT (8), NTHICK (6), DOMAINS (updown|up|down).

set -uo pipefail
REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
if   [ -f "${REPO_DIR}/simulate_stem_eels.py" ];      then cd "${REPO_DIR}"        # sbatched from eels/
elif [ -f "${REPO_DIR}/eels/simulate_stem_eels.py" ]; then cd "${REPO_DIR}/eels"  # from repo root
else echo "cannot find eels/simulate_stem_eels.py from ${REPO_DIR}"; exit 1; fi

module purge
module load CUDA/12.6.0
PY="${CONDA_ENV:-/springbrook/share/physics/phucrh/envs/abtem}/bin/python"
echo "[$(date)] python: $PY"; "$PY" --version
"$PY" -c "import abtem, gpaw, cupy; print('abtem', abtem.__version__, 'gpaw', gpaw.__version__)"

"$PY" simulate_stem_eels.py --membrane \
  --scan "${SCAN:-5x5}" --n-lat "${NLAT:-8}" --n-thick "${NTHICK:-6}" \
  --domains "${DOMAINS:-updown}" --device gpu

echo "[$(date)] done -> eels/runs/stem_eels_membrane_*  (dp.npy + eels.npz)"
