#!/bin/bash
#SBATCH -J eels_preflight
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH -p compute            # <-- REAL partition name from `sinfo -s` (NOT the gpu queue)
#
# CASTEP-INDEPENDENT pre-flight for the PbTiO3 EELS subproject. Runs everything that does not
# need CASTEP/OptaDOS, on a Blythe compute node, so we validate the HPC Python env + our code
# before the CASTEP admin clears. Needs a conda env with: numpy scipy ase pymatgen spglib
# matplotlib (abtem/gpaw OPTIONAL -> the M6b elastic self-test runs only if abtem is present).
#
# One-time env (see the chat / README): create it, then point this script at it:
#   sbatch --export=ALL,CONDA_SH=$SHARE/phucrh/miniforge3/etc/profile.d/conda.sh,EELS_ENV=eels run_pyeels.slurm
# or hard-code CONDA_SH / EELS_ENV below.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

CONDA_SH="${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
EELS_ENV="${EELS_ENV:-eels}"
# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$EELS_ENV"
echo "[$(date)] python: $(which python)"; python --version

echo "== 1. unit tests (test_eels.py) =="
python test_eels.py                       # 10/10 expected; nonzero exit fails the job

echo "== 2. build + validate all cells (regenerate .cell inputs natively on the HPC) =="
python build_cells.py --corehole          # writes structures/ + structures/corehole/

echo "== 3. (optional) elastic multislice self-test -- only if abtem is in the env =="
if python -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('abtem') else 1)"; then
  python simulate_eels.py --selftest
else
  echo "abtem not installed in $EELS_ENV -> skipping M6b elastic self-test (expected for now)"
fi

echo "[$(date)] PREFLIGHT DONE -- if this job is COMPLETED/ExitCode 0:0, the HPC is ready for CASTEP."
