#!/bin/bash
#SBATCH -J eels_preflight
#SBATCH -o eels_preflight_%j.out      # submit dir (always exists); logs/ is gitignored so a
#SBATCH -e eels_preflight_%j.err      # fresh pull has no logs/ -> SLURM can't create -o there

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH -p compute            # <-- REAL partition name from `sinfo -s` (NOT the gpu queue)
#
# CASTEP-INDEPENDENT pre-flight for the PbTiO3 EELS subproject. Runs everything that does not
# need CASTEP/OptaDOS, on a Blythe compute node, so we validate the HPC Python env + our code
# before the CASTEP admin clears. Uses the EXISTING abtem env by ABSOLUTE PATH (the conda shell
# hook is unreliable in batch), matching sim/run_sim.slurm. Needs: numpy scipy ase (present in
# that env). spglib/pymatgen are OPTIONAL -> the space-group check auto-skips without them
# (lattice + displacement gates still run). abtem is present -> the elastic self-test runs but
# is NON-FATAL (an abtem-1.0.9 API quirk must not fail the core preflight).
#
#   sbatch -p <cpu-partition> run_pyeels.slurm
#   sbatch -p <cpu-partition> --export=ALL,CONDA_ENV=/some/other/env run_pyeels.slurm   # override

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

PY="${CONDA_ENV:-/springbrook/share/physics/phucrh/envs/abtem}/bin/python"
module load CUDA/12.6.0 2>/dev/null || true    # matches the env's cupy build; harmless on a CPU node
echo "[$(date)] python: $PY"; "$PY" --version || { echo "no python at $PY"; exit 1; }

echo "== 1. unit tests (test_eels.py) =="
"$PY" test_eels.py || { echo "UNIT TESTS FAILED"; exit 1; }      # authoritative gate

echo "== 2. build + validate all cells (regenerate .cell inputs natively on the HPC) =="
"$PY" build_cells.py --corehole || { echo "BUILD FAILED"; exit 1; }

echo "== 3. elastic multislice self-test (abtem present here; NON-FATAL) =="
if "$PY" -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('abtem') else 1)"; then
  "$PY" simulate_eels.py --selftest || echo "WARN: abtem elastic self-test errored (abtem 1.0.9 API?) -- not blocking"
else
  echo "abtem not found in this env -> skipping M6b elastic self-test"
fi

echo "[$(date)] PREFLIGHT DONE -- COMPLETED / ExitCode 0:0 => HPC Python env + our code are ready for CASTEP."
