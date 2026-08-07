#!/bin/bash
#SBATCH -J castep
#SBATCH -o castep_%j.out
#SBATCH -e castep_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH -p compute
#
# Generic single-CASTEP-job runner for the EELS subproject: M2(a) plain SCF now, and later the
# core-hole SCF (M3/M4) once the OTFG string is locked. Needs <SEED>.cell + <SEED>.param in the
# submit dir. Auto-detects castep.mpi/.serial/castep and loads CASTEP 24.1 (foss-2023b).
#
#   sbatch -p compute --export=ALL,SEED=smoke run_castep.slurm
#   # override module if the name differs:
#   sbatch -p compute --export=ALL,SEED=smoke,CASTEP_MODULE=CASTEP/24.1-foss-2023b run_castep.slurm

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
: "${SEED:?set SEED=<name>  (needs <name>.cell + <name>.param in this dir)}"
CASTEP_MODULE="${CASTEP_MODULE:-phys/CASTEP/24.1-foss-2023b}"

module purge 2>/dev/null || true
module load "$CASTEP_MODULE"
BIN="$(command -v castep.mpi || command -v castep.serial || command -v castep || true)"
[ -n "$BIN" ] || { echo "no castep binary after 'module load $CASTEP_MODULE'"; exit 1; }
echo "[$(date)] binary: $BIN"; "$BIN" --version 2>&1 | head -6

echo "[$(date)] running CASTEP on seed '$SEED' ..."
if [[ "$BIN" == *.mpi ]]; then srun "$BIN" "$SEED"; else "$BIN" "$SEED"; fi
rc=$?

echo "[$(date)] CASTEP exit $rc -> ${SEED}.castep"
if grep -qiE 'Total time|Peak Memory|finished' "${SEED}.castep" 2>/dev/null; then
  echo "LOOKS COMPLETE"; grep -iE 'Final energy|NB est' "${SEED}.castep" | tail -3
else
  echo "did NOT finish cleanly -- inspect ${SEED}.castep and castep_${SLURM_JOB_ID}.err"
fi
exit $rc
