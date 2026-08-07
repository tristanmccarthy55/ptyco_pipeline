#!/bin/bash
#SBATCH -J coreloss
#SBATCH -o coreloss_%j.out
#SBATCH -e coreloss_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH -p compute
#
# Core-loss ELNES: CASTEP core-hole SCF (spectral/coreloss) then OptaDOS (task:core).
# Needs <SEED>.cell (core-hole, with SPECIES_POT), <SEED>.param (coreloss.param), and
# <SEED>.odi in the submit dir. OptaDOS is serial; CASTEP is MPI.
#
#   sbatch -p compute --export=ALL,SEED=tio2_O run_coreloss.slurm
#   # optional second OptaDOS pass with a different .odi (e.g. the q-perp direction for M4):
#   sbatch -p compute --export=ALL,SEED=tet_Pz_Oap,ODI2=coreloss_qperp.odi run_coreloss.slurm

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
: "${SEED:?set SEED=<name>  (needs <name>.cell + <name>.param + <name>.odi here)}"
CASTEP_MODULE="${CASTEP_MODULE:-phys/CASTEP/24.1-foss-2023b}"
OPTADOS="${OPTADOS:-/springbrook/share/physics/phucrh/optados/optados/optados.x}"

module purge 2>/dev/null || true
module load "$CASTEP_MODULE"
BIN="$(command -v castep.mpi || command -v castep.serial || command -v castep || true)"
[ -n "$BIN" ] || { echo "no castep binary after 'module load $CASTEP_MODULE'"; exit 1; }
[ -x "$OPTADOS" ] || { echo "optados not executable at $OPTADOS"; exit 1; }

echo "[$(date)] CASTEP core-hole SCF: $SEED   ($BIN)"
if [[ "$BIN" == *.mpi ]]; then srun "$BIN" "$SEED"; else "$BIN" "$SEED"; fi
grep -qiE 'Total time|finished' "${SEED}.castep" 2>/dev/null \
  || { echo "CASTEP did not finish -- see ${SEED}.castep / coreloss_${SLURM_JOB_ID}.err"; exit 1; }
echo "[$(date)] CASTEP done. Binary files written for OptaDOS:"; ls -1 "${SEED}".*bin* 2>/dev/null || true

echo "[$(date)] OptaDOS core-loss: $SEED   ($OPTADOS)"
"$OPTADOS" "$SEED" || { echo "OptaDOS failed -- ENDIANNESS? rebuild CASTEP or OptaDOS to match; see ${SEED}.odo"; exit 1; }

# optional second q-direction pass (M4): swap the .odi and re-run OptaDOS (CASTEP not repeated)
if [ -n "${ODI2:-}" ] && [ -f "$ODI2" ]; then
  echo "[$(date)] OptaDOS second pass with $ODI2"
  cp "${SEED}.odi" "${SEED}.odi.bak"; cp "$ODI2" "${SEED}.odi"
  "$OPTADOS" "$SEED" && for f in "${SEED}"*dat "${SEED}"*.agr; do [ -e "$f" ] && cp "$f" "${SEED}.qperp.$(basename "$f")"; done
  mv "${SEED}.odi.bak" "${SEED}.odi"
fi

echo "[$(date)] CORELOSS DONE -> spectrum in ${SEED}*.dat / .agr (see ${SEED}.odo for the OptaDOS log)"
