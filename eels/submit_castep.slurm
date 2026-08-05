#!/bin/bash
#SBATCH -J pto_elnes
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=32              # CASTEP is CPU MPI; size from `sinfo`
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH -p compute              # <-- REAL partition name from M0 (`sinfo`); NOT the gpu queue
#
# Blythe (Warwick) core-hole ELNES driver. Work under $SHARE, never $HOME (tiny quota) --
# same rule as HPC_COMMANDS.md. Runs the expensive CASTEP SCF ONCE, then the cheap OptaDOS
# TWICE (q||c and q_perp) by swapping the .odi -> the two spectra whose difference is the
# along-beam dichroism.
#
# Submit one seed per (structure, excited-site) .cell:
#   sbatch --export=ALL,SEED=tet_Pz_Oap submit_castep.slurm
# The matching <SEED>.cell (with the SPECIES_POT core-hole block appended) and the two
# .odi files (coreloss_qc.odi, coreloss_qperp.odi) must sit in $WORKDIR alongside the .param.

set -euo pipefail
: "${SEED:?set SEED, e.g. --export=ALL,SEED=tet_Pz_Oap}"
WORKDIR="${WORKDIR:-$PWD}"
cd "$WORKDIR"
mkdir -p logs

module load CASTEP OptaDOS      # <-- REAL module names from M0 (`module avail | grep -i castep`)

# CASTEP wants <SEED>.cell + <SEED>.param present. Link the shared coreloss.param if needed.
[ -f "${SEED}.param" ] || cp templates/coreloss.param "${SEED}.param"

echo "[$(date)] CASTEP SCF: $SEED"
srun castep.mpi "$SEED"

# --- OptaDOS, both momentum directions (cheap; no CASTEP rerun) ---
for Q in qc qperp; do
  echo "[$(date)] OptaDOS core-loss: $SEED ($Q)"
  cp "coreloss_${Q}.odi" "${SEED}.odi"
  optados "$SEED"
  # collect the core-loss spectrum for this q into a q-tagged copy.
  # VERIFY the exact OptaDOS output name for your version (e.g. <SEED>_core_edge.dat)
  for f in "${SEED}"*core*.dat "${SEED}"*.agr; do
    [ -e "$f" ] && cp "$f" "${SEED}.${Q}.$(basename "$f")"
  done
done

echo "[$(date)] done: $SEED  -> pull *.qc.* and *.qperp.* back to the Mac for analyze_elnes.py"
