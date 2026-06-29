#!/usr/bin/env bash
# Tiled abTEM 4D-STEM sim for jobs too big for the 2-day wall.
#
# Splits the scan into TILES independent GPU jobs (a SLURM array, capped at
# MAXPARALLEL concurrent so we stay inside the 15-GPU limit), then a dependent merge
# job reassembles them into one dataset. Each tile is small (low memory, short wall)
# and independently retryable; together they can simulate something far larger than a
# single job could. Everything lands in one self-describing folder:
#
#   sim_out_<tag>/
#       tiles/tile_<i>of<N>/01/   per-tile data + tile_<i>.out log
#       01/                       merged data_dp / positions / probe / meta
#       merge_<jobid>.out
#
# Usage (repo root, Blythe login node):
#   SCAN_STEP=0.05 SLICE_THICKNESS=0.5 PHONONS=8 TILES=16 WALLTIME=12:00:00 \
#       bash sim/run_sim_tiled.sh
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"
SIM_DIR="${REPO_DIR}/sim"
CONDA_ENV="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}"; PYBIN="${CONDA_ENV}/bin/python"

SCAN_STEP="${SCAN_STEP:-0.05}"
SLICE_THICKNESS="${SLICE_THICKNESS:-0.5}"
PHONONS="${PHONONS:-8}"
PHONON_SIGMA="${PHONON_SIGMA:-0.08}"
TILES="${TILES:-16}"
MAXPARALLEL="${MAXPARALLEL:-15}"        # <= the gpu per-user gres/gpu limit

if [ "${PHONONS}" -gt 0 ]; then PHTAG="ph${PHONONS}s${PHONON_SIGMA}"; else PHTAG="coherent"; fi
TAG="${SIM_TAG:-step${SCAN_STEP}_slice${SLICE_THICKNESS}_${PHTAG}}"
JOB_DIR="${REPO_DIR}/sim_out_${TAG}"
TILES_BASE="${JOB_DIR}/tiles"
mkdir -p "${TILES_BASE}"

if [ -e "${JOB_DIR}/01/data_dp.hdf5" ] && [ "${OVERWRITE:-0}" != "1" ]; then
    echo "ERROR: ${JOB_DIR}/01/data_dp.hdf5 already exists — OVERWRITE=1 or a new SIM_TAG." >&2
    exit 1
fi

TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")

# 1) tile array — each task simulates one scan band on one GPU
AID=$(sbatch --parsable \
    --job-name="sim_${TAG}_tile" \
    --array=0-$((TILES-1))%"${MAXPARALLEL}" \
    --output="${TILES_BASE}/tile_%a.out" \
    --error="${TILES_BASE}/tile_%a.err" \
    "${TIME_ARG[@]}" \
    --export=ALL,N_TILES="${TILES}",TILES_BASE="${TILES_BASE}",SCAN_STEP="${SCAN_STEP}",SLICE_THICKNESS="${SLICE_THICKNESS}",PHONONS="${PHONONS}",PHONON_SIGMA="${PHONON_SIGMA}" \
    sim/run_sim.slurm)
echo "tiles : array ${AID}  (${TILES} tiles, <=${MAXPARALLEL} concurrent)  -> ${TILES_BASE}/"

# 2) merge — runs only if ALL tiles succeed; streams one tile at a time (low memory)
MID=$(sbatch --parsable \
    --job-name="sim_${TAG}_merge" \
    --partition=gpu --account=physics --gres=gpu:lovelace_l40:1 \
    --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G --time=02:00:00 \
    --dependency=afterok:"${AID}" \
    --output="${JOB_DIR}/merge_%j.out" \
    --error="${JOB_DIR}/merge_%j.err" \
    --wrap="'${PYBIN}' '${SIM_DIR}/merge_tiles.py' --tiles-dir '${TILES_BASE}' --out-dir '${JOB_DIR}'")
echo "merge : job ${MID} (afterok:${AID})  -> sim_out_${TAG}/01/"
echo
echo "next, once merged:"
echo "  python sim/add_poisson_noise.py --in-dir sim_out_${TAG} --dose 1e8"
echo "  SIM_SRC=sim_out_${TAG}_dose1e8/01 REGLAYER=0 NITER=120 WALLTIME=2-00:00:00 bash run_recon_multi.sh 105"
