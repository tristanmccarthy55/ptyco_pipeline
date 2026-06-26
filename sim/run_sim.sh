#!/usr/bin/env bash
# Launch an abTEM 4D-STEM sim into its OWN self-describing directory:
#   sim_out_step<step>_slice<slice>_<coherent|ph<N>s<sigma>>/
# containing the dataset (01/) AND this job's SLURM log+err, so everything about
# the run lives in one folder with the params in the name. Never clobbers a
# validated sim (guard below); sim_out_*/ is gitignored.
#
# Usage (from the repo root on a Blythe login node):
#   bash sim/run_sim.sh                                   # uses worker defaults
#   SCAN_STEP=0.1 SLICE_THICKNESS=2 bash sim/run_sim.sh   # production coherent
#   SCAN_STEP=0.05 PHONONS=8 bash sim/run_sim.sh          # finer scan + TDS
#   SIM_TAG=myname bash sim/run_sim.sh                    # force an explicit tag
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

SCAN_STEP="${SCAN_STEP:-0.1}"
SLICE_THICKNESS="${SLICE_THICKNESS:-2}"
PHONONS="${PHONONS:-0}"
PHONON_SIGMA="${PHONON_SIGMA:-0.08}"

if [ "${PHONONS}" -gt 0 ]; then PHTAG="ph${PHONONS}s${PHONON_SIGMA}"; else PHTAG="coherent"; fi
TAG="${SIM_TAG:-step${SCAN_STEP}_slice${SLICE_THICKNESS}_${PHTAG}}"
JOB_DIR="${REPO_DIR}/sim_out_${TAG}"
mkdir -p "${JOB_DIR}/01"

if [ -e "${JOB_DIR}/01/data_dp.hdf5" ] && [ "${OVERWRITE:-0}" != "1" ]; then
    echo "ERROR: ${JOB_DIR}/01/data_dp.hdf5 exists — refusing to overwrite a finished sim." >&2
    echo "       Use a different SIM_TAG, or OVERWRITE=1 to force." >&2
    exit 1
fi

# phonon sims are ~N_PHONONS x the coherent runtime; bump walltime past the 12 h
# worker default with e.g. WALLTIME=1-00:00:00
TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")

JID=$(sbatch --parsable \
    --job-name="sim_${TAG}" \
    --output="${JOB_DIR}/slurm_%j.out" \
    --error="${JOB_DIR}/slurm_%j.err" \
    "${TIME_ARG[@]}" \
    --export=ALL,JOB_DIR="${JOB_DIR}",SCAN_STEP="${SCAN_STEP}",SLICE_THICKNESS="${SLICE_THICKNESS}",PHONONS="${PHONONS}",PHONON_SIGMA="${PHONON_SIGMA}",OVERWRITE="${OVERWRITE:-0}" \
    sim/run_sim.slurm)
echo "sim '${TAG}'  ->  sim_out_${TAG}/   job ${JID}   (log: sim_out_${TAG}/slurm_${JID}.out)"
echo "after it finishes, add dose:   python sim/add_poisson_noise.py --in-dir sim_out_${TAG} --dose 1e5"
echo "then recon:                    SIM_SRC=sim_out_${TAG}_dose1e5/01 bash run_recon_multi.sh 70"
