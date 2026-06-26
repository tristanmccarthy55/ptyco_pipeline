#!/usr/bin/env bash
# Launch one multislice recon per NLAYERS value, each in its OWN self-describing
# job folder so runs never clobber each other AND you can tell them apart at a
# glance. The folder name encodes the meaningful metadata, not just Nlayers:
#
#   recon_<simref>_NL<n>_reg<r>_p<modes>_b<beta>/
#       01/                         <- symlinked sim inputs + recon outputs (Niter*.mat)
#       slurm_<jobid>.out / .err    <- this job's log lives WITH its outputs
#
#   <simref> = the sim source dir (sim_out_step0.05_ph8_dose1e5 -> step0.05_ph8_dose1e5)
#   reg/p/beta default to the driver's defaults (D=[1,0.5], 1, 0.1) when unset.
#
# Usage (from repo root on a Blythe login node):
#   bash run_recon_multi.sh 42 70                                   # baseline sim
#   REGLAYER=0 PROBE_MODES=1 bash run_recon_multi.sh 70             # reg off, 1 mode
#   SIM_SRC=sim_out_step0.05_ph8_dose1e5/01 REGLAYER=0 NITER=120 \
#       WALLTIME=2-00:00:00 bash run_recon_multi.sh 70              # off a tagged sim
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

SIM_SRC="${SIM_SRC:-sim_out/01}"
SRC="${REPO_DIR}/${SIM_SRC}"
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)
for f in "${INPUTS[@]}"; do
    [ -e "${SRC}/${f}" ] || { echo "ERROR: missing ${SRC}/${f} — run the sim first." >&2; exit 1; }
done
[ "$#" -ge 1 ] || { echo "usage: [SIM_SRC=..] [REGLAYER=..] [PROBE_MODES=..] [BETA_LSQ=..] [NITER=..] bash run_recon_multi.sh <Nlayers> [...]" >&2; exit 1; }

# --- sim reference for the folder name (strip the sim_out[_] prefix) ----------
SIMREF="$(basename "$(dirname "${SRC}")")"     # sim_out_step0.05_ph8_dose1e5 | sim_out
SIMREF="${SIMREF#sim_out_}"; SIMREF="${SIMREF#sim_out}"
[ -z "${SIMREF}" ] && SIMREF="baseline"

# --- labels for the name (value or driver-default marker) ---------------------
REG="${REGLAYER:-D}"; PM="${PROBE_MODES:-1}"; BETA="${BETA_LSQ:-0.1}"

for NL in "$@"; do
    TAG="${SIMREF}_NL${NL}_reg${REG}_p${PM}_b${BETA}"
    JOB_DIR="${REPO_DIR}/recon_${TAG}"
    DST="${JOB_DIR}/01"
    mkdir -p "${DST}"
    for f in "${INPUTS[@]}"; do ln -sf "${SRC}/${f}" "${DST}/${f}"; done   # share input data

    # pass only the params the user actually set, so unset ones use driver defaults
    PARAMS="NLAYERS=${NL},SIM_BASE=${JOB_DIR}/"
    [ -n "${REGLAYER:-}" ]    && PARAMS="${PARAMS},REGLAYER=${REGLAYER}"
    [ -n "${PROBE_MODES:-}" ] && PARAMS="${PARAMS},PROBE_MODES=${PROBE_MODES}"
    [ -n "${BETA_LSQ:-}" ]    && PARAMS="${PARAMS},BETA_LSQ=${BETA_LSQ}"
    [ -n "${NITER:-}" ]       && PARAMS="${PARAMS},NITER=${NITER}"
    [ -n "${PROBE_START:-}" ] && PARAMS="${PARAMS},PROBE_START=${PROBE_START}"

    TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")

    JID=$(sbatch --parsable \
        --job-name="ptyML_${TAG}" \
        --output="${JOB_DIR}/slurm_%j.out" \
        --error="${JOB_DIR}/slurm_%j.err" \
        "${TIME_ARG[@]}" \
        --export=ALL,"${PARAMS}" \
        run_recon_synthetic_ML.slurm)
    echo "  ${TAG}  ->  recon_${TAG}/   job ${JID}   (log: recon_${TAG}/slurm_${JID}.out)"
done

echo
echo "Monitor: squeue -u \$USER"
