#!/usr/bin/env bash
# Chained multislice recon at several doses, ALL DOSES IN PARALLEL — for a full
# many-iteration reconstruction that outlives the 2-day wall.
#
# Per dose: a Poisson-noise JOB (off the login node) -> seg0 (fresh: presolve + full)
# -> seg1..seg(NSEG-1), each RESTARTING from the previous segment's last checkpoint.
# Segments are chained with `afterany`, so even a wall-timed-out segment hands its
# checkpoint to the next one. Total iterations ≈ NSEG * SEG_ITERS.
#
# Usage (repo root, Blythe login node):
#   SIM_SRC=sim_out_step0.05_slice0.5_ph16s0.08/01 DOSES="1e10 1e8 1e6 1e4" \
#     NL=105 SEG_ITERS=40 NSEG=5 SAVE_EVERY=10 REGLAYER=0 WALLTIME=1-18:00:00 \
#     bash run_recon_chain.sh
#
# IMPORTANT — budget SEG_ITERS first: run one short recon (NITER=10) on one dose,
# check `sacct` Elapsed, and set SEG_ITERS so a segment COMPLETES inside WALLTIME
# (else it times out and the chain crawls). Keep SAVE_EVERY < SEG_ITERS.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"
mkdir -p logs
PYBIN="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}/bin/python"

SIM_SRC="${SIM_SRC:-sim_out/01}"
SIM_DIRNAME="$(dirname "${SIM_SRC}")"                 # dir CONTAINING the scan folder 01/
DOSES="${DOSES:-1e10 1e8 1e6 1e4}"
NL="${NL:-105}"
SEG_ITERS="${SEG_ITERS:-40}"
NSEG="${NSEG:-5}"
SAVE_EVERY="${SAVE_EVERY:-10}"
REG="${REGLAYER:-0}"
WALLTIME="${WALLTIME:-1-18:00:00}"
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)

[ -e "${REPO_DIR}/${SIM_SRC}/data_dp.hdf5" ] || { echo "ERROR: no data_dp.hdf5 under ${SIM_SRC}" >&2; exit 1; }

for D in ${DOSES}; do
    NOISY="${SIM_DIRNAME}_dose${D}"                   # the per-dose noisy dataset
    SIMREF="$(basename "${NOISY}")"; SIMREF="${SIMREF#sim_out_}"; SIMREF="${SIMREF#sim_out}"
    [ -z "${SIMREF}" ] && SIMREF="baseline"
    CHAIN="recon_${SIMREF}_NL${NL}_reg${REG}_chain"

    # 1) Poisson dose as a JOB (streamed; never on the login node)
    DJID=$(sbatch --parsable \
        --job-name="dose_${D}" \
        --partition=gpu --account=physics --gres=gpu:lovelace_l40:1 \
        --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G --time=03:00:00 \
        --output="${REPO_DIR}/logs/dose_${D}_%j.out" \
        --error="${REPO_DIR}/logs/dose_${D}_%j.err" \
        --wrap="'${PYBIN}' '${REPO_DIR}/sim/add_poisson_noise.py' --in-dir '${SIM_DIRNAME}' --dose ${D} --out-dir '${NOISY}'")
    echo "dose ${D}: job ${DJID}  -> ${NOISY}/"

    # 2) segment chain: seg0 afterok the dose job; seg k afterany seg k-1, restarting
    PREV=""; DEP="afterok:${DJID}"
    for ((k=0; k<NSEG; k++)); do
        SEG="${REPO_DIR}/${CHAIN}/seg${k}"
        mkdir -p "${SEG}/01"
        for f in "${INPUTS[@]}"; do ln -sf "${REPO_DIR}/${NOISY}/01/${f}" "${SEG}/01/${f}"; done
        RESTART_EXPORT=""; [ -n "${PREV}" ] && RESTART_EXPORT=",RESTART_DIR=${PREV}"
        SJID=$(sbatch --parsable \
            --job-name="ch_d${D}_s${k}" \
            --time="${WALLTIME}" \
            --output="${SEG}/slurm_%j.out" --error="${SEG}/slurm_%j.err" \
            --dependency="${DEP}" \
            --export=ALL,NLAYERS="${NL}",SIM_BASE="${SEG}/",REGLAYER="${REG}",NITER="${SEG_ITERS}",SAVE_EVERY="${SAVE_EVERY}"${RESTART_EXPORT} \
            run_recon_synthetic_ML.slurm)
        echo "  seg${k}: job ${SJID}  [${DEP}]  -> ${CHAIN}/seg${k}/"
        PREV="${SEG}/"; DEP="afterany:${SJID}"
    done
done

echo
echo "Launched ${NSEG}-segment chains for doses: ${DOSES}  (all doses run in parallel)"
echo "Final object per dose is the newest Niter*.mat in .../seg$((NSEG-1))/01/*/"
echo "Monitor: squeue -u \$USER"
