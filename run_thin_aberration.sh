#!/usr/bin/env bash
# [thin-ab] Tractable aberration-retrieval test — de-risks the thick BIN=1 experiment that
# timed out / NaN'd. A THIN PTO/STO slab (few slices) + a ROUND aberration (Cs+C5) balanced by
# defocus to a COMPACT probe that fits standard BIN=4. Few layers + a compact probe = fast AND
# well-conditioned, so blind probe retrieval has a real chance.
#
# Physics: at 70 mrad, defocus (+165 Å) cancels the round C3+C5 ray spread -> Cs=1µm + C5=1mm
# gives a ~10 Å compact-but-aberrated probe (vs ~35 Å unbalanced). Non-round terms (C56/C12)
# are OUT (defocus can't cancel them). See run_aberration_experiment.sh for the thick version.
#
# Legs (all thin, BIN=4, NL few, reg off):
#   perfect    : aberration-free compact probe (df=30 -> ~4 Å)        -> resolution reference
#   ab_known   : aberrated data, TRUE balanced probe fixed            -> ceiling (probe known)
#   ab_fitprobe: aberrated data, starts from the PERFECT probe, UPDATE -> the blind-retrieval test
# (fit-probe starts from the compact aberration-free probe, not the 35 Å df=165 disc which would
#  wrap the BIN=4 window; the recon grows/structures it toward the true aberrated probe.)
#
#   bash run_thin_aberration.sh
#   C5=5e6 THIN=2 STEP=0.4 bash run_thin_aberration.sh   # milder C5 / thinner / finer scan
# Remove later by deleting this file + the [thin-ab] blocks in simulate_4dstem.py & run_sim.slurm.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)          # probe chosen per-leg
THIN="${THIN:-3}"; SLICE=2; CONV=70; BIN=4; STEP="${STEP:-0.5}"
CS="${CS:-1e4}"; C5="${C5:-1e7}"; DF="${DF:-165}"; DF_PERF="${DF_PERF:-30}"
NL="${NL:-6}"; NITER="${NITER:-100}"; PSTART="${PSTART:-30}"; BETA="${BETA:-0.05}"

sim_job () {  # $1 tag  $2 aberrated(0/1)
    local tag="$1" ab="$2"
    local dir="${REPO_DIR}/sim_out_thin_${tag}"
    local exp="ALL,JOB_DIR=${dir},THIN_CELLS=${THIN},SLICE_THICKNESS=${SLICE},SCAN_STEP=${STEP},CONVERGENCE=${CONV},BIN_FACTOR=${BIN}"
    if [ "${ab}" = "1" ]; then exp="${exp},ABERRATED=1,PROBE_INITIAL=nominal,CS=${CS},C5=${C5},DEFOCUS=${DF}"
    else                       exp="${exp},DEFOCUS=${DF_PERF}"; fi
    sbatch --parsable --job-name="thin_sim_${tag}" --time=06:00:00 \
        --output="logs/thin_sim_${tag}_%j.out" --error="logs/thin_sim_${tag}_%j.err" \
        --export="${exp}" sim/run_sim.slurm
}
recon_job () {  # $1 name  $2 data_sim_dir  $3 probe_path  $4 PROBE_START("" = fixed)  $5 dep
    local name="$1" datadir="$2" probe="$3" pstart="$4" dep="$5"
    local rdir="${REPO_DIR}/recon_thin_${name}_NL${NL}_c${CONV}_reg0_p1"
    mkdir -p "${rdir}/01"
    for f in "${INPUTS[@]}"; do ln -sf "${datadir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${probe}" "${rdir}/01/probe_initial.mat"           # start probe (may be another sim's)
    local psx=""; [ -n "${pstart}" ] && psx=",PROBE_START=${pstart},BETA_LSQ=${BETA}"
    sbatch --parsable --job-name="thin_rec_${name}" --time=08:00:00 \
        --dependency="afterok:${dep}" \
        --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${NL}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",SAVE_EVERY=25${psx} \
        run_recon_synthetic_ML.slurm
}
PDIR="${REPO_DIR}/sim_out_thin_perfect"; ADIR="${REPO_DIR}/sim_out_thin_aberrated"
SP=$(sim_job perfect 0);   echo "sim perfect   : ${SP} -> sim_out_thin_perfect/01/ (df ${DF_PERF} Å, ~4 Å probe)"
SA=$(sim_job aberrated 1); echo "sim aberrated : ${SA} -> sim_out_thin_aberrated/01/ (Cs=${CS} C5=${C5}, df ${DF} Å, +probe_initial_true.mat)"

R1=$(recon_job perfect     "${PDIR}" "${PDIR}/01/probe_initial.mat"      ""          "${SP}")
R2=$(recon_job ab_fitprobe "${ADIR}" "${PDIR}/01/probe_initial.mat"      "${PSTART}" "${SA}:${SP}")
R3=$(recon_job ab_known    "${ADIR}" "${ADIR}/01/probe_initial_true.mat" ""          "${SA}")
echo "recon perfect ${R1} | ab_fitprobe ${R2} (probe update from iter ${PSTART}) | ab_known ${R3}"
echo
echo "Thin + compact-probe -> minutes-to-~1h each, BIN=4 (no 512-mask/parity/OOM). Compare depth"
echo "+ in-plane resolution: ab_fitprobe ~ ab_known ~ perfect => ptycho retrieved the aberration."
