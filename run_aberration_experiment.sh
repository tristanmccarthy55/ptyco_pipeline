#!/usr/bin/env bash
# Aberration-retrieval experiment: can ptychography FIT out the residual aberrations of a
# Cs(C3)-corrected probe that is flat only to ~30 mrad but opened to 100 mrad? If yes, a
# modest corrector could deliver large-aperture depth resolution — cheaper/older scopes win.
#
# Noiseless, NO phonons (perfect coherent data) so the ONLY variable is the probe aberration.
# NL70 geometry (0.15 Å step, slice 2, 20 Å window, reg off) — matches the clean depth run.
# Aberrations (abTEM Cnm, in sim/simulate_4dstem.py ABERRATIONS): Cs=0, C5=60 µm, six-fold
# astig C56=40 µm (hexapole signature), + parasitic astig/coma/trefoil/quadrafoil. Flat
# (|chi|<pi/4) to ~30 mrad, ~7 waves by 100 mrad; probe delocalises 1.9 -> 4.1 Å.
#
# Three-leg comparison (all NL70, reg off, 1 probe mode, same NITER):
#   1 PERFECT      : no aberration, probe fixed          -> the reference depth resolution
#   2 AB + FITPROBE: aberrated data, NOMINAL probe start, PROBE UPDATE ON -> the experiment
#   3 AB + KNOWN   : aberrated data, TRUE probe, fixed   -> control (data is fine IF probe known)
# Read depth resolution off each: 2 ~ 1  => ptycho retrieved the aberrations. 2 << 1 but
# 3 ~ 1 => the info is there but the fit failed (try PMODES>1 / earlier PSTART).
#
#   bash run_aberration_experiment.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)   # probe_initial handled per-leg
NL=70; STEP=0.15; SLICE=2; NITER="${NITER:-100}"; PSTART="${PSTART:-10}"

sim_job () {  # $1 tag  $2 ABERRATED(0/1)
    local tag="$1" ab="$2" dir="${REPO_DIR}/sim_out_${tag}"
    sbatch --parsable --job-name="ab_sim_${tag}" --time=10:00:00 \
        --output="logs/ab_sim_${tag}_%j.out" --error="logs/ab_sim_${tag}_%j.err" \
        --export=ALL,JOB_DIR="${dir}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}",ABERRATED="${ab}",PROBE_INITIAL=nominal \
        sim/run_sim.slurm
}
recon_job () {  # $1 name  $2 sim_dir  $3 probe_file(in sim 01/)  $4 PROBE_START("" = fixed)  $5 dep
    local name="$1" simdir="$2" probe="$3" pstart="$4" dep="$5"
    local rdir="${REPO_DIR}/recon_${name}_NL${NL}_reg0_p1"
    mkdir -p "${rdir}/01"
    for f in "${INPUTS[@]}"; do ln -sf "${simdir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${simdir}/01/${probe}" "${rdir}/01/probe_initial.mat"     # chosen starting probe
    local psexport=""; [ -n "${pstart}" ] && psexport=",PROBE_START=${pstart}"
    sbatch --parsable --job-name="ab_rec_${name}" --time=1-00:00:00 \
        --dependency="afterok:${dep}" \
        --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${NL}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}"${psexport} \
        run_recon_synthetic_ML.slurm
}

SP=$(sim_job perfect 0);   echo "sim perfect   : ${SP} -> sim_out_perfect/01/"
SA=$(sim_job aberrated 1); echo "sim aberrated : ${SA} -> sim_out_aberrated/01/ (+ probe_initial_true.mat)"

R1=$(recon_job perfect      "${REPO_DIR}/sim_out_perfect"   probe_initial.mat      ""        "${SP}")
R2=$(recon_job ab_fitprobe  "${REPO_DIR}/sim_out_aberrated" probe_initial.mat      "${PSTART}" "${SA}")
R3=$(recon_job ab_knownprobe "${REPO_DIR}/sim_out_aberrated" probe_initial_true.mat ""        "${SA}")
echo "recon 1 PERFECT       : ${R1} -> recon_perfect_NL70_reg0_p1/"
echo "recon 2 AB+FITPROBE    : ${R2} -> recon_ab_fitprobe_NL70_reg0_p1/   (probe update from iter ${PSTART})"
echo "recon 3 AB+KNOWNPROBE  : ${R3} -> recon_ab_knownprobe_NL70_reg0_p1/ (control)"
echo
echo "Pull the newest Niter*.mat from each; compare depth resolution (kz plane peak /"
echo "column cross-sections). The recovered probe lives in outputs.probe of recon 2 —"
echo "compare it to sim_out_aberrated/01/probe_initial_true.mat to see what ptycho retrieved."
