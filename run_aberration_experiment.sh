#!/usr/bin/env bash
# Aberration-retrieval experiment on a TRUTHFUL JEOL ARM (hexapole C3/Cs corrector).
# Question: can ptychography FIT out a real ARM's residual aberrations when you open the
# aperture far past its ~28 mrad spec — squeezing large-aperture depth resolution from a
# modest corrector? (cf. Nguyen et al., Science 383, 865, "uncorrected" ptychography.)
#
# TRUTHFUL residual (sim/simulate_4dstem.py ABERRATIONS): C3 (Cs)=1 µm residual (perfect
# nulling isn't real) + C5=C56=1 mm UNcorrected 5th order (a C3-corrector nulls Cs, so the
# 5th order is what pops out ∝α⁶) + 0.5 nm residual astig. Flat sweet-spot ~24 mrad (= the
# ARM's real max). At 100 mrad the 5th order is ~85 waves / ~10 nm probe (bigger than the
# sample — un-simulable, and un-usable: why an ARM caps at ~28 mrad). At 70 mrad it's 24
# waves / 21 Å (fits the 70 Å box, edge wrap 0.7%) — 2.5× past spec, still a big widening.
#
# Noiseless, NO phonons; the ONLY variable is the probe. All legs: 70 mrad, BIN=1 (full 70 Å
# real-space window for the delocalised probe), NL70, reg off, step 0.3 (~4500 pos; big
# overlap given the wide probe; BIN=1 patterns are 1424² so the scan is coarsened to fit RAM).
#
# The reference is a PERFECT 70-mrad recon (NOT the 100-mrad NL70_new_vol: different NA sets
# depth resolution λ/NA² by itself — 2.0 Å at 100 vs 4.0 Å at 70). Three legs:
#   1 PERFECT-70   : no aberration, probe fixed            -> the 70-mrad reference
#   2 AB + FITPROBE: aberrated, NOMINAL start, PROBE UPDATE-> the experiment
#   3 AB + KNOWN   : aberrated, TRUE probe, fixed          -> control (data OK iff probe known)
# 2 ~ 1 => ptycho retrieved the ARM residual. 2 << 1 but 3 ~ 1 => info there, fit failed
# (raise PMODES / lower PSTART).
#
#   bash run_aberration_experiment.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)   # probe_initial chosen per-leg
NL=70; STEP="${STEP:-0.3}"; SLICE=2; BIN=1; CONV=70; NITER="${NITER:-100}"; PSTART="${PSTART:-30}"
GROUP="${GROUP:-16}"   # full-engine GPU batch: 32 needs 52 GB > 47 GB L40 at Ndp=1426; 16 -> ~29 GB
# fit-probe leg only: blind probe retrieval (nominal start vs a 24-wave aberrated truth) NaN'd
# at full res with beta 0.1 / release at 10; release later + smaller LSQ step to stabilise.
BETA="${BETA:-0.05}"
# host RAM: the BIN=1 full engine holds the 36 GB patterns + ~36 GB amplitudes; the slurm
# default 64 GB OOM-kills it (exit 9). 128 GB clears the ~80 GB peak.
MEM="${MEM:-128G}"

sim_job () {  # $1 tag  $2 ABERRATED(0/1)
    local tag="$1" ab="$2"
    local dir="${REPO_DIR}/sim_out_${tag}"
    sbatch --parsable --job-name="ab_sim_${tag}" --time=1-00:00:00 \
        --output="logs/ab_sim_${tag}_%j.out" --error="logs/ab_sim_${tag}_%j.err" \
        --export=ALL,JOB_DIR="${dir}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}",CONVERGENCE="${CONV}",BIN_FACTOR="${BIN}",ABERRATED="${ab}",PROBE_INITIAL=nominal \
        sim/run_sim.slurm
}
recon_job () {  # $1 name  $2 sim_dir  $3 probe_file  $4 PROBE_START("" = fixed)  $5 dep
    local name="$1" simdir="$2" probe="$3" pstart="$4" dep="$5"
    local rdir="${REPO_DIR}/recon_${name}_NL${NL}_c${CONV}_reg0_p1"
    mkdir -p "${rdir}/01"
    for f in "${INPUTS[@]}"; do ln -sf "${simdir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${simdir}/01/${probe}" "${rdir}/01/probe_initial.mat"       # chosen starting probe
    local psx=""; [ -n "${pstart}" ] && psx=",PROBE_START=${pstart},BETA_LSQ=${BETA}"   # probe-fit legs: stabilised
    local dep_arg=(); [ -n "${dep}" ] && dep_arg=(--dependency="afterok:${dep}")   # empty dep -> run now (RECON_ONLY)
    sbatch --parsable --job-name="ab_rec_${name}" --time=2-00:00:00 --mem="${MEM}" \
        ${dep_arg[@]+"${dep_arg[@]}"} \
        --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${NL}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",GROUPING="${GROUP}"${psx} \
        run_recon_synthetic_ML.slurm
}
# RECON_ONLY=1 reuses existing sim_out_{perfect70,aberrated70}/01/ (sims already done) and
# resubmits just the recons with no afterok dependency -> they start immediately.
if [ "${RECON_ONLY:-0}" = "1" ]; then
    SP=""; SA=""
    for t in perfect70 aberrated70; do
        [ -e "${REPO_DIR}/sim_out_${t}/01/data_dp.hdf5" ] || {
            echo "ERROR: sim_out_${t}/01/data_dp.hdf5 missing — run the sims first (unset RECON_ONLY)." >&2; exit 1; }
    done
    echo "RECON_ONLY: reusing existing sim_out_{perfect70,aberrated70}/01/ (no sims)"
else
    SP=$(sim_job perfect70 0);   echo "sim perfect-70   : ${SP} -> sim_out_perfect70/01/"
    SA=$(sim_job aberrated70 1); echo "sim aberrated-70 : ${SA} -> sim_out_aberrated70/01/ (+ probe_initial_true.mat)"
fi

R1=$(recon_job perfect70     "${REPO_DIR}/sim_out_perfect70"   probe_initial.mat      ""          "${SP}")
R2=$(recon_job ab_fitprobe   "${REPO_DIR}/sim_out_aberrated70" probe_initial.mat      "${PSTART}" "${SA}")
R3=$(recon_job ab_knownprobe "${REPO_DIR}/sim_out_aberrated70" probe_initial_true.mat ""          "${SA}")
echo "recon 1 PERFECT-70    : ${R1} -> recon_perfect70_NL70_c70_reg0_p1/"
echo "recon 2 AB+FITPROBE    : ${R2} -> recon_ab_fitprobe_NL70_c70_reg0_p1/   (probe update from iter ${PSTART})"
echo "recon 3 AB+KNOWNPROBE  : ${R3} -> recon_ab_knownprobe_NL70_c70_reg0_p1/ (control)"
echo
echo "Pull the newest Niter*.mat from each; compare depth resolution (kz plane peak / column"
echo "cross-sections) across the three. The recovered probe is outputs.probe of recon 2 —"
echo "compare it to sim_out_aberrated70/01/probe_initial_true.mat to see what ptycho retrieved."
