#!/usr/bin/env bash
# Empirical PSFs for the ABERRATED (JEOL ARM @70 mrad) experiment — PREP so atom-finding on
# those recons has matched kernels when we get to it. Pushes a sparse Pb GRID (the atomfind
# grid approach: a lone atom is under-constrained, a 4 A grid converges while each blob stays
# isolated -> the blob IS the system PSF) through the BYTE-IDENTICAL pipeline as
# run_aberration_experiment.sh, so each kernel matches the recon it will deconvolve.
#
# Same knobs as run_aberration_experiment.sh: 70 mrad, BIN=1, NL70, reg off, 1 probe mode,
# GROUPING=16, step 0.3. Object: Pb grid at z=37, 14 A window, 4 A spacing (16 atoms).
# Three legs, matched one-to-one to the labyrinth experiment's three recons:
#   psf_perfect70   : no aberration, probe fixed             (matches recon_perfect70)
#   psf_ab_fitprobe : aberrated, NOMINAL start, probe update (matches recon_ab_fitprobe)
#   psf_ab_known    : aberrated, TRUE probe, fixed           (matches recon_ab_knownprobe)
#
# When they land, extract each kernel (existing tool, no new code):
#   python analysis/atomfind/extract_psf.py <recon_dir> <tag>   -> psf_<tag>_vol.npy on ~/Desktop
# then point config.py's single_atom_vol at the one matching the recon being atom-found.
#
#   bash run_psf_aberrated.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)      # probe_initial chosen per-leg
EL="${EL:-Pb}"; ATOMZ=37; GRID="${GRID_SPACING:-4}"; WIN="${WIN:-14}"
NL=70; STEP="${STEP:-0.3}"; SLICE=2; BIN=1; CONV=70
GROUP="${GROUP:-16}"; NITER="${NITER:-100}"; PSTART="${PSTART:-30}"
BETA="${BETA:-0.05}"   # fit-probe leg: stabilise blind probe retrieval (see run_aberration_experiment.sh)
MEM="${MEM:-192G}"     # BIN=1 full engine needs >64 GB host RAM (36 GB patterns + amplitudes)

sim_job () {  # $1 tag  $2 ABERRATED(0/1)
    local tag="$1" ab="$2"
    local dir="${REPO_DIR}/sim_out_psfab_${tag}"
    sbatch --parsable --job-name="psfab_sim_${tag}" --time=12:00:00 \
        --output="logs/psfab_sim_${tag}_%j.out" --error="logs/psfab_sim_${tag}_%j.err" \
        --export=ALL,JOB_DIR="${dir}",SINGLE_ATOM="${EL}",ATOM_Z="${ATOMZ}",GRID_SPACING="${GRID}",SCAN_WINDOW="${WIN}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}",CONVERGENCE="${CONV}",BIN_FACTOR="${BIN}",ABERRATED="${ab}",PROBE_INITIAL=nominal \
        sim/run_sim.slurm
}
recon_job () {  # $1 name  $2 sim_dir  $3 probe_file  $4 PROBE_START("" = fixed)  $5 dep
    local name="$1" simdir="$2" probe="$3" pstart="$4" dep="$5"
    local rdir="${REPO_DIR}/recon_psfab_${name}_${EL}_NL${NL}_c${CONV}_reg0_p1"
    mkdir -p "${rdir}/01"
    for f in "${INPUTS[@]}"; do ln -sf "${simdir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${simdir}/01/${probe}" "${rdir}/01/probe_initial.mat"
    local psx=""; [ -n "${pstart}" ] && psx=",PROBE_START=${pstart},BETA_LSQ=${BETA}"   # probe-fit legs: stabilised
    local dep_arg=(); [ -n "${dep}" ] && dep_arg=(--dependency="afterok:${dep}")
    sbatch --parsable --job-name="psfab_rec_${name}" --time=1-12:00:00 --mem="${MEM}" \
        ${dep_arg[@]+"${dep_arg[@]}"} \
        --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${NL}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",GROUPING="${GROUP}"${psx} \
        run_recon_synthetic_ML.slurm
}

# RECON_ONLY=1 reuses existing sim_out_psfab_{perfect70,aberrated70}/01/ and resubmits just the
# recons (no dependency) — e.g. to redo the fit-probe leg with the stabilised PSTART/BETA.
if [ "${RECON_ONLY:-0}" = "1" ]; then
    SP=""; SA=""
    for t in perfect70 aberrated70; do
        [ -e "${REPO_DIR}/sim_out_psfab_${t}/01/data_dp.hdf5" ] || {
            echo "ERROR: sim_out_psfab_${t}/01/data_dp.hdf5 missing — run the sims first (unset RECON_ONLY)." >&2; exit 1; }
    done
    echo "RECON_ONLY: reusing existing sim_out_psfab_{perfect70,aberrated70}/01/ (no sims)"
else
    SP=$(sim_job perfect70 0);   echo "sim perfect-70 grid   : ${SP} -> sim_out_psfab_perfect70/01/"
    SA=$(sim_job aberrated70 1); echo "sim aberrated-70 grid : ${SA} -> sim_out_psfab_aberrated70/01/ (+ probe_initial_true.mat)"
fi

R1=$(recon_job perfect70     "${REPO_DIR}/sim_out_psfab_perfect70"   probe_initial.mat      ""          "${SP}")
R2=$(recon_job ab_fitprobe   "${REPO_DIR}/sim_out_psfab_aberrated70" probe_initial.mat      "${PSTART}" "${SA}")
R3=$(recon_job ab_knownprobe "${REPO_DIR}/sim_out_psfab_aberrated70" probe_initial_true.mat ""          "${SA}")
echo "recon 1 psf_perfect70   : ${R1} -> recon_psfab_perfect70_${EL}_NL70_c70_reg0_p1/"
echo "recon 2 psf_ab_fitprobe : ${R2} -> recon_psfab_ab_fitprobe_${EL}_NL70_c70_reg0_p1/   (probe update from iter ${PSTART})"
echo "recon 3 psf_ab_known    : ${R3} -> recon_psfab_ab_knownprobe_${EL}_NL70_c70_reg0_p1/"
echo
echo "Pb is the reliable single-atom kernel (isolated O/Ti fall below the recon noise floor —"
echo "for those use the in-situ vacancy method). Extract each with extract_psf.py when they land."
