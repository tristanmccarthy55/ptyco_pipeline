#!/usr/bin/env bash
# Aberration-retrieval experiment: can ptychography FIT out the residual Cs of a corrector
# tuned for ~30 mrad but opened to 100 mrad? If yes, a modest corrector buys large-aperture
# depth resolution — the "cheaper/older scope" story (cf. Nguyen et al., Science 383, 865).
#
# Noiseless, NO phonons (perfect coherent data) so the ONLY variable is the probe aberration.
# NL70 geometry (0.15 Å step, slice 2, 20 Å window, reg off, defocus -20) — IDENTICAL to the
# clean NL70 depth run, so ~/Desktop/NL70_new_vol.npy IS the aberration-free reference (no
# need to re-run "perfect"). Aberrations (sim/simulate_4dstem.py ABERRATIONS): Cs=0.7 µm
# dominant (∝α⁴: flat to ~30 mrad, ~9 waves / 95% of χ at 100 mrad — "Cs pops back out"),
# plus C5/C56 + small parasitics. A Cs probe delocalises ~7 Å, past the default 4x-binned
# 17.5 Å window, so this runs at BIN=2 (35 Å window / 712 px, finer k-sampling — Nguyen Fig 3B).
#
# Legs (both aberrated recons at BIN=2, reg off, NL70; compare to NL70_new_vol):
#   AB + FITPROBE : NOMINAL probe start, PROBE UPDATE ON -> the experiment
#   AB + KNOWN    : TRUE probe, fixed                    -> control (data OK iff probe known)
# Read depth resolution off each vs NL70_new_vol: FITPROBE ~ NL70 => ptycho retrieved the
# aberrations. FITPROBE << NL70 but KNOWN ~ NL70 => info is there, the fit failed (raise
# PMODES / lower PSTART). NL70 (BIN=4) vs BIN=2 differ only in k-sampling, not object pixel
# (0.049 Å) — depth resolution is object-domain, so the comparison holds.
#
#   bash run_aberration_experiment.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)   # probe_initial chosen per-leg
NL=70; STEP=0.15; SLICE=2; BIN=2; NITER="${NITER:-100}"; PSTART="${PSTART:-10}"
SIMDIR="${REPO_DIR}/sim_out_aberrated"

# 1) ONE aberrated sim (BIN=2): writes data + NOMINAL probe_initial.mat + probe_initial_true.mat
SA=$(sbatch --parsable --job-name="ab_sim" --time=12:00:00 \
    --output="logs/ab_sim_%j.out" --error="logs/ab_sim_%j.err" \
    --export=ALL,JOB_DIR="${SIMDIR}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}",ABERRATED=1,PROBE_INITIAL=nominal,BIN_FACTOR="${BIN}" \
    sim/run_sim.slurm)
echo "sim aberrated (BIN=${BIN}) : ${SA} -> sim_out_aberrated/01/ (+ probe_initial_true.mat)"

recon_job () {  # $1 name  $2 probe_file  $3 PROBE_START("" = fixed)
    local name="$1" probe="$2" pstart="$3"
    local rdir="${REPO_DIR}/recon_${name}_NL${NL}_reg0_p1"
    mkdir -p "${rdir}/01"
    for f in "${INPUTS[@]}"; do ln -sf "${SIMDIR}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${SIMDIR}/01/${probe}" "${rdir}/01/probe_initial.mat"       # chosen starting probe
    local psx=""; [ -n "${pstart}" ] && psx=",PROBE_START=${pstart}"
    sbatch --parsable --job-name="ab_rec_${name}" --time=1-12:00:00 \
        --dependency="afterok:${SA}" \
        --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${NL}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}"${psx} \
        run_recon_synthetic_ML.slurm
}
R2=$(recon_job ab_fitprobe   probe_initial.mat      "${PSTART}")   # nominal start, probe update on
R3=$(recon_job ab_knownprobe probe_initial_true.mat "")            # true probe, fixed (control)
echo "recon AB+FITPROBE   : ${R2} -> recon_ab_fitprobe_NL70_reg0_p1/   (probe update from iter ${PSTART})"
echo "recon AB+KNOWNPROBE : ${R3} -> recon_ab_knownprobe_NL70_reg0_p1/ (control)"
echo
echo "Reference (aberration-free) = ~/Desktop/NL70_new_vol.npy (already have it)."
echo "Pull the newest Niter*.mat from each recon; compare depth resolution (kz plane peak /"
echo "column cross-sections) against NL70_new_vol. The recovered probe is outputs.probe of"
echo "AB+FITPROBE — compare it to sim_out_aberrated/01/probe_initial_true.mat to see what"
echo "ptycho pulled back."
