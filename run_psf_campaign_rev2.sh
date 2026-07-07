#!/usr/bin/env bash
# Reviewer-2 (NL105) empirical PSFs, KERNEL-MATCHED to the 16-frozen-phonon production
# recon — proper TDS + released variable probes:
#   * 16 frozen phonons with PER-SPECIES sigma (Pb/Sr/Ti/O RT displacements), not one
#     sigma for all — Pb & O vibrate more than Ti.
#   * 0.1 A scan step: the production data was combined to 0.1 A from the 0.05 A sim, so
#     match it (also -> 80x80 = 6400 positions, so the 16-phonon sim fits ONE job, no
#     tiling needed).
#   * recon with RELEASED, mixed-state probes (PROBE_MODES>1 + PROBE_START + multimodal
#     update; optional variable/OPR probe) so it fits the noisier data instead of
#     corrupting the object with the incoherent TDS background. Saved in probe-tagged
#     folders, separate from the fixed-probe recons.
#
# Per element: ONE 16-phonon per-species grid sim -> BOTH doses (1e10, 1e8) from the same
# noiseless patterns -> a released-probe recon each. All chained with afterok.
#
#   bash run_psf_campaign_rev2.sh
#   PMODES=8 PSTART=20 VARPROBE=2 bash run_psf_campaign_rev2.sh   # also turn on OPR
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
SIM_DIR="${REPO_DIR}/sim"
PYBIN="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}/bin/python"
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)

ELEMENTS="${ELEMENTS:-Pb O Ti}"
DOSES="${DOSES:-1e10 1e8}"
WIN=8; GRID=2.5; STEP=0.1; SLICE=0.5; NL=105; PH=16; ATOMZ=37   # 8 A @ 0.1 = 6400 pos, 3x3 grid
PMODES="${PMODES:-8}"        # mixed-state probe modes (released)
PSTART="${PSTART:-20}"       # iteration the probe is released (updates from here)
VARPROBE="${VARPROBE:-0}"    # >0 = orthogonal-probe-relaxation modes (variable probe); 0 = off

PTAG="p${PMODES}_pbst${PSTART}"; [ "${VARPROBE}" != "0" ] && PTAG="${PTAG}_vp${VARPROBE}"

for el in ${ELEMENTS}; do
    SIMDIR="${REPO_DIR}/sim_out_psfgrid_${el}_rev2ph${PH}ps"      # ps = per-species sigma

    # 1) one 16-phonon per-species grid sim (0.1 A -> single job)
    SJID=$(sbatch --parsable --job-name="psf2_sim_${el}" --time=1-00:00:00 \
        --output="logs/psf2_sim_${el}_%j.out" --error="logs/psf2_sim_${el}_%j.err" \
        --export=ALL,JOB_DIR="${SIMDIR}",SINGLE_ATOM="${el}",ATOM_Z="${ATOMZ}",GRID_SPACING="${GRID}",SCAN_WINDOW="${WIN}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}",PHONONS="${PH}",PER_SPECIES_SIGMA=1 \
        sim/run_sim.slurm)
    echo "${el}: sim ${SJID}  -> $(basename "${SIMDIR}")/01/"

    # 2) both doses from the same noiseless sim -> a released-probe recon each
    for D in ${DOSES}; do
        NOISY="${SIMDIR}_dose${D}"
        DJID=$(sbatch --parsable --job-name="psf2_dose_${el}_${D}" \
            --partition=gpu --account=physics --gres=gpu:lovelace_l40:1 --mem=32G --time=03:00:00 \
            --dependency="afterok:${SJID}" \
            --output="logs/psf2_dose_${el}_${D}_%j.out" --error="logs/psf2_dose_${el}_${D}_%j.err" \
            --wrap="'${PYBIN}' '${SIM_DIR}/add_poisson_noise.py' --in-dir '${SIMDIR}' --dose ${D} --out-dir '${NOISY}'")
        RDIR="${REPO_DIR}/recon_psfgrid_${el}_rev2ph${PH}ps_dose${D}_NL${NL}_reg0_${PTAG}"
        mkdir -p "${RDIR}/01"
        for f in "${INPUTS[@]}"; do ln -sf "${NOISY}/01/${f}" "${RDIR}/01/${f}"; done
        RJID=$(sbatch --parsable --job-name="psf2_rec_${el}_${D}" --time=1-00:00:00 \
            --dependency="afterok:${DJID}" \
            --output="${RDIR}/slurm_%j.out" --error="${RDIR}/slurm_%j.err" \
            --export=ALL,NLAYERS="${NL}",SIM_BASE="${RDIR}/",REGLAYER=0,PROBE_MODES="${PMODES}",PROBE_START="${PSTART}",VARIABLE_PROBE="${VARPROBE}",NITER=100 \
            run_recon_synthetic_ML.slurm)
        echo "    dose ${D}: dose ${DJID} -> recon ${RJID}  -> recon_psfgrid_${el}_rev2ph${PH}ps_dose${D}_NL${NL}_reg0_${PTAG}/"
    done
done
echo
echo "Launched per-species 16-phonon rev2 PSF campaign: {${ELEMENTS}} x doses {${DOSES}}"
echo "  probe: ${PMODES} mixed-state modes, released at iter ${PSTART}, OPR=${VARPROBE}."
echo "Extract each: python analysis/atomfind/extract_psf.py <recon_dir> <el>_rev2_d<dose>"
