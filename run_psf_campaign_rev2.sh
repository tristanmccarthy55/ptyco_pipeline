#!/usr/bin/env bash
# Reviewer-2 (0.05 A / slice 0.5 / NL105) empirical PSFs, KERNEL-MATCHED to the
# 16-frozen-phonon production recon. A coherent PSF is the wrong kernel for a TDS
# recon — thermal diffuse scattering attenuates the coherent signal (Debye-Waller)
# and dumps an incoherent background the single-mode recon can't fit, so the real
# blobs are broader/weaker than a coherent single atom predicts. This pushes the PSF
# atoms through the IDENTICAL forward model: 16 phonons, same step/slice/NL/dose.
#
# Per element: ONE scan-tiled 16-phonon grid sim (16x the cost -> tiled + merged, like
# the mega-sim) -> then BOTH doses (1e10, 1e8) fanned out from the same noiseless
# patterns (sim once, dose twice) -> a recon each. All chained with afterok.
#
#   bash run_psf_campaign_rev2.sh
# Extract when a recon lands:
#   python analysis/atomfind/extract_psf.py <recon_dir> <el>_rev2_d<dose>
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
SIM_DIR="${REPO_DIR}/sim"
PYBIN="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}/bin/python"
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)

ELEMENTS="${ELEMENTS:-Pb O Ti}"
DOSES="${DOSES:-1e10 1e8}"
TILES="${TILES:-15}"; MAXPARALLEL="${MAXPARALLEL:-15}"   # <= the 15-GPU per-user limit
WIN=8; GRID=2.5; STEP=0.05; SLICE=0.5; NL=105            # 8 A window @ 0.05 = 25.6k pos, 3x3 grid
PH=16; PHSIG=0.08; ATOMZ=37

for el in ${ELEMENTS}; do
    SIMDIR="${REPO_DIR}/sim_out_psfgrid_${el}_rev2ph${PH}"
    TILES_BASE="${SIMDIR}/tiles"; mkdir -p "${TILES_BASE}"

    # 1) tiled 16-phonon grid sim (each tile = one scan band, one GPU)
    AID=$(sbatch --parsable --job-name="psf2_sim_${el}_tile" \
        --array=0-$((TILES-1))%"${MAXPARALLEL}" --time=08:00:00 \
        --output="${TILES_BASE}/tile_%a.out" --error="${TILES_BASE}/tile_%a.err" \
        --export=ALL,N_TILES="${TILES}",TILES_BASE="${TILES_BASE}",SINGLE_ATOM="${el}",ATOM_Z="${ATOMZ}",GRID_SPACING="${GRID}",SCAN_WINDOW="${WIN}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}",PHONONS="${PH}",PHONON_SIGMA="${PHSIG}" \
        sim/run_sim.slurm)

    # 2) merge (only if ALL tiles succeed; streams one tile at a time)
    MID=$(sbatch --parsable --job-name="psf2_sim_${el}_merge" \
        --partition=gpu --account=physics --gres=gpu:lovelace_l40:1 --mem=32G --time=02:00:00 \
        --dependency=afterok:"${AID}" \
        --output="${SIMDIR}/merge_%j.out" --error="${SIMDIR}/merge_%j.err" \
        --wrap="'${PYBIN}' '${SIM_DIR}/merge_tiles.py' --tiles-dir '${TILES_BASE}' --out-dir '${SIMDIR}'")
    echo "${el}: tiles ${AID} -> merge ${MID}  -> sim_out_psfgrid_${el}_rev2ph${PH}/01/"

    # 3) both doses from the same noiseless sim -> a recon each
    for D in ${DOSES}; do
        NOISY="${SIMDIR}_dose${D}"
        DJID=$(sbatch --parsable --job-name="psf2_dose_${el}_${D}" \
            --partition=gpu --account=physics --gres=gpu:lovelace_l40:1 --mem=32G --time=03:00:00 \
            --dependency=afterok:"${MID}" \
            --output="logs/psf2_dose_${el}_${D}_%j.out" --error="logs/psf2_dose_${el}_${D}_%j.err" \
            --wrap="'${PYBIN}' '${SIM_DIR}/add_poisson_noise.py' --in-dir '${SIMDIR}' --dose ${D} --out-dir '${NOISY}'")
        RDIR="${REPO_DIR}/recon_psfgrid_${el}_rev2ph${PH}_dose${D}_NL${NL}_reg0_p1_b0.1"
        mkdir -p "${RDIR}/01"
        for f in "${INPUTS[@]}"; do ln -sf "${NOISY}/01/${f}" "${RDIR}/01/${f}"; done
        RJID=$(sbatch --parsable --job-name="psf2_rec_${el}_${D}" --time=1-12:00:00 \
            --dependency=afterok:"${DJID}" \
            --output="${RDIR}/slurm_%j.out" --error="${RDIR}/slurm_%j.err" \
            --export=ALL,NLAYERS="${NL}",SIM_BASE="${RDIR}/",REGLAYER=0,PROBE_MODES=1,NITER=100 \
            run_recon_synthetic_ML.slurm)
        echo "    dose ${D}: dose ${DJID} -> recon ${RJID}  -> recon_psfgrid_${el}_rev2ph${PH}_dose${D}_NL${NL}_reg0_p1_b0.1/"
    done
done
echo
echo "Launched 16-phonon rev2 PSF campaign: {${ELEMENTS}} x doses {${DOSES}}."
echo "Sanity-check tile 0 once it lands (like the mega-sim) before trusting the merge."
echo "Extract each: python analysis/atomfind/extract_psf.py <recon_dir> <el>_rev2_d<dose>"
