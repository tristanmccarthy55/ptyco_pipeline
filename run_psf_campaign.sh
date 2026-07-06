#!/usr/bin/env bash
# Generate all empirical PSFs (for analysis/atomfind). For each config: a GRID sim ->
# [dose] -> recon, chained with afterok, all configs in parallel. Grid (not a lone atom)
# because a single atom is too sparse to reconstruct (comes back ~84% noise); a 4 A grid
# at one depth converges like the real sample while each blob stays isolated in-plane and
# axially -> a clean, neighbour-free PSF (validated on Pb/NL70). After the recons land,
# extract each kernel with:  python analysis/atomfind/extract_psf.py <recon_dir> <name>
#
#   bash run_psf_campaign.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
PYBIN="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}/bin/python"
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)

# name           el  z  step  slice NL  win grid dose   (dose=0 -> coherent, no noise)
CONFIGS=(
  "Pb_NL70       Pb 37 0.15  2     70  14  4    0"     # the agent's current data (validated)
  "O_NL70        O  37 0.15  2     70  14  4    0"
  "Ti_NL70       Ti 37 0.15  2     70  14  4    0"
  "Pb_NL70_z10   Pb 10 0.15  2     70  14  4    0"     # depth series (entrance)
  "Pb_NL70_z64   Pb 64 0.15  2     70  14  4    0"     # depth series (exit)
  "Pb_rev2       Pb 37 0.05  0.5   105 8   2.5  1e8"   # reviewer-2 geometry (VERIFY first)
  "O_rev2        O  37 0.05  0.5   105 8   2.5  1e8"
  "Ti_rev2       Ti 37 0.05  0.5   105 8   2.5  1e8"
)

for cfg in "${CONFIGS[@]}"; do
    read -r name el z step slice nl win grid dose <<< "${cfg}"
    SIMDIR="sim_out_psfgrid_${name}"
    # 0.05 A over an 8 A window = 25.6k positions x 148 slices -> heavy sim + recon;
    # the coarse NL70 grids (8.8k pos, 37 slices) finish in a couple of hours.
    SIMWALL="08:00:00"; RECWALL="12:00:00"
    [ "${step}" = "0.05" ] && { SIMWALL="1-12:00:00"; RECWALL="1-12:00:00"; }

    # 1) grid sim
    SJID=$(sbatch --parsable --job-name="psf_sim_${name}" --time="${SIMWALL}" \
        --output="logs/psf_sim_${name}_%j.out" --error="logs/psf_sim_${name}_%j.err" \
        --export=ALL,JOB_DIR="${REPO_DIR}/${SIMDIR}",SINGLE_ATOM="${el}",ATOM_Z="${z}",GRID_SPACING="${grid}",SCAN_WINDOW="${win}",SCAN_STEP="${step}",SLICE_THICKNESS="${slice}" \
        sim/run_sim.slurm)

    SRC="${SIMDIR}/01"; DEP="afterok:${SJID}"
    # 2) optional dose (kernel-matched to a noisy production recon)
    if [ "${dose}" != "0" ]; then
        NOISY="${SIMDIR}_dose${dose}"
        DJID=$(sbatch --parsable --job-name="psf_dose_${name}" \
            --partition=gpu --account=physics --gres=gpu:lovelace_l40:1 --mem=32G --time=02:00:00 \
            --dependency="afterok:${SJID}" --output="logs/psf_dose_${name}_%j.out" --error="logs/psf_dose_${name}_%j.err" \
            --wrap="'${PYBIN}' '${REPO_DIR}/sim/add_poisson_noise.py' --in-dir '${SIMDIR}' --dose ${dose} --out-dir '${NOISY}'")
        SRC="${NOISY}/01"; DEP="afterok:${DJID}"
    fi

    # 3) recon (identical production params: reg off, 1 fixed probe mode, beta 0.1)
    RDIR="${REPO_DIR}/recon_psfgrid_${name}_NL${nl}_reg0_p1_b0.1"
    mkdir -p "${RDIR}/01"
    for f in "${INPUTS[@]}"; do ln -sf "${REPO_DIR}/${SRC}/${f}" "${RDIR}/01/${f}"; done
    RJID=$(sbatch --parsable --job-name="psf_rec_${name}" --time="${RECWALL}" \
        --dependency="${DEP}" --output="${RDIR}/slurm_%j.out" --error="${RDIR}/slurm_%j.err" \
        --export=ALL,NLAYERS="${nl}",SIM_BASE="${RDIR}/",REGLAYER=0,PROBE_MODES=1,NITER=100 \
        run_recon_synthetic_ML.slurm)
    echo "${name}: sim ${SJID} -> ${DEP##*:} -> recon ${RJID}  (-> recon_psfgrid_${name}_NL${nl}_reg0_p1_b0.1/)"
done
echo
echo "All PSF sim+recon chains launched in parallel. When a recon lands, pull its"
echo "  01/*step02* and run: python analysis/atomfind/extract_psf.py <recon_dir> <name>"
