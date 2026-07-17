#!/usr/bin/env bash
# In-situ vacancy-difference kernels — analysis/atomfind/PSF_SIM_REQUEST.md REQUEST 2.
#
#   PSF_insitu(X) = angle(recon_full) - angle(recon_minus_X)
#
# the reconstructed contribution of ONE atom WITH its crystalline support. This is the
# matched filter for weak in-crystal species (O), which the isolated-atom grid CANNOT give:
# isolated O reconstructs at SNR 0.39 (below the noise floor) precisely because isolating it
# strips the neighbour support that makes it visible at all.
#
# Everything runs the NL70 production pipeline (0.15 A step, coherent, slice 2 A, 20 A window
# centred (40,20), reg off, 1 fixed probe mode, NL70) so full and vacancy differ ONLY by the
# deleted atom. Each sim prints + writes 01/vacancy_info.json recording exactly which atom
# went (element, x, y, z in the prepared-cell frame, and what its column holds).
#
# WHY full IS RECONSTRUCTED TWICE — the CONTROL (not in the request; the method needs it):
# deleting 1 atom of ~19,440 is a minute perturbation, and two independent NONLINEAR recons of
# the same data need not converge identically. full_1 - full_2 measures the reconstruction
# REPRODUCIBILITY FLOOR. If |full_1 - full_2| is not << |full_1 - vac_X|, then a "vacancy
# difference" is convergence noise, not the atom — and every kernel here is meaningless.
# Check the control FIRST. Also: difference against full_1 from THIS campaign (same code,
# params, node), not the older NL70_new_vol.npy.
#
#   bash run_vacancy_campaign.sh
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)

NL=70; STEP=0.15; SLICE=2; VZ="${VACANCY_Z:-35}"; NITER="${NITER:-100}"

# name        el  column   (el/column = "-" for the full reference)
CONFIGS=(
  "full        -   -"       # reference + control (reconstructed twice)
  "vac_O_BO    O   BO"      # 1 — the O matched filter: an O sharing its column with Ti
  "vac_O_pure  O   pure"    # 2 — an O on a pure-O column
  "vac_Ti      Ti  any"     # 3
  "vac_Pb      Pb  any"     # 4 — sanity: should ~ the isolated-Pb kernel
)

for cfg in "${CONFIGS[@]}"; do
    read -r name el col <<< "${cfg}"
    SIMDIR="${REPO_DIR}/sim_out_${name}_NL70"
    VAC_EXPORT=""
    [ "${el}" != "-" ] && VAC_EXPORT=",VACANCY=${el},VACANCY_COLUMN=${col},VACANCY_Z=${VZ}"

    SJID=$(sbatch --parsable --job-name="vac_sim_${name}" --time=08:00:00 \
        --output="logs/vac_sim_${name}_%j.out" --error="logs/vac_sim_${name}_%j.err" \
        --export=ALL,JOB_DIR="${SIMDIR}",SCAN_STEP="${STEP}",SLICE_THICKNESS="${SLICE}"${VAC_EXPORT} \
        sim/run_sim.slurm)

    REPS=1; [ "${name}" = "full" ] && REPS=2      # the full run gets a control re-recon
    for ((k=1; k<=REPS; k++)); do
        SUF=""; [ "${name}" = "full" ] && SUF="_${k}"
        RDIR="${REPO_DIR}/recon_${name}${SUF}_NL${NL}_reg0_p1"
        mkdir -p "${RDIR}/01"
        for f in "${INPUTS[@]}"; do ln -sf "${SIMDIR}/01/${f}" "${RDIR}/01/${f}"; done
        RJID=$(sbatch --parsable --job-name="vac_rec_${name}${SUF}" --time=1-00:00:00 \
            --dependency="afterok:${SJID}" \
            --output="${RDIR}/slurm_%j.out" --error="${RDIR}/slurm_%j.err" \
            --export=ALL,NLAYERS="${NL}",SIM_BASE="${RDIR}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}" \
            run_recon_synthetic_ML.slurm)
        echo "${name}${SUF}: sim ${SJID} -> recon ${RJID}  -> recon_${name}${SUF}_NL${NL}_reg0_p1/"
    done
done
echo
echo "Launched the in-situ vacancy campaign (5 sims, 6 recons, all parallel)."
echo "Deliverables (pull the newest Niter*.mat from each recon, export to complex64 .npy):"
echo "  full_1 -> the REFERENCE to difference against   full_2 -> the CONTROL"
echo "  vac_O_BO / vac_O_pure / vac_Ti / vac_Pb"
echo "Each sim's 01/vacancy_info.json says exactly which atom was deleted."
echo
echo "CHECK THE CONTROL FIRST:  |full_1 - full_2|  must be << |full_1 - vac_X|,"
echo "else the vacancy differences are convergence noise rather than the atom."
