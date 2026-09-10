#!/usr/bin/env bash
# Launch one multislice HOLLOW ptychography reconstruction per hollow semi-angle, each in its
# own self-describing job folder. The simulated dataset is SHARED by symlink -- the hole is a
# reconstruction-time mask, so sweeping it costs reconstructions, not simulations.
#
#   bash fusion/run_fusion_recon.sh 0.75                 # the recommended operating point
#   bash fusion/run_fusion_recon.sh 0 0.50 0.75 0.95     # the whole sweep (4 jobs)
#   NLAYERS=24 NITER=200 bash fusion/run_fusion_recon.sh 0.75
#
# Env passed through to ptycho/run_fusion_hollow.m: NLAYERS NITER GROUPING REGLAYER
# PROBE_MODES BETA_LSQ PROBE_START WALLTIME.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

SIM_SRC="${SIM_SRC:-fusion/runs/fusion/01}"
SRC="${REPO_DIR}/${SIM_SRC}"
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)
for f in "${INPUTS[@]}"; do
    [ -e "${SRC}/${f}" ] || { echo "ERROR: missing ${SRC}/${f} — run run_fusion_sim.slurm first." >&2; exit 1; }
done
[ "$#" -ge 1 ] || { echo "usage: bash fusion/run_fusion_recon.sh <HSA> [HSA ...]   (HSA as a fraction of alpha; 0 = full detector)" >&2; exit 1; }

# depth sectioning is the whole point here: default to ~1 A layers over the reconstructed depth
BT=$(python3 -c "import json;print(json.load(open('fusion/runs/fusion/hollow_budget.json'))['beam_thickness_A'])" 2>/dev/null || echo 25.4)
NLAYERS="${NLAYERS:-$(python3 -c "print(max(1,round(${BT}/1.05)))")}"

for HSA in "$@"; do
    TAG="hsa${HSA}_NL${NLAYERS}"
    JOB_DIR="${REPO_DIR}/fusion/runs/recon_${TAG}"
    DST="${JOB_DIR}/01"
    mkdir -p "${DST}"
    for f in "${INPUTS[@]}"; do ln -sf "${SRC}/${f}" "${DST}/${f}"; done
    # the masks live with the data; symlink every one so the driver finds its own
    for m in "${SRC}"/mask_hsa*.mat; do [ -e "$m" ] && ln -sf "$m" "${DST}/$(basename "$m")"; done

    PARAMS="SIM_BASE=${JOB_DIR}/,HSA=${HSA},NLAYERS=${NLAYERS}"
    for v in NITER GROUPING REGLAYER PROBE_MODES BETA_LSQ PROBE_START; do
        [ -n "${!v:-}" ] && PARAMS="${PARAMS},${v}=${!v}"
    done

    TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")
    JID=$(sbatch --parsable \
        --job-name="fus_${TAG}" \
        --output="${JOB_DIR}/slurm_%j.out" \
        --error="${JOB_DIR}/slurm_%j.err" \
        "${TIME_ARG[@]}" \
        --export=ALL,"${PARAMS}" \
        fusion/run_fusion_recon.slurm)
    echo "HSA ${HSA}  ->  fusion/runs/recon_${TAG}/   job ${JID}   (NLAYERS=${NLAYERS})"
done

echo
echo "when they finish:"
echo "  ~/hyperspy-bundle/bin/python fusion/make_figure.py \\"
echo "      --recon fusion/runs/recon_hsa0.75_NL${NLAYERS}/01/<...>/Niter*.mat \\"
echo "      --budget fusion/runs/fusion/hollow_budget.json --out fusion/fusion_headline.png"
