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
#
# Defaults that matter here and differ from the parent driver:
#   REGLAYER    0    the depth regulariser is a low-pass in kz and would blur the very
#                    thing being measured. Only raise it (<=0.05) to rescue a divergent run.
#   PROBE_MODES 1    the simulated probe is fully coherent; extra modes just add background.
#   BETA_LSQ    0.1  Yu's value for a 30-layer hollow recon. Drop to 0.05 then 0.02 if a
#                    deep run NaNs -- slower, but the depth solve is what is at risk.
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

# KNOWN_OBJECT: start from the true object (fusion/known_object.py). The job cd's into ptycho/, so
# the path must be absolute. OBJECT_START=inf freezes it: the run then only EVALUATES the truth
# under the engine's own forward model. Both get their own folder so nothing is overwritten.
SUFFIX=""
if [ -n "${KNOWN_OBJECT:-}" ]; then
    [ -f "${KNOWN_OBJECT}" ] || { echo "ERROR: KNOWN_OBJECT not found: ${KNOWN_OBJECT}" >&2; exit 1; }
    KNOWN_OBJECT="$(cd "$(dirname "${KNOWN_OBJECT}")" && pwd)/$(basename "${KNOWN_OBJECT}")"
    SUFFIX="_known"
    case "${OBJECT_START:-1}" in inf|Inf|INF) SUFFIX="${SUFFIX}_frozen" ;; esac
fi

for HSA in "$@"; do
    TAG="hsa${HSA}_NL${NLAYERS}${SUFFIX}"
    JOB_DIR="${REPO_DIR}/fusion/runs/recon_${TAG}"
    DST="${JOB_DIR}/01"
    mkdir -p "${DST}"
    for f in "${INPUTS[@]}"; do ln -sf "${SRC}/${f}" "${DST}/${f}"; done
    # the masks live with the data; symlink every one so the driver finds its own
    for m in "${SRC}"/mask_hsa*.mat; do [ -e "$m" ] && ln -sf "$m" "${DST}/$(basename "$m")"; done

    PARAMS="SIM_BASE=${JOB_DIR}/,HSA=${HSA},NLAYERS=${NLAYERS}"
    for v in NITER GROUPING REGLAYER PROBE_MODES BETA_LSQ PROBE_START KNOWN_OBJECT OBJECT_START; do
        [ -n "${!v:-}" ] && PARAMS="${PARAMS},${v}=${!v}"
    done

    TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")
    echo "  reg=${REGLAYER:-0 (driver default: depth regulariser OFF)}  probes=${PROBE_MODES:-1}  beta_LSQ=${BETA_LSQ:-0.1}"
    [ -n "${KNOWN_OBJECT:-}" ] && echo "  start = TRUE object ${KNOWN_OBJECT}  object_change_start=${OBJECT_START:-1}"
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
echo "      --recon fusion/runs/recon_<TAG above>/01/*step02*/Niter*.mat \\"
echo "      --budget fusion/runs/fusion/hollow_budget.json --out fusion/fusion_headline.png"
