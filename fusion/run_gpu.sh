#!/usr/bin/env bash
# Submit any fusion GPU experiment to the queue, log beside the result. The Blythe login node
# has no CUDA driver, so nothing that touches abtem may run there.
#
#   bash fusion/run_gpu.sh depth_constraint --thickness 3 5 10 20
#   bash fusion/run_gpu.sh sign_encoding --thickness 3 5 10 20 40 --n-scan 8
set -euo pipefail
[ "$#" -ge 1 ] || { echo "usage: bash fusion/run_gpu.sh <module> [args...]" >&2; exit 1; }

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"
MOD="$1"; shift
[ -f "fusion/${MOD}.py" ] || { echo "no such experiment: fusion/${MOD}.py" >&2; exit 1; }
mkdir -p fusion/runs

TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")
JID=$(sbatch --parsable \
    --job-name="fus_${MOD}" \
    --output="${REPO_DIR}/fusion/runs/${MOD}_%j.out" \
    --error="${REPO_DIR}/fusion/runs/${MOD}_%j.err" \
    "${TIME_ARG[@]}" \
    --export=ALL,FUSION_MOD="${MOD}",FUSION_ARGS="$*" \
    fusion/run_gpu.slurm)
echo "${MOD} -> job ${JID}   (log: fusion/runs/${MOD}_${JID}.out)"
