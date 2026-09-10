#!/usr/bin/env bash
# Submit the sign-encoding thickness sweep to the GPU queue, with the log next to the result.
# Nothing in fusion/ that touches abtem may be run on a login node: it has no CUDA driver
# (cudaErrorInsufficientDriver), so a GPU job must go through SLURM.
#
#   bash fusion/run_sign_encoding.sh
#   THICKNESS="5 10 20 40 60" N_SCAN=12 WALLTIME=12:00:00 bash fusion/run_sign_encoding.sh
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"
JOB_DIR="${REPO_DIR}/fusion/runs"
mkdir -p "${JOB_DIR}"

PARAMS=""
for v in THICKNESS N_LAT N_SCAN CONVERGENCE HSA THETA SLICE_THICKNESS OUT; do
    [ -n "${!v:-}" ] && PARAMS="${PARAMS},${v}=${!v}"
done

TIME_ARG=(); [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")
JID=$(sbatch --parsable \
    --job-name="sign_enc" \
    --output="${JOB_DIR}/sign_encoding_%j.out" \
    --error="${JOB_DIR}/sign_encoding_%j.err" \
    "${TIME_ARG[@]}" \
    --export=ALL"${PARAMS}" \
    fusion/run_sign_encoding.slurm)
echo "sign-encoding sweep -> job ${JID}"
echo "  log:    fusion/runs/sign_encoding_${JID}.out"
echo "  result: fusion/runs/sign_encoding.json"
