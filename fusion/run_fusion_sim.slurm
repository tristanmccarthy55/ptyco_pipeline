#!/usr/bin/env bash
#SBATCH --job-name=fusion_sim
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
# NOTE: no -o/-e here on purpose. SLURM will NOT create a missing log directory and the
# job dies at 00:00:00 with no output if you point it at one. Pass --output/--error on
# the sbatch command line (run_fusion.sh does), or let it default to the submit dir.

# 4D-STEM simulation of the fusion toy membrane, then the hollow masks + dose split.
# One simulation serves EVERY hollow semi-angle: the hole is a reconstruction-time mask.
#
#   sbatch fusion/run_fusion_sim.slurm
#   SCAN_STEP=0.2 SCAN_WINDOW=32 PHONONS=8 sbatch fusion/run_fusion_sim.slurm

set -uo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -d "${REPO_DIR}/fusion" ] || REPO_DIR="$(cd "${REPO_DIR}/.." && pwd)"
cd "${REPO_DIR}" || exit 1

module purge
module load CUDA/12.6.0

SHARE="${SHARE:-/springbrook/share/physics}"
PY="${CONDA_ENV:-${SHARE}/phucrh/envs/abtem}/bin/python"

OUT_DIR="${OUT_DIR:-${REPO_DIR}/fusion/runs/fusion}"
SCAN_WINDOW="${SCAN_WINDOW:-24}"
SCAN_STEP="${SCAN_STEP:-0.3}"
CONVERGENCE="${CONVERGENCE:-100}"
SLICE_THICKNESS="${SLICE_THICKNESS:-1.0}"
PHONONS="${PHONONS:-0}"
HSA="${HSA:-0 0.50 0.75 0.95}"
BIN_FACTOR="${BIN_FACTOR:-4}"
TOY_VACUUM="${TOY_VACUUM:-2}"     # campaign uses 4: vacuum layers are where surface artifacts dump
# The specimen gets its OWN directory per leg. Changing the vacuum changes box_z, so rebuilding in
# place would silently invalidate the ground truth that the already-validated runs are scored
# against -- give a variant its own TOY_DIR and the existing sample stays untouched.
TOY_DIR="${TOY_DIR:-${REPO_DIR}/fusion/sample}"
# Probe shaping, to match the aberration campaign's depth-capable legs (campaign/round_sweep.tsv).
# a100 there is C3 -110000 A, C5 1e7 A, C1 -238 A -> d90 = 11 A at bin 2.
PROBE_ARGS=""
[ -n "${CS:-}" ]      && PROBE_ARGS="${PROBE_ARGS} --cs ${CS}"
[ -n "${C5:-}" ]      && PROBE_ARGS="${PROBE_ARGS} --c5 ${C5}"
[ -n "${DEFOCUS:-}" ] && PROBE_ARGS="${PROBE_ARGS} --defocus ${DEFOCUS}"

echo "[$(date)] job ${SLURM_JOB_ID:-local} on $(hostname)"
echo "  repo    : ${REPO_DIR}"
echo "  python  : ${PY}"
echo "  out     : ${OUT_DIR}"
echo "  scan    : ${SCAN_WINDOW} A window, ${SCAN_STEP} A step"
echo "  optics  : alpha ${CONVERGENCE} mrad, slice ${SLICE_THICKNESS} A, bin ${BIN_FACTOR}, phonons ${PHONONS}"
echo "  probe   : ${PROBE_ARGS:-default (aberration-free, -20 A overfocus)}"
echo "  vacuum  : ${TOY_VACUUM} A each side   (specimen in ${TOY_DIR})"
echo "  HSA     : ${HSA}"
"${PY}" --version

# the specimen is cheap and deterministic -- rebuild it so the sim and the ground truth
# can never drift apart
"${PY}" fusion/toy_sample.py --out-dir "${TOY_DIR}" --vacuum "${TOY_VACUUM}" || {
    echo "[$(date)] TOY BUILD FAILED"; exit 1; }

"${PY}" fusion/simulate_fusion.py \
    --out-dir "${OUT_DIR}" \
    --vasp "${TOY_DIR}/toy_membrane.vasp" --truth "${TOY_DIR}/toy_truth.npz" \
    --scan-window "${SCAN_WINDOW}" --scan-step "${SCAN_STEP}" \
    --convergence "${CONVERGENCE}" --slice-thickness "${SLICE_THICKNESS}" \
    --phonons "${PHONONS}" --hsa ${HSA} --bin-factor "${BIN_FACTOR}" \
    ${PROBE_ARGS} --device gpu || {
    echo "[$(date)] FUSION SIM FAILED (exit $?)"; exit 1; }

echo "[$(date)] done -> ${OUT_DIR}"
echo "next:  bash fusion/run_fusion_recon.sh 0.75          # hollow reconstruction"
