#!/usr/bin/env bash
# Launch one multislice recon per NLAYERS value, each in its OWN output dir, so
# they run concurrently without clobbering each other's prepared-data / analysis /
# O_phase_roi outputs.
#
# Each run gets a private base_path (recon_NL<n>/) with the shared sim_out/01 input
# data SYMLINKED in — PtychoShelves writes all outputs under base_path, so separate
# base_path == separate outputs. No driver changes needed.
#
# Usage (from the repo root, on a Blythe login node):
#   bash run_recon_multi.sh 41 82            # two jobs, 41- and 82-layer, at once
#   SIM_SRC=sim_out/01 bash run_recon_multi.sh 41 74 82
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"
mkdir -p logs

SRC="${REPO_DIR}/${SIM_SRC:-sim_out/01}"     # shared input (the sim output)
INPUTS=(data_dp.hdf5 data_position.hdf5 probe_initial.mat sim_meta.mat)
for f in "${INPUTS[@]}"; do
    [ -e "${SRC}/${f}" ] || { echo "ERROR: missing ${SRC}/${f} — run the sim first." >&2; exit 1; }
done

if [ "$#" -eq 0 ]; then echo "usage: bash run_recon_multi.sh <Nlayers> [<Nlayers> ...]" >&2; exit 1; fi

for NL in "$@"; do
    TAG="NL${NL}"
    DST="${REPO_DIR}/recon_${TAG}/01"
    mkdir -p "${DST}"
    for f in "${INPUTS[@]}"; do ln -sf "${SRC}/${f}" "${DST}/${f}"; done   # symlink shared data

    TIME_ARG=()
    [ -n "${WALLTIME:-}" ] && TIME_ARG=(--time="${WALLTIME}")   # e.g. WALLTIME=1-18:00:00 for 82 layers
    JID=$(sbatch --parsable \
        --job-name="ptyML_${TAG}" \
        "${TIME_ARG[@]}" \
        --export=ALL,NLAYERS="${NL}",SIM_BASE="${REPO_DIR}/recon_${TAG}/" \
        run_recon_synthetic_ML.slurm)
    echo "  NLAYERS=${NL}  ->  recon_${TAG}/   job ${JID}   (log: logs/ptycho_synthML_${JID}.out)"
done

echo
echo "Monitor: squeue -u \$USER"
