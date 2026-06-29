#!/usr/bin/env bash
# Dose series: from ONE noiseless sim, make a Poisson-noised copy at each dose and
# launch an identical recon for each (dose is the only variable) -> for a recon-vs-dose
# figure (analysis/dose_compare.py). Each recon lands in its own dose-tagged folder.
#
#   SIM_SRC=sim_out/01 DOSES="1e10 1e8 1e6 1e4" NL=70 bash run_dose_series.sh
#
# Notes:
# - The Poisson step is streamed (low memory) and runs inline here; fine for the
#   ~10 GB baseline. For the 81 GB 0.05 A set, run add_poisson_noise.py as a job instead.
# - The recons queue behind anything already using your 15-GPU allocation (e.g. a
#   running tiled sim), then run as GPUs free up.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${REPO_DIR}"
SIM_SRC="${SIM_SRC:-sim_out/01}"
SIM_DIRNAME="$(dirname "${SIM_SRC}")"              # dir CONTAINING the scan folder (01/)
DOSES="${DOSES:-1e10 1e8 1e6 1e4}"
NL="${NL:-70}"
PYBIN="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}/bin/python"

[ -e "${REPO_DIR}/${SIM_SRC}/data_dp.hdf5" ] || { echo "ERROR: no data_dp.hdf5 under ${SIM_SRC}" >&2; exit 1; }

for D in ${DOSES}; do
    NOISY_DIR="${SIM_DIRNAME}_dose${D}"           # explicit name so we know the recon source
    echo "=== dose ${D} -> ${NOISY_DIR}/ ==="
    "${PYBIN}" sim/add_poisson_noise.py \
        --in-dir "${SIM_DIRNAME}" --dose "${D}" --out-dir "${NOISY_DIR}"
    SIM_SRC="${NOISY_DIR}/01" \
        REGLAYER="${REGLAYER:-0}" PROBE_MODES="${PROBE_MODES:-1}" \
        NITER="${NITER:-120}" WALLTIME="${WALLTIME:-1-00:00:00}" \
        bash run_recon_multi.sh "${NL}"
done

echo
echo "Launched ${NL}-layer recons for doses: ${DOSES}"
echo "When done, pull each recon_*dose*/01/*step02* to ~/Desktop/dose_series/ and run:"
echo "  python analysis/dose_compare.py"
