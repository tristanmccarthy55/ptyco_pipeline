#!/usr/bin/env bash
#SBATCH --job-name=ptycho_synthML
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/ptycho_synthML_%j.out
#SBATCH --error=logs/ptycho_synthML_%j.err

# Multislice (depth-resolved) reconstruction of the abTEM-synthetic dataset.
# Uses the SAME sim_out/01/ as the single-slice check (no re-simulation needed).

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PTYCHO_DIR="${REPO_DIR}/ptycho"
mkdir -p "${REPO_DIR}/logs"

module purge
module load MATLAB/2024b_Update_3

SIM_BASE="${SIM_BASE:-${REPO_DIR}/sim_out/}"
# MATLAB starts in ${PTYCHO_DIR} (-sd below), so a RELATIVE SIM_BASE would resolve against ptycho/,
# not the directory sbatch was run from -- job 1272524 died that way. Anchor it to the submit dir.
case "${SIM_BASE}" in /*) ;; *) SIM_BASE="${REPO_DIR}/${SIM_BASE}";; esac
export SIM_BASE="${SIM_BASE%/}/"
export NLAYERS="${NLAYERS:-7}"          # recon depth slices (override per submission)

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "REPO_DIR : ${REPO_DIR}"
echo "SIM_BASE : ${SIM_BASE}"
echo "NLAYERS  : ${NLAYERS}"
echo "GPU      : ${CUDA_VISIBLE_DEVICES:-<unset>}"

# --- optional: write the TRIAL probe (C1 defocus search, aberration_experiment/NEXT_PHASE.md
# step 1). When PROBE_C1 is set, sim/make_probe.py builds ${SIM_BASE}01/probe_initial.mat for this
# trial's C1 (C3/C5 from the sim's aberrations.json unless PROBE_C3/PROBE_C5 override them) through
# the sim's OWN build_initial_probe on the abtem env's python -- so the probe grid is the data's by
# construction. It runs before the preflight: a writer failure fails the job here, loudly, instead
# of inside PtychoShelves' endless "File corrupt ... Retrying" loop.
if [ -n "${PROBE_C1:-}" ]; then
    CONDA_ENV="${CONDA_ENV:-$SHARE/phucrh/envs/abtem}"; PYBIN="${CONDA_ENV}/bin/python"
    [ -x "${PYBIN}" ] || { echo "PREFLIGHT FAIL: abtem env python not found at ${PYBIN}" >&2; exit 1; }
    # runtime caches off the over-quota $HOME (same as sim/run_sim.slurm)
    export MPLCONFIGDIR="${MPLCONFIGDIR:-$SHARE/phucrh/.cache/mpl}" XDG_CACHE_HOME="${XDG_CACHE_HOME:-$SHARE/phucrh/.cache/xdg}"
    export CUPY_CACHE_DIR="${CUPY_CACHE_DIR:-$SHARE/phucrh/.cache/cupy}" NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$SHARE/phucrh/.cache/numba}"
    mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}" "${CUPY_CACHE_DIR}" "${NUMBA_CACHE_DIR}"
    PW_ARG=(--sim-dir "${SIM_BASE}01" --c1="${PROBE_C1}" --out "${SIM_BASE}01/probe_initial.mat")
    [ -n "${PROBE_C3:-}" ] && PW_ARG+=(--c3="${PROBE_C3}")
    [ -n "${PROBE_C5:-}" ] && PW_ARG+=(--c5="${PROBE_C5}")
    echo "PROBE    : trial C1=${PROBE_C1} A  (sim/make_probe.py ${PW_ARG[*]})"
    "${PYBIN}" "${REPO_DIR}/sim/make_probe.py" "${PW_ARG[@]}" \
        || { echo "PREFLIGHT FAIL: probe writer failed (see above)" >&2; exit 1; }
fi

# Preflight: every input must resolve (-e follows symlinks, so a dangling link fails here). Without
# this, a missing sim_meta dies inside MATLAB with the reason in the unpacked .err, and a bad
# probe_initial.mat sends PtychoShelves' loader into an endless "File corrupt ... Retrying" loop
# that burns the whole walltime (the 2026-09-09 batch: 12 jobs). Errors go to stdout, which is packed.
missing=0
for f in sim_meta.mat data_dp.hdf5 data_position.hdf5 probe_initial.mat; do
    [ -e "${SIM_BASE}01/${f}" ] || { echo "PREFLIGHT FAIL: ${SIM_BASE}01/${f} missing or a dangling link"; missing=1; }
done
[ "${missing}" = 0 ] || { echo "aborting before MATLAB" >&2; exit 1; }

matlab -nodisplay -nosplash -nodesktop \
       -sd "${PTYCHO_DIR}" \
       -batch "run_synthetic_recon_ML"
