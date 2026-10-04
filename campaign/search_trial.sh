#!/usr/bin/env bash
# [campaign] One trial of an outer search (campaign/run_thin_atomfind.sh SEARCH_C1F / SEARCH_C3F), run as one task of a
# job array: row SLURM_ARRAY_TASK_ID of SEARCH_MANIFEST (search_af_<leg>/trials.tsv: idx f_c1 f_c3 c1_A c3_A recon_dir)
# names the trial's recon dir and its probe's C1 / C3. The reconstruction itself is run_recon_synthetic_ML.slurm,
# unchanged: PROBE_C1 / PROBE_C3 make sim/make_probe.py write the trial probe in the job (PROBE_SCALE, e.g. the kicked
# A1, rides in the environment), and with PROBE_START unset that probe stays FIXED. campaign/select_trial.py then reads
# every trial's final error.
set -euo pipefail
: "${SEARCH_MANIFEST:?set by run_thin_atomfind.sh}" "${SLURM_ARRAY_TASK_ID:?run as a job array task}"
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
row=$(awk -F'\t' -v i="${SLURM_ARRAY_TASK_ID}" 'NR > 1 && $1 == i' "${SEARCH_MANIFEST}")
[ -n "$row" ] || { echo "search trial: no row ${SLURM_ARRAY_TASK_ID} in ${SEARCH_MANIFEST}" >&2; exit 1; }
IFS=$'\t' read -r idx f1 f3 c1 c3 rdir <<<"$row"
export SIM_BASE="${REPO_DIR}/${rdir}/" PROBE_C1="$c1" PROBE_C3="$c3"
echo "search trial ${idx}: C1 x ${f1}, C3 x ${f3} of the kicked start -> C1 ${c1} A, C3 ${c3} A, probe fixed (${rdir})"
exec bash "${REPO_DIR}/run_recon_synthetic_ML.slurm"
