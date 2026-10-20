#!/usr/bin/env bash
#SBATCH --job-name=af_srch
#SBATCH --partition=gpu
#SBATCH --account=physics
#SBATCH --gres=gpu:lovelace_l40:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=88G
#SBATCH --time=04:00:00
# (The resource lines above mirror run_recon_synthetic_ML.slurm's. Slurm reads #SBATCH only from the script it is handed,
# so the recon script's own header does NOT apply when this wrapper is the job: without these lines the first search
# (2026-10-04, arrays 1305039 / 1305044) ran on CPU nodes with no GPU and 0 of 98 trials finished.)
# [campaign] One trial of an outer search (campaign/run_thin_atomfind.sh SEARCH_C1F / SEARCH_C3F), run as one task of a
# job array: row SLURM_ARRAY_TASK_ID of SEARCH_MANIFEST (search_af_<leg>/trials.tsv: idx f_c1 f_c3 c1_A c3_A recon_dir)
# names the trial's recon dir and its probe's C1 / C3. The reconstruction itself is run_recon_synthetic_ML.slurm,
# unchanged: PROBE_C1 / PROBE_C3 make sim/make_probe.py write the trial probe in the job (PROBE_SCALE, e.g. the kicked
# A1, rides in the environment), and with PROBE_START unset that probe stays FIXED. campaign/select_trial.py then reads
# every trial's final error.
set -euo pipefail
: "${SEARCH_MANIFEST:?set by run_thin_atomfind.sh}" "${SLURM_ARRAY_TASK_ID:?run as a job array task}"
# no GPU -> stop now: MATLAB would otherwise fall back to its CPU solver and run into the walltime
[ -n "${CUDA_VISIBLE_DEVICES:-}" ] || { echo "search trial: no GPU allocated (CUDA_VISIBLE_DEVICES empty) -- check the #SBATCH --gres line" >&2; exit 1; }
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# columns by NAME from the manifest's header: the C1/C3 search's (idx f_c1 f_c3 c1_A c3_A recon_dir) and a staged search's
# (campaign/psearch_plan.py: idx x y c1_A c3_A c12_A phi12 recon_dir; c1_A / c3_A "-" = stage 1's answer, SEARCH_BASE_FILE)
col(){ awk -F'\t' -v i="${SLURM_ARRAY_TASK_ID}" -v k="$1" 'NR == 1 {for (j = 1; j <= NF; j++) if ($j == k) c = j; next}
                                                         $1 == i {if (c) print $c; exit}' "${SEARCH_MANIFEST}"; }
rdir=$(col recon_dir); c1=$(col c1_A); c3=$(col c3_A); c12=$(col c12_A); phi=$(col phi12)
[ -n "$rdir" ] || { echo "search trial: no row ${SLURM_ARRAY_TASK_ID} in ${SEARCH_MANIFEST}" >&2; exit 1; }
if [ "$c1" = "-" ] || [ "$c12" = "-" ]; then     # "-" = the previous stage's answer ("C1 C3 C12 phi12")
    [ -s "${SEARCH_BASE_FILE:-}" ] || { echo "search trial: a term from the previous stage, but SEARCH_BASE_FILE '${SEARCH_BASE_FILE:-}' is empty" >&2; exit 1; }
    read -r b1 b3 b12 bphi _ <"${SEARCH_BASE_FILE}"
    if [ "$c1" = "-" ]; then c1="$b1"; c3="$b3"; fi
    if [ "$c12" = "-" ]; then
        [ -n "${bphi:-}" ] || { echo "search trial: A1 from the previous stage, but ${SEARCH_BASE_FILE} has no C12 / phi12" >&2; exit 1; }
        c12="$b12"; phi="$bphi"
    fi
fi
export SIM_BASE="${REPO_DIR}/${rdir}/" PROBE_C1="$c1" PROBE_C3="$c3"
[ -n "$c12" ] && export PROBE_SET="C12=${c12}:phi12=${phi}"     # absolute A1 (applied after any PROBE_SCALE)
echo "search trial ${SLURM_ARRAY_TASK_ID}: C1 ${c1} A, C3 ${c3} A${c12:+, A1 ${c12} A @ ${phi} rad}, probe fixed (${rdir})"
exec bash "${REPO_DIR}/run_recon_synthetic_ML.slurm"
