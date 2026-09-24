#!/usr/bin/env bash
# [analysis] Run analysis/analyse_sweep.py ON BLYTHE, as a CPU job next to the reconstructions, and pack only its
# small outputs (kernels, atomfind reports and figures, summary.csv). A BIN-1-class leg's h5 is ~1 GB and a sweep's
# raw data runs to ~150 GB: neither comes home. Every geometric constant is read from the legs (analyse_sweep.py).
#
#   LABELS="round_a040 ceosopt_a080" GT_REGION=210 bash campaign/run_analysis.sh     # region-box legs
#   LABELS="..." GT=/path/to/dir_with_gt_prepared.npz bash campaign/run_analysis.sh   # any other geometry
#   DEP=<job id> ...     start only once that job has ended (e.g. a sweep's pack job: afterany)
#
# GT_REGION=S uses $SHARE/$USER/gt_region<S>/gt_prepared.npz, building it in the job if it is missing
# (atomfind.make_gt_cache --thin-cells THIN --z-vacuum ZVAC --region-side S, from the simulator's own builder).
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
LABELS="${LABELS:?set LABELS to the sweep labels to analyse}"
THIN="${THIN:-5}"; ZVAC="${ZVAC:-4}"
PYBIN="${CONDA_ENV:-${SHARE:-}/phucrh/envs/abtem}/bin/python"
[ -x "${PYBIN}" ] || { echo "no python at ${PYBIN} (set CONDA_ENV)" >&2; exit 1; }
BASE="${SHARE:+$SHARE/$USER}"; BASE="${BASE:-$REPO_DIR}"
if [ -n "${GT:-}" ]; then GTDIR="$GT"; BUILD_GT=""
elif [ -n "${GT_REGION:-}" ]; then GTDIR="${BASE}/gt_region${GT_REGION}"
    BUILD_GT="[ -f '${GTDIR}/gt_prepared.npz' ] || (cd '${REPO_DIR}/analysis' && '${PYBIN}' -m atomfind.make_gt_cache --thin-cells ${THIN} --z-vacuum ${ZVAC} --region-side ${GT_REGION} --out '${GTDIR}/gt_prepared.npz')"
else echo "set GT=<dir> or GT_REGION=<side>" >&2; exit 1; fi
TS="$(date +%Y%m%d_%H%M%S)"; TAG="$(echo ${LABELS} | tr ' ' '-')"
NAME="analysis_${TAG}_${TS}"; OUT="${BASE}/${NAME}"
JOB="logs/${NAME}.sh"
cat >"${JOB}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
mkdir -p '${GTDIR}'
${BUILD_GT}
'${PYBIN}' '${REPO_DIR}/analysis/analyse_sweep.py' --root '${REPO_DIR}' --labels ${LABELS} --gt '${GTDIR}' \\
    --out '${OUT}' --z-vacuum ${ZVAC} --python '${PYBIN}' || echo "analyse_sweep: some labels did not complete (see summary.csv)"
tar czf '${OUT}.tgz' -C '${BASE}' '${NAME}'
du -h '${OUT}.tgz'
EOF
dep=(); [ -n "${DEP:-}" ] && dep=(--dependency="afterany:${DEP}")
J=$(sbatch --parsable --job-name=af_analysis --time="${ATIME:-06:00:00}" --mem="${AMEM:-64G}" --cpus-per-task=4 \
    ${dep[@]+"${dep[@]}"} --output="logs/${NAME}_%j.out" --error="logs/${NAME}_%j.err" "${JOB}")
echo "analysis ${J} -> ${OUT}.tgz  (labels: ${LABELS}; GT ${GTDIR})"
echo "scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${BASE}/analysis_${TAG}_*.tgz' ."
