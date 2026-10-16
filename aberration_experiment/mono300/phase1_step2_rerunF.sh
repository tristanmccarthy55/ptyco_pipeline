#!/usr/bin/env bash
# Phase 1 step 2, RERUN WITH THE ENGINE FIX (2026-10-16). Every mode leg that NaN'd -- at 4;2 and worse at 16;8 -- had the
# same cause: LSQML built the object's illumination from probe mode 1 alone (init_solver.m, ptycho_solver.m: aprobe2 =
# |probe{1}|^2) while the object update sums all modes, so the step was too large by ~1 / mode-1 share, locally much more.
# Fixed: the illumination is now the sum over all modes (single-mode runs unchanged). So ALL FOUR mode points -- fs 7.29 A
# (8 modes), fs 5 (7), src 0.2 and 0.4 (8) -- rerun with every leg at the original "4;2", 250 it, in new dirs (RTAG _mixF)
# on the kept step-1 sims; fs 7.29 included so every mode leg compared ran the same engine. Cancels the _mix16 legs still
# in the queue (their points cannot complete). Analyses SUFFIX _<pt>_mixF.
#   bash aberration_experiment/mono300/phase1_step2_rerunF.sh
# Submits nothing if the pulled engine lacks the fix, if the shell exports a job variable, if a dry run fails, or if a
# point's sims are not all there. Output, job ids included: logs/phase1_step2_rerunF_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; LOG="logs/phase1_step2_rerunF_$(date +%Y%m%d_%H%M%S).txt"
BASE="${SHARE:+$SHARE/$USER}"; BASE="${BASE:-$(pwd)}"
RERUN=("FOCAL_SPREAD=7.29 _fs7p29 8" "FOCAL_SPREAD=5 _fs5 7" "SOURCE_FWHM=0.2 _src0p2 8" "SOURCE_FWHM=0.4 _src0p4 8")
{
for f in ptycho/+engines/+GPU/+initialize/init_solver.m ptycho/+engines/+GPU/ptycho_solver.m; do   # the engine fix is pulled
  LC_ALL=C grep -a -q "for ll_m = 2:numel(self.probe)" "$f" || { echo "NOT SUBMITTED: ${f} lacks the all-modes illumination fix (git pull)"; exit 1; }
done
LEAK=""
for v in $( { grep -ohE '\$\{[A-Z][A-Z0-9_]*:[-+=?]' campaign/run_thin_atomfind.sh campaign/run_analysis.sh run_recon_synthetic_ML.slurm \
                | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/'
              LC_ALL=C grep -a -ohE "getenv\('[A-Z0-9_]+'\)" ptycho/run_synthetic_recon_ML.m | sed -E "s/getenv\('([A-Z0-9_]+)'\)/\1/"; } | sort -u); do
  case "$v" in SHARE|CONDA_ENV|TMPDIR|SLURM_*|CUDA_*|CUPY_CACHE_DIR|MATLAB_*|MPLCONFIGDIR|NUMBA_CACHE_DIR|XDG_CACHE_HOME) continue;; esac
  printenv "$v" >/dev/null && LEAK="${LEAK} ${v}=$(printenv "$v")"   # EXPORTED only: that is what reaches the jobs
done
[ -z "$LEAK" ] || { echo "NOT SUBMITTED: your shell exports variables the jobs read:${LEAK} -- unset them or open a new shell"; exit 1; }
T=$(bash campaign/test_coh_names.sh); echo "$T" | grep -v '^ok'
[ "$(echo "$T" | tail -1)" = "all 26 cases passed" ] || { echo "NOT SUBMITTED: the naming dry run failed (above)"; exit 1; }
PY="${CONDA_ENV:-${SHARE:-}/phucrh/envs/abtem}/bin/python"; [ -x "$PY" ] || PY=python3
T=$(PYTHON="$PY" bash campaign/test_driver_legacy.sh HEAD); echo "$T" | grep -v '^same'
[ "$(echo "$T" | tail -1)" = "all launch lines unchanged" ] || { echo "NOT SUBMITTED: the old-launch dry run failed (above)"; exit 1; }
for X in "${RERUN[@]}"; do                    # every leg's step-1 sim must still be there (the cleanup keeps them)
  read -r P S M <<<"$X"
  for m in lab Pb Ti; do
    [ -e "sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5" ] || { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${S} has no data"; exit 1; }
  done
done
IDS=$(squeue -h -u "${USER}" -o "%i %100j" 2>/dev/null | awk '$2 ~ /_mix16$/ {print $1}' | paste -sd' ' -)
if [ -n "$IDS" ]; then scancel $IDS && echo "cancelled the _mix16 legs still queued: ${IDS}"; else echo "(no _mix16 legs in the queue)"; fi
AIDS=""
echo "== reruns with the engine fix at 4;2 (RTAG _mixF)"
for X in "${RERUN[@]}"; do
  read -r P S M <<<"$X"
  OUT=$(env "$P" RECON_ONLY=1 LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        RTAG=_mixF KICK_ALL=1 KICK_START=mixed KICK_MODES="$M" KICK_NITER=250 KICK_GROUPING="4;2" RMEM=96G \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at ${S}"; exit 1; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ "$SFX" = "$S" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: suffix '${SFX}' (want ${S}) / pack '${PJ}'"; exit 1; }
  echo "$OUT" | sed -n '2p;/^pack /p'
  A=$(LABELS=armf_a080 SUFFIX="${SFX}_mixF" DEP="$PJ" GT_REGION=210 bash campaign/run_analysis.sh | head -1)
  echo "$A"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
done
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
