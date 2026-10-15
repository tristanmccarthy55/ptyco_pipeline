#!/usr/bin/env bash
# Phase 1 step 2, RETRY (2026-10-15): the four mode legs that died of "error contains NaNs" (LSQML line 224) in their
# first 1-3 presolve iterations -- Pb + Ti at fs 5 A, Ti at src 0.2 A, lab at src 0.4 A -- resubmitted with EXACTLY their
# point's settings, so lab and kernels keep one reconstruction operator (PSF_KERNELS.md). The NaN is stochastic (the
# engine shuffles its batches; 8 of 12 identical-recipe legs run): LSQML takes the step length from mode 1 alone and divides
# the object step by the number of modes, which a real mixed state (mode 1 at 9-43 %) leaves ill-conditioned for sparse grid
# objects in small batches. A repeat failure shows within minutes. The original analyses of the three points (1323791 /
# 1323801 / 1323806) are cancelled and replaced by ones that wait on BOTH the old pack (its surviving legs) and the new one.
#   bash aberration_experiment/mono300/phase1_step2_retry.sh
# Submits nothing if a dead leg's job is still in the queue, if the shell exports a job variable, or if a dry run fails.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; LOG="logs/phase1_step2_retry_$(date +%Y%m%d_%H%M%S).txt"
# point setting, suffix, modes, the dead legs (their old job ids), the old pack, the old analysis
RETRY=("FOCAL_SPREAD=5 _fs5 7 Pb:1323788,Ti:1323789 1323790 1323791"
       "SOURCE_FWHM=0.2 _src0p2 8 Ti:1323799 1323800 1323801"
       "SOURCE_FWHM=0.4 _src0p4 8 lab:1323802 1323805 1323806")
{
for X in "${RETRY[@]}"; do
  read -r P S M LEGS OP OA <<<"$X"
  for lj in ${LEGS//,/ }; do
    [ -z "$(squeue -h -j "${lj#*:}" 2>/dev/null)" ] || { echo "NOT SUBMITTED: ${S} ${lj%%:*} leg (job ${lj#*:}) is still in the queue"; exit 1; }
  done
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
AIDS=""
for X in "${RETRY[@]}"; do
  read -r P S M LEGS OP OA <<<"$X"
  MS=$(echo "${LEGS//,/ }" | tr ' ' '\n' | sed 's/:.*//' | paste -sd' ' -)
  scancel "$OA" 2>/dev/null && echo "cancelled the original ${S} analysis ${OA}" || echo "(the original ${S} analysis ${OA} was not in the queue)"
  OUT=$(env "$P" MODES="$MS" RECON_ONLY=1 LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        RTAG=_mix KICK_ALL=1 KICK_START=mixed KICK_MODES="$M" KICK_NITER=250 KICK_GROUPING="4;2" RMEM=96G \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at ${S}"; exit 1; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ "$SFX" = "$S" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: suffix '${SFX}' (want ${S}) / pack '${PJ}'"; exit 1; }
  echo "$OUT" | sed -n '2p;/^pack /p'; echo "  ^ retried legs: ${MS}"
  A=$(LABELS=armf_a080 SUFFIX="${SFX}_mix" DEP="${PJ}:${OP}" GT_REGION=210 bash campaign/run_analysis.sh | head -1)
  echo "$A  (waits on packs ${PJ} and ${OP})"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
done
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
