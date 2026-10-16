#!/usr/bin/env bash
# Phase 1 step 2, RERUN at a larger grouping (the user's go, 2026-10-16). Three mode legs NaN at their start twice at "4;2"
# (Pb fs5, Ti src0p2, lab src0p4: close to deterministic -- LSQML fits its step on mode 1 alone, ill-conditioned when a
# real mixed state leaves mode 1 at 9-43 %). So every leg of the three incomplete points -- fs 5 A (7 modes), src 0.2 and
# 0.4 A (8) -- reruns at "16;8" (4x the positions per batch; 4;2 used ~6-7 GB of 47.7), one operator per point, in NEW
# dirs (RTAG _mix16) so the 4;2 legs stay for comparison; on the kept step-1 sims (RECON_ONLY). Analyses SUFFIX _<pt>_mix16.
# Also: cancels the 4;2 Ti fs5 retry (1325697: its point can never complete at 4;2), and reads the two finished 4;2 lab legs
# (fs5, src0p2) at once with their step-1 FIXED-PROBE kernels -- a labelled fallback (tag _mix_kern-armf_a080_<pt>_fixed),
# provisional until the matched _mix16 kernels land.
#   bash aberration_experiment/mono300/phase1_step2_rerun16.sh
# Submits nothing if the shell exports a job variable, if a dry run fails, if a rerun point's sims are not all there, or
# if a fallback kernel is missing. Output, job ids included: logs/phase1_step2_rerun16_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; LOG="logs/phase1_step2_rerun16_$(date +%Y%m%d_%H%M%S).txt"
BASE="${SHARE:+$SHARE/$USER}"; BASE="${BASE:-$(pwd)}"
RERUN=("FOCAL_SPREAD=5 _fs5 7" "SOURCE_FWHM=0.2 _src0p2 8" "SOURCE_FWHM=0.4 _src0p4 8")
FALLBACK=(_fs5 _src0p2)                       # finished 4;2 lab legs, read with the step-1 fixed-probe kernels
CANCEL=1325697                                # the 4;2 Ti fs5 retry
{
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
KDIR=()                                       # per FALLBACK entry: the step-1 analysis whose kernels it reads
for i in "${!FALLBACK[@]}"; do
  S="${FALLBACK[$i]}"; KDIR[$i]=$(ls -d "${BASE}"/analysis_armf_a080${S}_fixed_2026*/psf 2>/dev/null | sort | tail -1)
  for el in Pb Ti; do
    [ -f "${KDIR[$i]}/psf_${el}_armf_a080${S}_fixed_vol.npy" ] \
      || { echo "NOT SUBMITTED: no ${el} kernel psf_${el}_armf_a080${S}_fixed_vol.npy under ${BASE}/analysis_armf_a080${S}_fixed_*/psf"; exit 1; }
  done
done
if [ -n "$(squeue -h -j "$CANCEL" 2>/dev/null)" ]; then scancel "$CANCEL" && echo "cancelled the 4;2 Ti fs5 retry ${CANCEL}"
else echo "(the 4;2 Ti fs5 retry ${CANCEL} was not in the queue)"; fi
AIDS=""
echo "== reruns at 16;8 (RTAG _mix16)"
for X in "${RERUN[@]}"; do
  read -r P S M <<<"$X"
  OUT=$(env "$P" RECON_ONLY=1 LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        RTAG=_mix16 KICK_ALL=1 KICK_START=mixed KICK_MODES="$M" KICK_NITER=250 KICK_GROUPING="16;8" RMEM=96G \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at ${S}"; exit 1; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ "$SFX" = "$S" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: suffix '${SFX}' (want ${S}) / pack '${PJ}'"; exit 1; }
  echo "$OUT" | sed -n '2p;/^pack /p'
  A=$(LABELS=armf_a080 SUFFIX="${SFX}_mix16" DEP="$PJ" GT_REGION=210 bash campaign/run_analysis.sh | head -1)
  echo "$A"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
done
echo "== provisional: the finished 4;2 lab legs with their step-1 fixed-probe kernels"
for i in "${!FALLBACK[@]}"; do
  S="${FALLBACK[$i]}"
  A=$(LABELS=armf_a080 SUFFIX="${S}_mix" PSF_DIR="${KDIR[$i]}" PSF_TAG="armf_a080${S}_fixed" GT_REGION=210 \
      bash campaign/run_analysis.sh | head -1)
  echo "$A  (kernels ${KDIR[$i]})"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
done
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
