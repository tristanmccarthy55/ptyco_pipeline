#!/usr/bin/env bash
# Phase 1, step 2 (the user's choice, 2026-10-14), in one go:
#  A. JOINT points, true probe fixed, new sims (kept): focal spread 2.5 and 3.5 A, each with a 0.1 and a 0.2 A source.
#  B. PROBE-MODE legs where the step-1 control failed or weakened -- focal spread 5 A (7 modes) and 7.29 A (8), source 0.2
#     and 0.4 A (8) -- each started from its sim's own mixed state (KICK_START=mixed: sim/make_mixed_probe.py), probe
#     released (PSTART 40 / 20), 250 it, "4;2", 96G; lab, Pb and Ti alike, so the kernels come out of the same operator.
#     Run on the step-1 sims (RECON_ONLY); a point whose sims were deleted is re-simulated in the same submission.
# A goes in first, so its sims and recons queue ahead of the ~30-34 h mode legs. Each submission's analysis is chained on
# its pack (SUFFIX _<point>_fixed / _<point>_mix). Run from anywhere on Blythe:
#   bash aberration_experiment/mono300/phase1_step2_submit.sh
# Submits nothing if the shell exports a variable the driver, the analysis, the recon job or its MATLAB script reads, if a
# dry run fails, or if a mode leg's sims are half there. Output, job ids included: logs/phase1_step2_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; LOG="logs/phase1_step2_$(date +%Y%m%d_%H%M%S).txt"
JOINT=("FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 _fs2p5_src0p1" "FOCAL_SPREAD=2.5 SOURCE_FWHM=0.2 _fs2p5_src0p2"
       "FOCAL_SPREAD=3.5 SOURCE_FWHM=0.1 _fs3p5_src0p1" "FOCAL_SPREAD=3.5 SOURCE_FWHM=0.2 _fs3p5_src0p2")
MIX=("FOCAL_SPREAD=5 _fs5 7" "FOCAL_SPREAD=7.29 _fs7p29 8" "SOURCE_FWHM=0.2 _src0p2 8" "SOURCE_FWHM=0.4 _src0p4 8")
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
for X in "${MIX[@]}"; do                     # a mode leg's sims must be whole (reuse) or gone (re-simulate), never half
  read -r P S M <<<"$X"; n=0
  for m in lab Pb Ti; do [ -e "sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5" ] && n=$((n + 1)); done
  case $n in 0|3) ;; *) echo "NOT SUBMITTED: ${S}: ${n} of 3 legs' sim data present -- restore or delete the rest"; exit 1;; esac
done
AIDS=""
submit(){ # $1 expected dir suffix, $2 analysis tag (_fixed / _mix), then VAR=value ... for the driver
  local want="$1" at="$2"; shift 2
  OUT=$(env "$@" LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at $*"; return 1; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ "$SFX" = "$want" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: suffix '${SFX}' (want ${want}) / pack '${PJ}' for $*"; return 1; }
  echo "$OUT" | sed -n '2p;/^pack /p'
  A=$(LABELS=armf_a080 SUFFIX="${SFX}${at}" DEP="$PJ" GT_REGION=210 bash campaign/run_analysis.sh | head -1)
  echo "$A"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
}
echo "== A. joint points, true probe fixed (new sims, kept)"
for X in "${JOINT[@]}"; do
  read -r P1 P2 S <<<"$X"; submit "$S" _fixed "$P1" "$P2" RTAG=_fixed || exit 1
done
echo "== B. probe-mode legs from each sim's own mixed state"
for X in "${MIX[@]}"; do
  read -r P S M <<<"$X"; RO=()
  [ -e "sim_out_af_armf_a080_lab${S}/01/data_dp.hdf5" ] && RO=(RECON_ONLY=1) || echo "  (${S}: sims gone -- re-simulating them first)"
  submit "$S" _mix "$P" ${RO[@]+"${RO[@]}"} RTAG=_mix KICK_ALL=1 KICK_START=mixed KICK_MODES="$M" KICK_NITER=250 \
         KICK_GROUPING="4;2" RMEM=96G || exit 1
done
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
