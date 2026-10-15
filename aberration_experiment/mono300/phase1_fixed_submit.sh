#!/usr/bin/env bash
# Phase 1, step 1 (the user's choice, 2026-10-13): at the ten phase-1 points, the lab, Pb and Ti sims plus their
# true-probe-fixed legs (_<point>_fixed, 1 mode, 200 it, "32;16"), each point's analysis chained on its pack. The sims are
# KEPT (CLEANDATA=0) for step 2, the probe-mode legs (RECON_ONLY on these sims). Run from anywhere on Blythe:
#   bash aberration_experiment/mono300/phase1_fixed_submit.sh
# Submits nothing if the shell exports a variable the driver or the analysis reads (it would change every leg: MODES=lab
# drops the kernels) or if the naming dry run fails. Output, job ids included: logs/phase1_fixed_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; LOG="logs/phase1_fixed_$(date +%Y%m%d_%H%M%S).txt"
{
LEAK=""
for v in $(grep -ohE '\$\{[A-Z][A-Z0-9_]*:[-+=?]' campaign/run_thin_atomfind.sh campaign/run_analysis.sh \
           | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/' | sort -u); do
  case "$v" in SHARE|CONDA_ENV) continue;; esac              # where the share and the python env are: meant to be set
  printenv "$v" >/dev/null && LEAK="${LEAK} ${v}=$(printenv "$v")"   # EXPORTED only: that is what reaches the jobs
done
[ -z "$LEAK" ] || { echo "NOT SUBMITTED: your shell exports driver variables:${LEAK} -- unset them or open a new shell"; exit 1; }
T=$(bash campaign/test_coh_names.sh); echo "$T" | grep -v '^ok'
[ "$(echo "$T" | tail -1)" = "all 26 cases passed" ] || { echo "NOT SUBMITTED: the naming dry run failed (above)"; exit 1; }
for P in FOCAL_SPREAD=7.29 FOCAL_SPREAD=5 FOCAL_SPREAD=3.5 FOCAL_SPREAD=2.5 FOCAL_SPREAD=1.5 FOCAL_SPREAD=0.75 \
         SOURCE_FWHM=0.05 SOURCE_FWHM=0.1 SOURCE_FWHM=0.2 SOURCE_FWHM=0.4; do
  OUT=$(env "$P" RTAG=_fixed LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at $P"; break; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ -n "$SFX" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: no suffix / pack id for $P"; break; }
  echo "$OUT" | sed -n '2p;/^pack /p'
  LABELS=armf_a080 SUFFIX="${SFX}_fixed" DEP="$PJ" GT_REGION=210 bash campaign/run_analysis.sh | head -1
done
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
