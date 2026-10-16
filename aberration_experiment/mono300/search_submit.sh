#!/usr/bin/env bash
# C. THE REALISTIC START (the user's plan, 2026-10-16): the same physics as the ceiling A (ARM tableau, round Gaussian focal
# spread 2.5 A + source 0.1 A, noiseless, no phonons) but the probe is NOT given: C1, C3 and A1 errors are DRAWN from the
# corrector's quoted precision (ARM run 3, 95 %: C1 39.5 A, C3 0.804 um, A1 59 A; seed 1), an outer search bounded by
# that precision finds them (stage 1 C1 x C3 7 x 7, stage 2 the A1 vector 5 x 5: 74 fixed-probe trials of 50 it), and every
# leg (lab, Pb, Ti: one probe per instrument, the lab leg's answer) refines from the searched probe with the mixed state built
# on it -- probe update with 8 modes, 250 it, "4;2". Higher orders stay at their true values (stated). On the kept joint sims
# (RECON_ONLY). Dirs recon_af_armf_a080_<leg>_fs2p5_src0p1_ps_NL18, search_af_armf_a080_lab_fs2p5_src0p1/; analysis
# SUFFIX _fs2p5_src0p1_ps. Extend to src 0.2 if the ceiling A shows it survives.
#   bash aberration_experiment/mono300/search_submit.sh
# Submits nothing if the pulled engine lacks the fix, if the shell exports a job variable, if a dry run fails, or if the
# sims are not all there. Output, job ids included: logs/search_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; LOG="logs/search_$(date +%Y%m%d_%H%M%S).txt"
PTS=("_fs2p5_src0p1")
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
for S in "${PTS[@]}"; do for m in lab Pb Ti; do
  [ -e "sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5" ] || { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${S} has no data (C reuses step 2's joint sims)"; exit 1; }
done; done
AIDS=""
submit(){ # $1 expected dir suffix, then the analysis suffixes (comma list), then -- and VAR=value ... for the driver
  local want="$1" ats="$2"; shift 3
  OUT=$(env "$@" LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at ${want}"; return 1; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ "$SFX" = "$want" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: suffix '${SFX}' (want ${want}) / pack '${PJ}'"; return 1; }
  echo "$OUT" | sed -n '2p;/^pack /p'
  local at
  for at in ${ats//,/ }; do
    A=$(LABELS=armf_a080 SUFFIX="${SFX}${at}" DEP="$PJ" GT_REGION=210 bash campaign/run_analysis.sh | head -1)
    echo "$A"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
  done
}
echo "== C. drawn start (seed 1) -> staged outer search (C1 x C3, A1) -> refine with the mixed state on the searched probe"
submit _fs2p5_src0p1 _ps -- FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 RECON_ONLY=1 PSEARCH=1 PSEARCH_SEED=1 \
       KICK_ALL=1 KICK_START=mixed KICK_MODES=8 KICK_NITER=250 KICK_GROUPING="4;2" RMEM=96G || exit 1
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
