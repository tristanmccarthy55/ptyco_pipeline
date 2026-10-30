#!/usr/bin/env bash
# The CEILING runs for the Spectra 300 Mono (the user's plan, 2026-10-16): ARM tableau, round Gaussian focal spread and
# source, noiseless, probe update with modes from the mixed-state start (true central probe), engine with the all-modes
# illumination fix. Two parts, one submission:
#  A. coherence only, at the realistic Spectra Mono settings: focal spread 2.5 A (~90 meV through Cc ~1.6 mm) with a
#     0.1 and a 0.2 A source -- mode legs on the kept joint sims (RECON_ONLY; their fixed controls are step 2's _fixed).
#  B. + phonons, liquid nitrogen (100 K) and liquid helium (20 K: the Debye model's zero-point floor, ~4 K values), at
#     fs 2.5 + src 0.1 -- new sims (16 configurations sample phonons and coherence jointly), each with a fixed-probe
#     control (FIXED_TWIN -> _mixF_fixed) and the mode legs (_mixF). Extend to src 0.2 if A shows it survives.
# 8 modes (the joint states need more than 8 for 99 %: capped, as before), 250 it, "4;2", 96G, 48 h.
#   bash aberration_experiment/mono300/ceiling_submit.sh
# Submits nothing if the pulled engine lacks the fix, if the shell exports a job variable, if a dry run fails, if A's sims
# are not all there, or if B's sims already exist. Output, job ids included: logs/ceiling_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
# RETIRED 2026-10-30. This ran once (2026-10-16): A's mode legs, and B's phonon sims, which moved every atom sqrt(3) too far
# (abTEM 1.0.9; sim/simulate_4dstem.py a8b9264: "100 K" was Pb ~310 K, "20 K" Pb ~110 K / Ti ~295 K). Run again it would
# repeat A and rebuild B under the driver's new _ph16v2 names. The true-temperature rerun of B gets its own script.
echo "NOT SUBMITTED: ceiling_submit.sh is retired (2026-10-30, the phonon fix): see aberration_experiment/mono300/PHASE1.md"; exit 1
mkdir -p logs; LOG="logs/ceiling_$(date +%Y%m%d_%H%M%S).txt"
A_PTS=("_fs2p5_src0p1" "_fs2p5_src0p2")
B_PTS=("100 _ph16_T100_fs2p5_src0p1" "20 _ph16_T20_fs2p5_src0p1")
MODEOPTS=(KICK_ALL=1 KICK_START=mixed KICK_MODES=8 KICK_NITER=250 KICK_GROUPING="4;2" RMEM=96G RTAG=_mixF)
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
for S in "${A_PTS[@]}"; do for m in lab Pb Ti; do
  [ -e "sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5" ] || { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${S} has no data (A reuses step 2's joint sims)"; exit 1; }
done; done
for X in "${B_PTS[@]}"; do read -r K S <<<"$X"; for m in lab Pb Ti; do
  [ -e "sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5" ] && { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${S} already has data"; exit 1; }
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
echo "== A. coherence only, mode legs at the realistic Spectra Mono settings (kept joint sims)"
submit _fs2p5_src0p1 _mixF -- FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 RECON_ONLY=1 "${MODEOPTS[@]}" || exit 1
submit _fs2p5_src0p2 _mixF -- FOCAL_SPREAD=2.5 SOURCE_FWHM=0.2 RECON_ONLY=1 "${MODEOPTS[@]}" || exit 1
echo "== B. + phonons at 100 K and 20 K, fs 2.5 + src 0.1 (new sims; fixed control + mode legs)"
for X in "${B_PTS[@]}"; do
  read -r K S <<<"$X"
  submit "$S" _mixF,_mixF_fixed -- PHONONS=16 PER_SPECIES=1 PHONON_T="$K" FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 FIXED_TWIN=1 \
         "${MODEOPTS[@]}" || exit 1
done
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
