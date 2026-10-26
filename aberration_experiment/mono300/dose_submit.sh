#!/usr/bin/env bash
# D. THE DOSE SERIES, 1e8 FIRST (the user, 2026-10-26: 1e8 as the gate; 1e5-1e7 only if it holds). The Spectra 300 Mono
# point with everything put back at once: ARM tableau, focal spread 2.5 A + source 0.1 A, phonons at 20 K (B's sims, kept),
# shot noise, and the probe NOT given -- C's seed-1 draw from the corrector's precision, the 3-stage search afresh on the
# noisy data (C1 x C3 7 x 7, the A1 vector 5 x 5, C1 x C3 again: 123 fixed-probe trials of 50 it), then every leg (lab, Pb,
# Ti; Pb and Ti take the lab leg's answer) refined from the search's answer with the mixed state built on it -- 8 modes,
# 250 it, "4;2". Noise: sim/add_poisson_noise.py on B's noiseless sims, a CPU job per leg and dose (seeds 0/1/2 for
# lab/Pb/Ti). Dirs sim_out_af_armf_a080_<leg>_ph16_T20_fs2p5_src0p1_dose<D>, search_af_armf_a080_lab_ph16_T20_fs2p5_src0p1_
# dose<D>/, recon_af_armf_a080_<leg>_ph16_T20_fs2p5_src0p1_dose<D>_ps_NL18; analysis SUFFIX _ph16_T20_fs2p5_src0p1_dose<D>_ps.
# ~230 GPU-h a dose (123 trials x ~1.3 h + 3 legs x 23.5 h, the measured C and B legs), ~1.5-2 days with the queue.
#   bash aberration_experiment/mono300/dose_submit.sh 1e8               # the gate
#   bash aberration_experiment/mono300/dose_submit.sh "1e7 1e6 1e5"     # the rest, if 1e8 holds
# Submits nothing if the pulled engine lacks the fix, if the shell exports a job variable, if a dry run fails, if B's 20 K
# sims are not all there, if a dose is not one of 1e5-1e8, or if a dose was searched before. Output: logs/dose_<ts>.txt.
# Its own dry run (stub sbatch, throwaway git copy): bash aberration_experiment/mono300/test_dose_submit.sh
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
DOSE_LIST="${1:-}"
mkdir -p logs; LOG="logs/dose_$(date +%Y%m%d_%H%M%S).txt"
S="_ph16_T20_fs2p5_src0p1"                                  # B's 20 K sims at the Mono point
{
[ -n "$DOSE_LIST" ] || { echo "NOT SUBMITTED: give the doses, e.g. dose_submit.sh 1e8"; exit 1; }
ATS=""
for d in $DOSE_LIST; do
  case "$d" in 1e5|1e6|1e7|1e8) ;; *) echo "NOT SUBMITTED: dose '${d}': one of 1e5 1e6 1e7 1e8 (written so)"; exit 1;; esac
  [ -e "search_af_armf_a080_lab${S}_dose${d}" ] && { echo "NOT SUBMITTED: search_af_armf_a080_lab${S}_dose${d} exists (D at ${d} was submitted before)"; exit 1; }
  ATS="${ATS:+${ATS},}_dose${d}_ps"
done
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
for m in lab Pb Ti; do
  [ -e "sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5" ] || { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${S} has no data (D noises B's 20 K sims)"; exit 1; }
done
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
echo "== D at ${DOSE_LIST} e/A^2: noise on B's 20 K sims -> seed-1 draw -> 3-stage search -> refine with the mixed state"
submit "$S" "$ATS" -- DOSES="$DOSE_LIST" PHONONS=16 PER_SPECIES=1 PHONON_T=20 FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 RECON_ONLY=1 \
       PSEARCH=1 PSEARCH_SEED=1 PSEARCH_STAGES=3 KICK_ALL=1 KICK_START=mixed KICK_MODES=8 KICK_NITER=250 KICK_GROUPING="4;2" \
       RMEM=96G || exit 1
echo "ANALYSES:${AIDS}"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
