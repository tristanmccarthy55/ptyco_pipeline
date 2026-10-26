#!/usr/bin/env bash
# detector256, the LIMITS wave (2026-10-26; README "Limits wave"). R2 (256 px at +-96 mrad, modelled by upsampling) kept
# the Pb depths but lost Ti. Each run below changes ONE thing against R2 (256 px, +-96 mrad, 0.5 A step, 200 + 200 it),
# true probe fixed, fs 2.5 A + src 0.1 A noiseless, lab + Pb + Ti alike, every analysis with 2.22 A kernel crops chained on
# its pack:
#   R2_n400   R2 at 400 + 400 iterations     reuses R2's copies             was R2 short of convergence?
#   R2_n1000  R2 at 1000 + 1000 iterations   reuses R2's copies
#   C512      512 px  (2 x 2 sums, upsampling 1)     new copies of the kept sims    how many pixels does Ti need?
#   C128      128 px  (8 x 8 sums, upsampling 3)     new copies                     an EMPAD-class camera: does Pb survive?
#   R4        354 px at +-133 mrad (4 x 4, up 2)     new copies                     does more angle help at this pixel size?
#   S04       R2 on a 0.4 A scan (2500 positions)    NEW sims, label armf_a080_s04  does a denser scan bring Ti back?
# Run from anywhere on Blythe:   bash aberration_experiment/detector256/submit_limits.sh
# Submits nothing if the shell exports a variable the jobs read, if a dry run or the frame unit test fails, if a kept sim or
# R2's copies are missing, or if the dense-scan sims are half there. Log (every job id): logs/detector256_limits_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
mkdir -p logs; TS0="$(date +%Y%m%d_%H%M%S)"; LOG="logs/detector256_limits_${TS0}.txt"
SFX=_fs2p5_src0p1
{
LEAK=""
for v in $( { grep -ohE '\$\{[A-Z][A-Z0-9_]*:[-+=?]' campaign/run_thin_atomfind.sh campaign/run_analysis.sh run_recon_synthetic_ML.slurm \
                | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/'
              LC_ALL=C grep -a -ohE "getenv\('[A-Z0-9_]+'\)" ptycho/run_synthetic_recon_ML.m | sed -E "s/getenv\('([A-Z0-9_]+)'\)/\1/"; } | sort -u); do
  case "$v" in SHARE|CONDA_ENV|TMPDIR|SLURM_*|CUDA_*|CUPY_CACHE_DIR|MATLAB_*|MPLCONFIGDIR|NUMBA_CACHE_DIR|XDG_CACHE_HOME) continue;; esac
  printenv "$v" >/dev/null && LEAK="${LEAK} ${v}=$(printenv "$v")"   # EXPORTED only: that is what reaches the jobs
done
[ -z "$LEAK" ] || { echo "NOT SUBMITTED: your shell exports variables the jobs read:${LEAK} -- unset them or open a new shell"; exit 1; }
PY="${CONDA_ENV:-${SHARE:-}/phucrh/envs/abtem}/bin/python"; [ -x "$PY" ] || { echo "NOT SUBMITTED: no python at ${PY}"; exit 1; }
T=$(bash campaign/test_coh_names.sh); echo "$T" | grep -v '^ok'
[ "$(echo "$T" | tail -1)" = "all 26 cases passed" ] || { echo "NOT SUBMITTED: the naming dry run failed (above)"; exit 1; }
T=$(PYTHON="$PY" bash campaign/test_driver_legacy.sh HEAD); echo "$T" | grep -v '^same'
[ "$(echo "$T" | tail -1)" = "all launch lines unchanged" ] || { echo "NOT SUBMITTED: the old-launch dry run failed (above)"; exit 1; }
T=$(bash aberration_experiment/detector256/test_det_driver.sh); echo "$T" | grep -v '^ok'
[ "$(echo "$T" | tail -1)" = "all 25 cases passed" ] || { echo "NOT SUBMITTED: the detector dry run failed (above)"; exit 1; }
T=$("$PY" aberration_experiment/detector256/test_upsampling_frame.py 2>&1); echo "$T" | grep -E '^FAIL|PASSED|FAILED'
[ "$(echo "$T" | tail -1)" = "ALL PASSED" ] || { echo "NOT SUBMITTED: the frame unit test failed (above)"; exit 1; }
for m in lab Pb Ti; do                        # the kept 0.5 A sims, and R2's copies (reused by R2_n400 / R2_n1000)
  for f in data_dp.hdf5 data_position.hdf5 sim_meta.mat probe_initial_true.mat; do
    [ -e "sim_out_af_armf_a080_${m}${SFX}/01/${f}" ] || { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${SFX}/01/${f} missing"; exit 1; }
  done
  [ -e "sim_out_af_armf_a080_${m}${SFX}_detc96b4/01/rebin_detector.json" ] \
      || { echo "NOT SUBMITTED: R2's copy sim_out_af_armf_a080_${m}${SFX}_detc96b4 is not whole"; exit 1; }
done
n=0; for m in lab Pb Ti; do [ -e "sim_out_af_armf_a080_s04_${m}${SFX}/01/data_dp.hdf5" ] && n=$((n + 1)); done
case $n in 0) S04_RO=();; 3) S04_RO=(RECON_ONLY=1); echo "(S04: its sims are there -- reused)";;
     *) echo "NOT SUBMITTED: ${n} of 3 dense-scan sims present (sim_out_af_armf_a080_s04_*${SFX}) -- restore or delete the rest"; exit 1;; esac
AIDS=""; PULL=()
submit(){ # $1 run name, $2 label, $3 recon tag after the suffix, then VAR=value ... for the driver
  local name="$1" lab="$2" rt="$3"; shift 3
  echo "== ${name}"
  OUT=$(env "$@" LABELS="$lab" TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at ${name}"; return 1; }
  local PJ; PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  echo "$OUT" | grep -q "recon dirs ${rt%%_n[0-9]*}$" && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: ${name}: recon dirs / pack"; return 1; }
  echo "$OUT" | sed -n '2,3p;/^pack /p'
  A=$(LABELS="$lab" SUFFIX="${SFX}${rt}" DEP="$PJ" GT_REGION=210 PSF_HALF_A=2.22 bash campaign/run_analysis.sh | head -1)
  echo "$A"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
  PULL+=("atomfind_results_${lab}${SFX}${rt}_2026*.tgz" "analysis_${lab}${SFX}${rt}_kh2p22_2026*.tgz")
}
submit R2_n400  armf_a080 _detc96b4_up2_n400  RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2 NITER=400 RTAG=_n400 || exit 1
submit R2_n1000 armf_a080 _detc96b4_up2_n1000 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2 NITER=1000 RTAG=_n1000 || exit 1
submit C512     armf_a080 _detc96b2_up1       RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=2 UPSAMPLING=1 || exit 1
submit C128     armf_a080 _detc96b8_up3       RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=8 UPSAMPLING=3 || exit 1
submit R4       armf_a080 _detc133b4_up2      RECON_ONLY=1 DET_CROP_MRAD=133 DET_BIN=4 UPSAMPLING=2 || exit 1
submit S04      armf_a080_s04 _detc96b4_up2   ${S04_RO[@]+"${S04_RO[@]}"} DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2 || exit 1
echo "ANALYSES:${AIDS}"
echo "done when: sacct -j $(echo ${AIDS} | tr ' ' ',') shows all COMPLETED (S04 last: sims ~7 h, then ~4 h recons). Then, on the Mac:"
echo "  mkdir -p ~/Desktop/detector256_limits && cd ~/Desktop/detector256_limits"
for g in "${PULL[@]}"; do echo "  scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/${g}' ."; done
echo "  for t in *.tgz; do tar xzf \"\$t\"; done"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
