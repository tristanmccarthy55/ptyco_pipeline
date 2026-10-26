#!/usr/bin/env bash
# detector256: is a 256 x 256 camera enough for the 80 mrad runs? (aberration_experiment/detector256/README.md)
# Every leg on the KEPT joint sims at focal spread 2.5 A + source 0.1 A (ARM tableau, 80 mrad, noiseless): true probe fixed,
# the reference control's settings (GROUPING 32;16, NITER 200 + 200, BETA_LSQ 0.05, NL 18), lab + Pb + Ti alike (the
# kernels see the same camera); each run's analysis chained on its pack. Run from anywhere on Blythe:
#   bash aberration_experiment/detector256/submit.sh        # R1 (+-96 mrad, full pixels) + R2 (+-96 mrad, 4x4 -> 256 px, upsampling 2)
#   bash aberration_experiment/detector256/submit.sh r3     # R3 (R2's 256 px, no upsampling): once R2's copies are whole
#   bash aberration_experiment/detector256/submit.sh r4     # R4 (+-133 mrad, 4x4 -> 354 px, upsampling 2): only if R2 works
# Submits nothing if the shell exports a variable the jobs read, if a dry run fails, if the frame unit test fails (it also
# checks this checkout carries the 2026-10-23 engine changes), or if a sim is missing. Log: logs/detector256_<what>_<ts>.txt.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
WHAT="${1:-r1r2}"
case "$WHAT" in r1r2|r3|r4) ;; *) echo "usage: $0 [r1r2|r3|r4]"; exit 1;; esac
mkdir -p logs; LOG="logs/detector256_${WHAT}_$(date +%Y%m%d_%H%M%S).txt"
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
for m in lab Pb Ti; do                        # the kept sims (disk cleanup 2026-10-14 kept the joint ones), whole
  for f in data_dp.hdf5 data_position.hdf5 sim_meta.mat probe_initial_true.mat; do
    [ -e "sim_out_af_armf_a080_${m}${SFX}/01/${f}" ] || { echo "NOT SUBMITTED: sim_out_af_armf_a080_${m}${SFX}/01/${f} missing"; exit 1; }
  done
done
"$PY" - "sim_out_af_armf_a080_lab${SFX}/01" <<'EOF' || { echo "NOT SUBMITTED: the lab sim's geometry is not the reference's"; exit 1; }
import sys, h5py
from scipy.io import loadmat
d = sys.argv[1]; m = loadmat(d + "/sim_meta.mat", simplify_cells=True)["meta"]
with h5py.File(d + "/data_dp.hdf5", "r") as f: sh = f["dp"].shape
print(f"lab sim: data {sh}, Ndpx {m['Ndpx']}, {m['d_alpha_rad'] * 1e3:.5f} mrad/px (+-{m['Ndpx'] // 2 * m['d_alpha_rad'] * 1e3:.1f}), "
      f"object pixel {m['dx_object_A']:.5f} A, fs {m['focal_spread_A']} A, src {m['source_size_fwhm_A']} A")
sys.exit(0 if (m["Ndpx"] == 1420 and sh[1:] == (1420, 1420) and m["focal_spread_A"] == 2.5 and m["source_size_fwhm_A"] == 0.1) else 1)
EOF
AIDS=""; PULL=()
submit(){ # $1 recon tag after the suffix (e.g. _detc96b4_up2), then VAR=value ... for the driver
  local rt="$1"; shift
  OUT=$(env "$@" RECON_ONLY=1 LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 \
        FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at $*"; return 1; }
  local PJ; PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  echo "$OUT" | grep -q "recon dirs ${rt}$" && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: recon dirs not ${rt} / no pack for $*"; return 1; }
  echo "$OUT" | sed -n '2,3p;/^pack /p'
  # kernels cropped to +-2.22 A (the reference's 30 px): at the 0.1025 A pixel 30 px reaches the grid's neighbours (README)
  A=$(LABELS=armf_a080 SUFFIX="${SFX}${rt}" DEP="$PJ" GT_REGION=210 PSF_HALF_A=2.22 bash campaign/run_analysis.sh | head -1)
  echo "$A"; AIDS="${AIDS} $(echo "$A" | awk '{print $2}')"
  PULL+=("atomfind_results_armf_a080${SFX}${rt}_2026*.tgz" "analysis_armf_a080${SFX}${rt}_kh2p22_2026*.tgz")
}
case "$WHAT" in
  r1r2)
    echo "== R1: +-96 mrad, full pixels (1420 -> 1024 px; window 105 A, object pixel 0.1025 A)"
    submit _detc96 DET_CROP_MRAD=96 || exit 1
    echo "== R2: +-96 mrad, 4x4 sums -> 256 px, the solver models the sums (upsampling 2: 1024 px model, 105 A window)"
    submit _detc96b4_up2 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2 || exit 1;;
  r3)
    for m in lab Pb Ti; do [ -e "sim_out_af_armf_a080_${m}${SFX}_detc96b4/01/rebin_detector.json" ] \
        || { echo "NOT SUBMITTED: R2's copy sim_out_af_armf_a080_${m}${SFX}_detc96b4 is not whole yet (its rebin job: sacct --name=af_rebin_armf_a080_${m}${SFX}_detc96b4)"; exit 1; }; done
    echo "== R3: R2's 256 px data, no upsampling (the solver's window 26 A: 56 % of the probe outside it)"
    submit _detc96b4 DET_CROP_MRAD=96 DET_BIN=4 || exit 1;;
  r4)
    echo "== R4: +-133 mrad, 4x4 sums -> 354 px (1420 -> 1416 -> 354), upsampling 2 (1416 px model, 105 A window)"
    submit _detc133b4_up2 DET_CROP_MRAD=133 DET_BIN=4 UPSAMPLING=2 || exit 1;;
esac
echo "ANALYSES:${AIDS}"
echo "done when: sacct -j $(echo ${AIDS} | tr ' ' ',') shows COMPLETED. Then, on the Mac, into a fresh folder:"
echo "  mkdir -p ~/Desktop/detector256_${WHAT} && cd ~/Desktop/detector256_${WHAT}"
for g in "${PULL[@]}"; do echo "  scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/${g}' ."; done
echo "  for t in *.tgz; do tar xzf \"\$t\"; done"
} 2>&1 | tee "$LOG"
echo "log: $(pwd)/${LOG}"
