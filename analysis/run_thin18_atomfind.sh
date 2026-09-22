#!/usr/bin/env bash
# Experiment B: kernels + atomfind for the 70 A (18-cell) slab, once its tarball is down.
# Triage FIRST -- a leg can report COMPLETED, write an h5 and carry no NaN and still be saturated
# junk (the 2026-09-21 thick batch was 0 of 6 that way, not 1 of 6).
#
#   bash run_thin18_atomfind.sh ~/Desktop/thin18_0922/<extracted dir>
set -uo pipefail
REPO=/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/ptychoshelves-clean
PY=~/hyperspy-bundle/bin/python
R="${1:?give the extracted tarball directory}"
OUT=~/Desktop/thin18_0922/analysis; mkdir -p "$OUT/psf"
cd "$REPO"

echo "=== TRIAGE (nothing is analysed until this is read)"
$PY analysis/triage_recon.py --root "$R" || echo "  ^ legs above are unusable; analyse only the 'ok' ones"

run(){   # $1 alpha | $2 NL | $3 zdrop
    local a="$1" nl="$2" zdrop="$3"
    local meta="$R/sim_out_af_a${a}_lab_thin18/01/sim_meta.mat"
    [ -f "$meta" ] || { echo "MISSING $meta"; return; }
    for el in Pb Ti; do
        local d="$R/recon_af_a${a}_${el}_thin18_NL${nl}"
        [ -d "$d" ] || { echo "MISSING $d"; continue; }
        echo "=== extract_psf $el a${a}"
        $PY analysis/atomfind/extract_psf.py "$d" "${el}_a${a}_thin18" --zdrop "$zdrop" --out "$OUT/psf" 2>&1 | grep -v Completed | tail -4
    done
    local lab; lab=$(ls "$R/recon_af_a${a}_lab_thin18_NL${nl}"/analysis/*/*/*_recons.h5 2>/dev/null | head -1)
    [ -n "$lab" ] || { echo "MISSING lab recon for a${a}"; return; }
    [ -f "$OUT/psf/psf_Pb_a${a}_thin18_vol.npy" ] && [ -f "$OUT/psf/psf_Ti_a${a}_thin18_vol.npy" ] || { echo "NO KERNELS for a${a}"; return; }
    echo "=== atomfind a${a} (thick preset; dz and the depth constants come from sim_meta)"
    $PY analysis/atomfind/run_atomfind.py --preset thick --recon "$lab" --sim-meta "$meta" \
        --data-dir ~/Desktop/thin18_gt \
        --single-atom-vol "$OUT/psf/psf_Pb_a${a}_thin18_vol.npy" \
        --ti-kernel-vol "$OUT/psf/psf_Ti_a${a}_thin18_vol.npy" \
        --out "$OUT/atomfind_a${a}_thin18" 2>&1 | grep -vE "Completed" | tail -8
}

# zdrop = round(4 A / dz): dz = 77.935/39 = 2.00 at a70, 77.935/64 = 1.22 at a90
run 70 39 2
run 90 64 3

echo
echo "=== ladder rows (thick at step 0.02, AND its thin control at the same step)"
echo "  \$PY analysis/relaxation_ladder.py --step 5 --label 'the 70 A sample, solver step 0.02' \\"
echo "      --atomfind $OUT/atomfind_a70_thin18 $OUT/atomfind_a90_thin18 \\"
echo "      --note 'known probe, noiseless; compare against the thin leg at the SAME step'"
