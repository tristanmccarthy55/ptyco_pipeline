#!/usr/bin/env bash
# Kernels + atomfind for the 2026-09-22 six-fold THRESHOLD ladder (0.1 / 0.2 / 0.3 / 0.45 waves).
# Per rung: extract the matched Pb and Ti kernels from that rung's OWN grid recons, then run the
# blind finder on that rung's lab recon with those kernels. Same pattern as run_0921_atomfind.sh.
set -uo pipefail
REPO=/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/ptychoshelves-clean
PY=~/hyperspy-bundle/bin/python
R=~/Desktop/relax_0922
NR=$R/atomfind_results_nr0p1_C56_0p1w-nr0p2_C56_0p2w-nr0p3_C56_0p3w-nr0p45_C56_0p45w_20260921_225614
OUT=$R/analysis_0922; mkdir -p $OUT/psf
cd "$REPO"

run(){   # $1 tag | $2 root | $3 dir prefix before the element | $4 NL | $5 zdrop
    local tag="$1" root="$2" pre="$3" nl="$4" zdrop="$5"
    local dz; dz=$($PY -c "print(27.525/$nl)")
    for el in Pb Ti; do
        local d="$root/${pre}_${el}_NL${nl}"
        [ -d "$d" ] || { echo "MISSING $d"; continue; }
        echo "=== extract_psf $el $tag"
        $PY analysis/atomfind/extract_psf.py "$d" "${el}_${tag}" --zdrop "$zdrop" --out "$OUT/psf" 2>&1 | grep -v Completed | tail -4
    done
    local lab; lab=$(ls "$root/${pre}_lab_NL${nl}"/analysis/*/*/*_recons.h5 2>/dev/null | head -1)
    [ -n "$lab" ] || { echo "MISSING lab recon for $tag"; return; }
    [ -f "$OUT/psf/psf_Pb_${tag}_vol.npy" ] && [ -f "$OUT/psf/psf_Ti_${tag}_vol.npy" ] || { echo "NO KERNELS for $tag"; return; }
    echo "=== atomfind $tag  (dz=$dz)"
    $PY analysis/atomfind/run_atomfind.py --preset thin --recon "$lab" --dz "$dz" \
        --data-dir ~/Desktop/thin_ab_af/gtdata \
        --single-atom-vol "$OUT/psf/psf_Pb_${tag}_vol.npy" --ti-kernel-vol "$OUT/psf/psf_Ti_${tag}_vol.npy" \
        --out "$OUT/atomfind_${tag}" 2>&1 | grep -v Completed | tail -6
}

for lbl in nr0p1_C56_0p1w nr0p2_C56_0p2w nr0p3_C56_0p3w nr0p45_C56_0p45w; do
    run "$lbl" "$NR" "recon_af_${lbl}" 14 2
done
echo "ALL DONE"
