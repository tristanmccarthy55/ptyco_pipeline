#!/usr/bin/env bash
# [coherence] Dry run of campaign/run_thin_atomfind.sh's partial-coherence dir names, with a stub sbatch (nothing is
# submitted; each case runs a copy of the driver in its own temporary root, so the repo gains no dirs). Checks that
#   - each coherence setting gets its own sim dir, named by its values (_fs / _dE_cc / _src / _cs, '.' -> 'p');
#   - the setting every run before 2026-10-12 used (0.3 eV, Cc 1.4 mm, 0.4 A source) keeps the old name _coh;
#   - the lab, Pb and Ti sims of a submission get the same coherence export (kernels carry the lab's physics);
#   - settings the sim would ignore are refused before anything is submitted.
#   bash campaign/test_coh_names.sh                 # the repo's driver
#   bash campaign/test_coh_names.sh <driver.sh>     # another copy, e.g. an older commit's
set -uo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DRIVER="${1:-${REPO_DIR}/campaign/run_thin_atomfind.sh}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/env/bin"; printf '#!/bin/sh\nexit 0\n' >"$TMP/env/bin/python"; chmod +x "$TMP/env/bin/python"
# The stub sbatch is an exported FUNCTION, not a file on PATH: it needs no executable file and wins over any real
# sbatch on PATH, so the test can never submit.
sbatch(){ local n; n=$(( $(cat "$STUB_N" 2>/dev/null || echo 1000) + 1 )); echo "$n" >"$STUB_N"
          printf '%s\n' "$*" >>"$STUB_LOG"; echo "$n"; }
export -f sbatch
# Every variable the driver reads (${NAME:-..} / ${NAME:+..} / ${NAME:=..}) is cleared, so nothing exported in the
# caller's shell (MODES, NL, KICK_LABELS, ...) reaches a case; BASH_ENV / ENV so no startup file runs in the driver.
DVARS=(); while IFS= read -r v; do DVARS+=("$v"); done < <(grep -oE '\$\{[A-Z][A-Z0-9_]*:[-+=]' "$DRIVER" \
    | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/' | sort -u)
FAIL=0; N=0
run_driver(){ # $1 root, then VAR=value ... -> the driver's exit status; the sbatch calls in $1/sbatch.log
    local root="$1"; shift; mkdir -p "$root/campaign"
    cp "$DRIVER" "$root/campaign/run_thin_atomfind.sh"; cp "$REPO_DIR/campaign/ceos_sweep.tsv" "$root/campaign/"
    local u=(-u BASH_ENV -u ENV) v; for v in "${DVARS[@]}"; do u+=(-u "$v"); done
    env "${u[@]}" STUB_LOG="$root/sbatch.log" STUB_N="$root/n" CONDA_ENV="$TMP/env" SHARE= \
        LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 "$@" \
        bash "$root/campaign/run_thin_atomfind.sh" >"$root/out.txt" 2>&1
}
case_ok(){ # $1 expected suffix after the mode, then VAR=value ...
    local want="$1"; shift; N=$((N + 1)); local root="$TMP/c$N" got exp
    if ! run_driver "$root" "$@"; then echo "FAIL  $* : driver refused: $(tail -1 "$root/out.txt")"; FAIL=1; return; fi
    got=$(grep -- '--job-name=af_sim_' "$root/sbatch.log" | sed -E 's/.*JOB_DIR=[^,]*\/(sim_out_af_[^,]*).*/\1/' | sort | tr '\n' ' ')
    exp=$(for m in Pb Ti lab; do echo "sim_out_af_armf_a080_${m}${want}"; done | sort | tr '\n' ' ')
    [ "$got" = "$exp" ] || { echo "FAIL  $* : sims [${got}] want [${exp}]"; FAIL=1; return; }
    # the coherence part of each sim's export: every setting must appear once per leg (lab, Pb, Ti). No paste here:
    # GNU paste -s turns empty input into one blank line, which failed every no-coherence case on Blythe (2026-10-13)
    local cohs; cohs=$(grep -- '--job-name=af_sim_' "$root/sbatch.log" \
        | grep -oE '(FOCAL_SPREAD|ENERGY_SPREAD_EV|CC_MM|SOURCE_FWHM|COHERENCE_SAMPLES)=[^, ]*' | sort | uniq -c | awk '$1 != 3')
    [ -z "$cohs" ] || { echo "FAIL  $* : lab / Pb / Ti sims got different coherence exports:"; echo "$cohs"; FAIL=1; return; }
    echo "ok    ${want:-<none>}  <-  $*"
}
case_refused(){ # VAR=value ... that must be refused before any sbatch
    N=$((N + 1)); local root="$TMP/c$N"
    if run_driver "$root" "$@"; then echo "FAIL  $* : accepted (dirs $(grep -c af_sim_ "$root/sbatch.log" 2>/dev/null))"; FAIL=1; return; fi
    [ -s "$root/sbatch.log" ] && { echo "FAIL  $* : refused only after submitting"; FAIL=1; return; }
    echo "ok    refused  <-  $*   ($(tail -1 "$root/out.txt"))"
}
echo "driver: ${DRIVER} (${#DVARS[@]} driver variables cleared per case)"
# the old cold-FEG setting keeps its name, however the numbers are written (per-species phonon sims carry v2 since the 2026-10-30 abTEM sigma fix)
case_ok _coh            ENERGY_SPREAD_EV=0.3 CC_MM=1.4 SOURCE_FWHM=0.4
case_ok _coh            ENERGY_SPREAD_EV=0.30 CC_MM=1.40 SOURCE_FWHM=0.40
case_ok _ph16v2_coh     PHONONS=16 PER_SPECIES=1 ENERGY_SPREAD_EV=0.3 CC_MM=1.4 SOURCE_FWHM=0.4
case_ok _ph16v2         PHONONS=16 PER_SPECIES=1
case_ok _ph16v2_T100    PHONONS=16 PER_SPECIES=1 PHONON_T=100
case_ok ""
# new: one dir per value
case_ok _fs7p3          FOCAL_SPREAD=7.3
case_ok _fs7p3          FOCAL_SPREAD=7.30
case_ok _fs3p5          FOCAL_SPREAD=3.5
case_ok _fs0p75         FOCAL_SPREAD=0.75
case_ok _src0p2         SOURCE_FWHM=0.2
case_ok _src0p8         SOURCE_FWHM=0.8
case_ok _fs2p5_src0p4   FOCAL_SPREAD=2.5 SOURCE_FWHM=0.4
case_ok _dE0p3_cc1p4    ENERGY_SPREAD_EV=0.3 CC_MM=1.4
case_ok _dE0p09_cc2p4   ENERGY_SPREAD_EV=0.09 CC_MM=2.4
case_ok _dE0p3_cc1p4_src0p2 ENERGY_SPREAD_EV=0.3 CC_MM=1.4 SOURCE_FWHM=0.2
case_ok _fs1p5_cs32     FOCAL_SPREAD=1.5 COHERENCE_SAMPLES=32
case_ok _ph16v2_T100_fs2p5 PHONONS=16 PER_SPECIES=1 PHONON_T=100 FOCAL_SPREAD=2.5
# refused: the sim would ignore part of it, or the value is not a spread
case_refused ENERGY_SPREAD_EV=0.3
case_refused CC_MM=1.4
case_refused FOCAL_SPREAD=3 ENERGY_SPREAD_EV=0.3 CC_MM=1.4
case_refused PHONONS=16 PER_SPECIES=1 FOCAL_SPREAD=2 COHERENCE_SAMPLES=32
case_refused COHERENCE_SAMPLES=32
case_refused FOCAL_SPREAD=0
case_refused FOCAL_SPREAD=abc
case_refused SOURCE_FWHM=-0.4
[ "$FAIL" = 0 ] && echo "all ${N} cases passed" || { echo "SOME CASES FAILED"; exit 1; }
