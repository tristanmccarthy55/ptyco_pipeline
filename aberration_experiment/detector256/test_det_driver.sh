#!/usr/bin/env bash
# [detector256] Dry run of campaign/run_thin_atomfind.sh's DETECTOR switches (DET_CROP_MRAD / DET_BIN / UPSAMPLING), with a
# stub sbatch (an exported function: nothing is ever submitted; each case runs a copy of the driver in its own temporary
# root). Checks, for R1-R4 of aberration_experiment/detector256/README.md:
#   - one rebin job per leg (lab, Pb, Ti) with the right sim/rebin_detector.py arguments, the legs chained on it (afterok);
#     with RECON_ONLY and the copy already there, no rebin job and the legs run at once;
#   - the recon dirs' names, their data links (the copy) and probe link (probe_initial_true_up<k>.mat of the copy);
#   - UPSAMPLING reaches the recon job's export only when the solver should model the sums; BETA_LSQ / REGLAYER / NITER /
#     GROUPING exactly as the reference control's (fs 2.5 + src 0.1, true probe fixed);
#   - settings that make no sense, or that the driver does not support with a detector yet, are refused before any sbatch.
#   bash aberration_experiment/detector256/test_det_driver.sh
set -uo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DRIVER="${1:-${REPO_DIR}/campaign/run_thin_atomfind.sh}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/env/bin"; printf '#!/bin/sh\nexit 0\n' >"$TMP/env/bin/python"; chmod +x "$TMP/env/bin/python"
sbatch(){ local n; n=$(( $(cat "$STUB_N" 2>/dev/null || echo 1000) + 1 )); echo "$n" >"$STUB_N"
          printf '%s\n' "$*" >>"$STUB_LOG"; echo "$n"; }
export -f sbatch
DVARS=(); while IFS= read -r v; do DVARS+=("$v"); done < <(grep -oE '\$\{[A-Z][A-Z0-9_]*:[-+=]' "$DRIVER" \
    | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/' | sort -u)
FAIL=0; N=0
SFX=_fs2p5_src0p1
run_driver(){ # $1 root, $2 "sims" | "sims+copy<tag>" | "none", then VAR=value ... -> the driver's exit status
    local root="$1" have="$2"; shift 2; mkdir -p "$root/campaign"
    cp "$DRIVER" "$root/campaign/run_thin_atomfind.sh"; cp "$REPO_DIR/campaign/ceos_sweep.tsv" "$root/campaign/"
    local m; for m in lab Pb Ti; do
        case "$have" in sims*) mkdir -p "$root/sim_out_af_armf_a080_${m}${SFX}/01"; : >"$root/sim_out_af_armf_a080_${m}${SFX}/01/data_dp.hdf5";; esac
        case "$have" in *+copy*) local c="${have#*+copy}"; mkdir -p "$root/sim_out_af_armf_a080_${m}${SFX}${c}/01"
                                 : >"$root/sim_out_af_armf_a080_${m}${SFX}${c}/01/data_dp.hdf5"
                                 [ "${PARTIAL:-0}" = 1 ] || : >"$root/sim_out_af_armf_a080_${m}${SFX}${c}/01/rebin_detector.json";; esac
    done
    local u=(-u BASH_ENV -u ENV) v; for v in "${DVARS[@]}"; do u+=(-u "$v"); done
    env "${u[@]}" STUB_LOG="$root/sbatch.log" STUB_N="$root/n" CONDA_ENV="$TMP/env" SHARE= \
        LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=0 FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 "$@" \
        bash "$root/campaign/run_thin_atomfind.sh" >"$root/out.txt" 2>&1
}
fail(){ echo "FAIL  $1"; FAIL=1; }
case_ok(){ # $1 have, $2 data tag, $3 recon tag, $4 expected rebin args ("" = reuse), $5 k (probe up<k>), $6 UPSAMPLING export ("" = none), then VAR=value ...
    local have="$1" dt="$2" rt="$3" rargs="$4" k="$5" up="$6"; shift 6; N=$((N + 1)); local root="$TMP/c$N" m ok=1
    if ! run_driver "$root" "$have" "$@"; then fail "$* : driver refused: $(tail -1 "$root/out.txt")"; return; fi
    local L="$root/sbatch.log"
    local nreb; nreb=$(grep -c -- '--job-name=af_rebin_' "$L"); local nrec; nrec=$(grep -c -- '--job-name=af_rec_' "$L")
    [ "$(grep -c -- '--job-name=af_sim_' "$L")" = 0 ] || { fail "$* : a sim was submitted"; ok=0; }
    [ "$nrec" = 3 ] || { fail "$* : ${nrec} recon jobs, want 3"; ok=0; }
    if [ -n "$rargs" ]; then [ "$nreb" = 3 ] || { fail "$* : ${nreb} rebin jobs, want 3"; ok=0; }
    else [ "$nreb" = 0 ] || { fail "$* : ${nreb} rebin jobs, want 0 (reuse)"; ok=0; }; fi
    for m in lab Pb Ti; do
        local src="sim_out_af_armf_a080_${m}${SFX}" rd="recon_af_armf_a080_${m}${SFX}${rt}_NL18"
        local reb; reb=$(grep -- "--job-name=af_rebin_armf_a080_${m}${SFX}${dt} " "$L")
        if [ -n "$rargs" ]; then
            [[ "$reb" == *"rebin_detector.py' --in-dir '${root}/${src}' --out-dir '${root}/${src}${dt}'${rargs}"* ]] \
                || { fail "$* : ${m} rebin call: ${reb:-<none>}"; ok=0; }
        fi
        local rec; rec=$(grep -- "--job-name=af_rec_armf_a080_${m}${SFX}${rt} " "$L")
        [ -n "$rec" ] || { fail "$* : no recon job af_rec_armf_a080_${m}${SFX}${rt}"; ok=0; continue; }
        local id; id=$(grep -n -- "--job-name=af_rebin_armf_a080_${m}${SFX}${dt} " "$L" | cut -d: -f1)   # stub ids: 1000 + call number
        [ -n "$id" ] && id=$((id + 1000))
        if [ -n "$rargs" ]; then [[ "$rec" == *"--dependency=afterok:${id} "* ]] || { fail "$* : ${m} recon not chained on its rebin ${id}"; ok=0; }
        else [[ "$rec" != *"--dependency"* ]] || { fail "$* : ${m} recon waits on something: ${rec}"; ok=0; }; fi
        local ni="${WANT_NITER:-200}"
        [[ "$rec" == *"NLAYERS=18,SIM_BASE=${root}/${rd}/,REGLAYER=0,PROBE_MODES=1,NITER=${ni},SAVE_EVERY=${ni},GROUPING=32;16,BETA_LSQ=0.05"* ]] \
            || { fail "$* : ${m} recon export: ${rec}"; ok=0; }
        if [ -n "$up" ]; then [[ "$rec" == *",UPSAMPLING=${up} "* ]] || { fail "$* : ${m} no UPSAMPLING=${up}"; ok=0; }
        else [[ "$rec" != *"UPSAMPLING"* ]] || { fail "$* : ${m} UPSAMPLING exported"; ok=0; }; fi
        local f; for f in data_dp.hdf5 data_position.hdf5 sim_meta.mat; do
            [ "$(readlink "$root/$rd/01/$f")" = "${root}/${src}${dt}/01/${f}" ] || { fail "$* : ${m} ${f} -> $(readlink "$root/$rd/01/$f")"; ok=0; }
        done
        [ "$(readlink "$root/$rd/01/probe_initial.mat")" = "${root}/${src}${dt}/01/probe_initial_true_up${k}.mat" ] \
            || { fail "$* : ${m} probe -> $(readlink "$root/$rd/01/probe_initial.mat")"; ok=0; }
    done
    grep -q -- "--job-name=af_pack " "$L" || { fail "$* : no pack job"; ok=0; }
    grep -q "atomfind_results_armf_a080${SFX}${rt}_TS\|atomfind_results_armf_a080${SFX}${rt}_[0-9]" "$root/out.txt" \
        || { fail "$* : tarball not tagged ${SFX}${rt}"; ok=0; }
    [ "$ok" = 1 ] && echo "ok    ${rt}  ($(grep -c . "$L") sbatch calls${rargs:+, rebin${rargs}})  <-  $*"
}
case_chain(){ # a label with no sims yet (the dense scan): sim -> rebin (afterok) -> recon (afterok), all named by the label
    local lab="$1" dt="$2" rt="$3" step="$4"; shift 4; N=$((N + 1)); local root="$TMP/c$N" m ok=1 sid rid
    if ! run_driver "$root" none LABELS="$lab" "$@"; then fail "$* : driver refused: $(tail -1 "$root/out.txt")"; return; fi
    local L="$root/sbatch.log"
    for m in lab Pb Ti; do
        local sd="sim_out_af_${lab}_${m}${SFX}"
        sid=$(grep -n -- "--job-name=af_sim_${lab}_${m}${SFX} " "$L" | cut -d: -f1); [ -n "$sid" ] && sid=$((sid + 1000))
        rid=$(grep -n -- "--job-name=af_rebin_${lab}_${m}${SFX}${dt} " "$L" | cut -d: -f1); [ -n "$rid" ] && rid=$((rid + 1000))
        [ -n "$sid" ] && [ -n "$rid" ] || { fail "$* : ${m} sim ${sid:-none} / rebin ${rid:-none}"; ok=0; continue; }
        local sim; sim=$(grep -- "--job-name=af_sim_${lab}_${m}${SFX} " "$L")
        [[ "$sim" == *"JOB_DIR=${root}/${sd},"*"SCAN_STEP=${step},"*"FOCAL_SPREAD=2.5,SOURCE_FWHM=0.1"* ]] || { fail "$* : ${m} sim export: ${sim}"; ok=0; }
        grep -q -- "--job-name=af_rebin_${lab}_${m}${SFX}${dt} .*--dependency=afterok:${sid} .*--in-dir '${root}/${sd}' --out-dir '${root}/${sd}${dt}'" "$L" \
            || { fail "$* : ${m} rebin not chained on sim ${sid}"; ok=0; }
        grep -q -- "--job-name=af_rec_${lab}_${m}${SFX}${rt} .*--dependency=afterok:${rid} .*SIM_BASE=${root}/recon_af_${lab}_${m}${SFX}${rt}_NL18/" "$L" \
            || { fail "$* : ${m} recon not chained on rebin ${rid}"; ok=0; }
        [ "$(readlink "$root/recon_af_${lab}_${m}${SFX}${rt}_NL18/01/data_dp.hdf5")" = "${root}/${sd}${dt}/01/data_dp.hdf5" ] \
            || { fail "$* : ${m} data link"; ok=0; }
    done
    grep -q "atomfind_results_${lab}${SFX}${rt}_[0-9]" "$root/out.txt" || { fail "$* : tarball not tagged ${lab}${SFX}${rt}"; ok=0; }
    [ "$ok" = 1 ] && echo "ok    ${lab}${rt}  (sim -> rebin -> recon, $(grep -c . "$L") sbatch calls)  <-  $*"
}
case_refused(){ # $1 have, then VAR=value ... that must be refused before any sbatch
    local have="$1"; shift; N=$((N + 1)); local root="$TMP/c$N"
    if run_driver "$root" "$have" "$@"; then fail "$* : accepted"; return; fi
    [ -s "$root/sbatch.log" ] && { fail "$* : refused only after submitting"; return; }
    echo "ok    refused  <-  $*   ($(tail -1 "$root/out.txt"))"
}
echo "driver: ${DRIVER} (${#DVARS[@]} driver variables cleared per case)"
case_ok sims _detc96 _detc96 " --crop-mrad 96 --bin 1" 0 "" RECON_ONLY=1 DET_CROP_MRAD=96                      # R1
case_ok sims _detc96b4 _detc96b4_up2 " --crop-mrad 96 --bin 4" 2 2 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2   # R2
case_ok sims+copy_detc96b4 _detc96b4 _detc96b4 "" 0 "" RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4                  # R3, copy reused
case_ok sims+copy_detc96b4 _detc96b4 _detc96b4_up2 "" 2 2 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2  # R2 again, reused
PARTIAL=1 case_ok sims+copy_detc96b4 _detc96b4 _detc96b4_up2 " --crop-mrad 96 --bin 4" 2 2 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2  # a half-written copy is redone
case_ok sims _detc133b4 _detc133b4_up2 " --crop-mrad 133 --bin 4" 2 2 RECON_ONLY=1 DET_CROP_MRAD=133 DET_BIN=4 UPSAMPLING=2  # R4
case_ok sims _detc96p5 _detc96p5 " --crop-mrad 96.5 --bin 1" 0 "" RECON_ONLY=1 DET_CROP_MRAD=96.5             # '.' -> 'p'
case_ok sims _detc96b4 _detc96b4 " --crop-mrad 96 --bin 4" 0 "" RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=0
# the limits wave (2026-10-26): 512 / 128 px cameras, longer R2 under a run tag on the reused copy, the dense-scan label from scratch
case_ok sims _detc96b2 _detc96b2_up1 " --crop-mrad 96 --bin 2" 1 1 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=2 UPSAMPLING=1   # C512
case_ok sims _detc96b8 _detc96b8_up3 " --crop-mrad 96 --bin 8" 3 3 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=8 UPSAMPLING=3   # C128
WANT_NITER=400 case_ok sims+copy_detc96b4 _detc96b4 _detc96b4_up2_n400 "" 2 2 RECON_ONLY=1 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2 NITER=400 RTAG=_n400
case_chain armf_a080_s04 _detc96b4 _detc96b4_up2 0.4 DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2                       # S04
case_refused sims RECON_ONLY=1 UPSAMPLING=2
case_refused sims RECON_ONLY=1 DET_BIN=2 UPSAMPLING=2
case_refused sims RECON_ONLY=1 DET_BIN=3
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=0
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=-96
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=abc
case_refused sims RECON_ONLY=1 DET_BIN=4 UPSAMPLING=5
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=96 DOSES=1e7
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=96 KICK_ALL=1
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=96 PSEARCH=1
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=96 SEARCH_C1F=1
case_refused sims RECON_ONLY=1 DET_CROP_MRAD=96 KICK_LABELS=armf_a080
case_refused none RECON_ONLY=1 DET_CROP_MRAD=96
[ "$FAIL" = 0 ] && echo "all ${N} cases passed" || { echo "SOME CASES FAILED"; exit 1; }
