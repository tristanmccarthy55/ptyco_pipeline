#!/usr/bin/env bash
# [mono300] Dry run of dose_submit.sh (D) with a stub sbatch (an exported function: nothing is ever submitted). Each case
# runs the submit script in its own throwaway git copy of the files it touches (git, because its guard runs
# campaign/test_driver_legacy.sh against HEAD), with B's 20 K sims as empty files. Checks:
#   - per leg and dose one noise job on B's noiseless sim (seed 0/1/2 for lab/Pb/Ti), no sim job;
#   - per dose the lab leg's 3-stage search (49 / 25 / 49 trials), stage 1 on the lab noise job (afterok), each select on
#     its array, each later stage on the select before it; the draw is C's seed-1 start (search_af_..._s3/start.json);
#   - the refine legs: 8 modes, 250 it, "4;2", 96G, the mixed state on the search's answer, each on its own noise job
#     (afterok) and its dose's last select (afterany); data linked from the noisy copy;
#   - one pack, one analysis per dose with SUFFIX _ph16_T20_fs2p5_src0p1_dose<D>_ps;
#   - refusals before any sbatch: no dose, a dose not written 1e5-1e8, a dose searched before, a sim missing, an exported
#     job variable.
#   bash aberration_experiment/mono300/test_dose_submit.sh      (PYTHON=<python with numpy>: default the hyperspy bundle's)
set -uo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY="${PYTHON:-$HOME/hyperspy-bundle/bin/python}"; [ -x "$PY" ] || PY="${SHARE:-}/phucrh/envs/abtem/bin/python"
"$PY" -c 'import numpy' 2>/dev/null || { echo "no python with numpy (set PYTHON=...): campaign/psearch_plan.py draws the start"; exit 1; }
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
sbatch(){ local n; n=$(( $(cat "$STUB_N" 2>/dev/null || echo 1000) + 1 )); echo "$n" >"$STUB_N"
          printf '%s\n' "$*" >>"$STUB_LOG"; echo "$n"; }
export -f sbatch
S=_ph16_T20_fs2p5_src0p1
FILES=(campaign/run_thin_atomfind.sh campaign/run_analysis.sh campaign/test_coh_names.sh campaign/test_driver_legacy.sh
       campaign/psearch_plan.py campaign/select_grid.py campaign/select_trial.py campaign/search_trial.sh
       run_recon_synthetic_ML.slurm ptycho/run_synthetic_recon_ML.m ptycho/+engines/+GPU/+initialize/init_solver.m
       ptycho/+engines/+GPU/ptycho_solver.m aberration_experiment/page/arm_results.json aberration_experiment/mono300/dose_submit.sh)
# every variable the driver, the analysis driver or the job script reads, cleared so the caller's shell cannot reach a case
DVARS=(); while IFS= read -r v; do DVARS+=("$v"); done < <( { grep -ohE '\$\{[A-Z][A-Z0-9_]*:[-+=?]' "$REPO_DIR/campaign/run_thin_atomfind.sh" \
    "$REPO_DIR/campaign/run_analysis.sh" "$REPO_DIR/run_recon_synthetic_ML.slurm" | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/'
    LC_ALL=C grep -a -ohE "getenv\('[A-Z0-9_]+'\)" "$REPO_DIR/ptycho/run_synthetic_recon_ML.m" | sed -E "s/getenv\('([A-Z0-9_]+)'\)/\1/"; } \
    | grep -vxE 'SHARE|CONDA_ENV|USER|PATH|HOME|TMPDIR' | sort -u)
FAIL=0; N=0
fail(){ echo "FAIL  $1"; FAIL=1; }
make_root(){ # $1 root, $2 legs whose sims exist -> a git copy with B's 20 K sims as empty files
    local root="$1" f m; mkdir -p "$root/env/bin" "$root/share"; ln -s "$PY" "$root/env/bin/python"
    for f in "${FILES[@]}" campaign/*.tsv; do
        (cd "$REPO_DIR" && for g in $f; do mkdir -p "$root/$(dirname "$g")"; cp "$g" "$root/$g"; done)
    done
    (cd "$root" && git init -q && git add -A && git -c user.email=t@t -c user.name=t commit -qm dry) || return 1
    for m in $2; do mkdir -p "$root/sim_out_af_armf_a080_${m}${S}/01"; : >"$root/sim_out_af_armf_a080_${m}${S}/01/data_dp.hdf5"; done
}
run_submit(){ # $1 root, $2 extra env (VAR=value or ""), then the script's arguments -> its screen output in $1/screen.txt
    local root="$1" extra="$2"; shift 2
    local u=(-u BASH_ENV -u ENV) v; for v in "${DVARS[@]}"; do u+=(-u "$v"); done
    (cd / && env "${u[@]}" STUB_LOG="$root/sbatch.log" STUB_N="$root/n" CONDA_ENV="$root/env" SHARE="$root/share" USER=phucrh \
        ${extra:+"$extra"} bash "$root/aberration_experiment/mono300/dose_submit.sh" "$@") >"$root/screen.txt" 2>&1
}
case_refused(){ # $1 what, $2 legs with sims, $3 a dir to make first ("" none), $4 extra env, then the script's arguments
    local what="$1" have="$2" pre="$3" extra="$4"; shift 4; N=$((N + 1)); local root="$TMP/c$N"
    make_root "$root" "$have" || { fail "${what}: could not build the copy"; return; }
    [ -n "$pre" ] && mkdir -p "$root/$pre"
    run_submit "$root" "$extra" "$@"
    local nc=0; [ -f "$root/sbatch.log" ] && nc=$(grep -c . "$root/sbatch.log")
    if grep -q '^NOT SUBMITTED' "$root/screen.txt" && [ "$nc" = 0 ]; then
        echo "ok    refused: ${what}  ($(grep -m1 '^NOT SUBMITTED' "$root/screen.txt" | cut -c1-110))"
    else fail "${what}: not refused (${nc} sbatch calls): $(tail -2 "$root/screen.txt" | tr '\n' ' ')"; fi
}
case_ok(){ # $1 the doses
    local doses="$1"; N=$((N + 1)); local root="$TMP/c$N" ok=1 d m L id sel3 nd
    make_root "$root" "lab Pb Ti" || { fail "${doses}: could not build the copy"; return; }
    run_submit "$root" "" "$doses"; L="$root/sbatch.log"
    grep -q '^ANALYSES:' "$root/screen.txt" || { fail "${doses}: did not reach the end: $(grep -m3 'NOT SUBMITTED\|STOPPED\|FAIL' "$root/screen.txt" | tr '\n' ' ')"; return; }
    nd=$(echo $doses | wc -w | tr -d ' ')
    [ "$(grep -c -- '--job-name=af_sim_' "$L")" = 0 ] || { fail "${doses}: a sim was submitted"; ok=0; }
    [ "$(grep -c -- '--job-name=af_noise' "$L")" = $((3 * nd)) ] || { fail "${doses}: $(grep -c -- '--job-name=af_noise' "$L") noise jobs, want $((3 * nd))"; ok=0; }
    [ "$(grep -c -- '--job-name=af_kick_' "$L")" = $((3 * nd)) ] || { fail "${doses}: $(grep -c -- '--job-name=af_kick_' "$L") refine legs, want $((3 * nd))"; ok=0; }
    [ "$(grep -c -- '--job-name=af_pack' "$L")" = 1 ] || { fail "${doses}: want one pack job"; ok=0; }
    [ "$(grep -c -- '--job-name=af_analysis' "$L")" = "$nd" ] || { fail "${doses}: want ${nd} analysis jobs"; ok=0; }
    for d in $doses; do
        local nm="armf_a080_lab${S}_dose${d}" sd="$root/search_af_armf_a080_lab${S}_dose${d}" so=0 nid
        NID_lab=""; NID_Pb=""; NID_Ti=""
        for m in lab Pb Ti; do   # noise: B's sim -> its dose copy, seed by leg
            nid=$(grep -n -- "add_poisson_noise.py' --in-dir '${root}/sim_out_af_armf_a080_${m}${S}' --out-dir '${root}/sim_out_af_armf_a080_${m}${S}_dose${d}' --dose ${d} --seed ${so}\"\{0,1\}$" "$L" | cut -d: -f1)
            [ -n "$nid" ] || { fail "${d} ${m}: no noise job seed ${so}"; ok=0; so=$((so + 1)); continue; }
            eval "NID_${m}=$((nid + 1000))"; so=$((so + 1))
        done
        local prev="" st arr want
        for st in 1 2 3; do      # the staged search: arrays and selects chained
            case $st in 2) want=24;; *) want=48;; esac
            arr=$(grep -n -- "--array=0-${want} --job-name=af_ps${st}_${nm} " "$L" | cut -d: -f1)
            [ -n "$arr" ] || { fail "${d}: no stage-${st} array of $((want + 1)) trials"; ok=0; continue; }
            if [ "$st" = 1 ]; then grep -q -- "--job-name=af_ps1_${nm} .*--dependency=afterok:${NID_lab} " "$L" || { fail "${d}: stage 1 not on the lab noise job ${NID_lab}"; ok=0; }
            else grep -q -- "--job-name=af_ps${st}_${nm} .*--dependency=afterany:${prev} " "$L" || { fail "${d}: stage ${st} not on select ${prev}"; ok=0; }; fi
            grep -q -- "--job-name=af_psel${st}_${nm} .*--dependency=afterany:$((arr + 1000)) " "$L" || { fail "${d}: select ${st} not on its array"; ok=0; }
            prev=$(( $(grep -n -- "--job-name=af_psel${st}_${nm} " "$L" | cut -d: -f1) + 1000 ))
        done
        sel3="$prev"
        "$PY" -c 'import json, sys; a = json.load(open(sys.argv[1]))["start"]
b = dict(c1_A=-951.6354206191025, c3_A=-403229.6888399238, c12_A=36.757308834245634, phi12=-0.8152211507777685)
sys.exit(0 if all(abs(a[k] - b[k]) <= 1e-6 * max(1.0, abs(b[k])) for k in b) else "draw " + str(a))' "$sd/start.json" \
            || { fail "${d}: the draw is not C's seed-1 start"; ok=0; }
        [ -s "$sd/stage1.tsv" ] && [ "$(( $(wc -l <"$sd/stage1.tsv") - 1 ))" = 49 ] || { fail "${d}: stage1.tsv not 49 trials"; ok=0; }
        local t0; t0=$(awk -F'\t' 'NR == 1 {for (j = 1; j <= NF; j++) if ($j == "recon_dir") c = j; next} NR == 2 {print $c}' "$sd/stage1.tsv")
        [ "$(readlink "$root/$t0/01/data_dp.hdf5")" = "${root}/sim_out_af_armf_a080_lab${S}_dose${d}/01/data_dp.hdf5" ] \
            || { fail "${d}: trial ${t0} reads $(readlink "$root/$t0/01/data_dp.hdf5")"; ok=0; }
        for m in lab Pb Ti; do   # the refine legs
            local rn="armf_a080_${m}${S}_dose${d}_ps" rec nidm; rec=$(grep -- "--job-name=af_kick_${rn} " "$L")
            [ -n "$rec" ] || { fail "${d} ${m}: no refine leg af_kick_${rn}"; ok=0; continue; }
            eval "nidm=\$NID_${m}"
            [[ "$rec" == *"--mem=96G "* ]] || { fail "${d} ${m}: not 96G"; ok=0; }
            [[ "$rec" == *"--dependency=afterok:${nidm},afterany:${sel3} "* ]] || { fail "${d} ${m}: deps $(echo "$rec" | grep -o -- '--dependency=[^ ]*')"; ok=0; }
            [[ "$rec" == *"PROBE_MODES=8,NITER=250,"*"GROUPING=4;2,PROBE_FROM=${root}/search_af_armf_a080_lab${S}_dose${d}/search_best.txt,PROBE_MIXED=8,"* ]] \
                || { fail "${d} ${m}: export $(echo "$rec" | grep -o -- '--export=[^ ]*')"; ok=0; }
            [ "$(readlink "$root/recon_af_${rn}_NL18/01/data_dp.hdf5")" = "${root}/sim_out_af_armf_a080_${m}${S}_dose${d}/01/data_dp.hdf5" ] \
                || { fail "${d} ${m}: leg reads $(readlink "$root/recon_af_${rn}_NL18/01/data_dp.hdf5")"; ok=0; }
            [ -e "$root/recon_af_${rn}_NL18/01/probe_initial.mat" ] && { fail "${d} ${m}: a probe_initial.mat in the leg (the job writes its own)"; ok=0; }
        done
        grep -q -- "--suffix '${S}_dose${d}_ps'" "$root"/logs/analysis_armf_a080${S}_dose${d}_ps_*.sh 2>/dev/null \
            || { fail "${d}: no analysis with suffix ${S}_dose${d}_ps"; ok=0; }
    done
    [ "$ok" = 1 ] && echo "ok    D at ${doses}: $(grep -c . "$L") sbatch calls; draw = C's start; tarball $(grep -o 'atomfind_results_[^ ]*\.tgz' "$root/screen.txt" | head -1 | sed 's/_[0-9]\{8\}_[0-9]\{6\}\.tgz$/_<ts>.tgz/')"
}
case_ok "1e8"
case_ok "1e7 1e6 1e5"
case_refused "no dose given" "lab Pb Ti" "" ""
case_refused "dose 1e9" "lab Pb Ti" "" "" 1e9
case_refused "dose written 1e+08" "lab Pb Ti" "" "" 1e+08
case_refused "1e8 searched before" "lab Pb Ti" "search_af_armf_a080_lab${S}_dose1e8" "" 1e8
case_refused "Ti's 20 K sim missing" "lab Pb" "" "" 1e8
case_refused "KICK_MODES exported" "lab Pb Ti" "" "KICK_MODES=3" 1e8
[ "$FAIL" = 0 ] && echo "all ${N} cases passed" || { echo "SOME CASES FAILED"; exit 1; }
