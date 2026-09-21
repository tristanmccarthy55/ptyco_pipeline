#!/usr/bin/env bash
# [aberration_experiment step 1] C1 (defocus) search with C3/C5 FIXED. Reconstruct the SAME sim
# through a ladder of fixed trial probes that differ only in C1, and read the solver's own
# Fourier-error trace (ptycho/run_synthetic_recon_ML.m writes it as *_error_trace.csv next to the
# h5) as the objective. Per trial: the sim's data + geometry are symlinked into a fresh recon dir,
# the trial probe is written IN THE JOB by sim/make_probe.py (through the sim's own
# build_initial_probe, so the probe grid is the data's by construction), and the recon runs with
# that probe FIXED, on the pinned engine constants (REGLAYER=0, BETA_LSQ=0.05, PROBE_MODES=1).
# Nothing is re-simulated: the sims come from campaign/run_thin_atomfind.sh (sim_out_af_a<A>_<mode>).
#
#   ALPHAS=70 DC1="0"                                   bash campaign/run_c1_search.sh   # smoke test
#   ALPHAS=70 DC1="-60 -50 ... -2 0 0 2 ... 60" NITER=50 bash campaign/run_c1_search.sh   # objective grid
#   ALPHAS=70 C1="-58" MODES="lab Pb Ti" NITER=200      bash campaign/run_c1_search.sh   # final + kernels
#
# DC1 = offsets [A] from the TSV's C1 (a value listed twice runs twice: suffix _r2 -- the objective's
# noise floor); C1 = absolute values instead. MODES lab|Pb|Ti -- the kernel legs get the SAME fitted
# probe (NEXT_PHASE rule 3). PACK_H5=0 packs sidecars/logs/probes but not the h5 (a90: 240 MB each).
# DRYRUN=1 prints every sbatch line and submits nothing. Every path lands under $SHARE/phucrh: recon
# dirs in the repo, logs in each recon dir, the tarball in $SHARE/$USER.
#
# PROBE UPDATE (stage 2.5 -- the experimentalist's route): PSTART=<iter> releases the probe in the
# presolve engine from that iteration, PSTART2=<iter> also in the full engine, with the TEM aperture
# constraint on (PSFFT=1, default when PSTART is set). The trial probe is then only the START: C3/C5 from
# the corrector tableau, C1 from DC1/C1, and the solver refines the whole probe. Campaign defaults to
# c1fit; dirs gain _ps<PSTART>[x<PSTART2>].
#   ALPHAS=70 DC1="-12 -6 0 6 12" PSTART=40 NITER=200 bash campaign/run_c1_search.sh
# The aperture constraint only works since the 2026-09-17 engine fix (load_from_p/init_solver): before
# it the presolve's mask was all zeros, erasing the probe and NaN-ing the next engine.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
TSV="${TSV:-campaign/round_sweep.tsv}"; [ -f "$TSV" ] || { echo "no $TSV" >&2; exit 1; }
ALPHAS="${ALPHAS:-70 90}"; MODES="${MODES:-lab}"; NITER="${NITER:-50}"
# one intermediate Niter*.mat per engine (at the end): at BIN 2 each holds a ~230 MB object, and a grid is 25 trials
SAVE="${SAVE_EVERY:-$NITER}"
THIN="${THIN:-5}"; ZVAC="${ZVAC:-4}"; BETA_LSQ="${BETA_LSQ:-0.05}"
PACK_H5="${PACK_H5:-1}"; DRYRUN="${DRYRUN:-0}"
PSTART="${PSTART:-}"; PSTART2="${PSTART2:-}"; PDFO="${PDFO:-0}"
# PDFO=1: DEFOCUS-ONLY probe update (engine option probe_defocus_only, 2026-09-21): C3/C5 fixed as in the
# start probe, the update moves defocus alone, so the aperture constraint is redundant (PSFFT default 0).
[ "$PDFO" = 1 ] && [ -z "$PSTART" ] && { echo "PDFO=1 needs PSTART (when the probe is released)" >&2; exit 1; }
PSFFT="${PSFFT:-$([ -n "$PSTART" ] && [ "$PDFO" != 1 ] && echo 1 || echo 0)}"
[ -z "$PSTART" ] && [ -n "$PSTART2" ] && { echo "PSTART2 needs PSTART" >&2; exit 1; }
# (no `$([ test ] && echo x)` in assignments: a false test fails the substitution and set -e exits the script)
DFO_TAG=""; [ "$PDFO" = 1 ] && DFO_TAG="dfo"
PS_TAG=""; [ -n "$PSTART" ] && PS_TAG="_ps${PSTART}${PSTART2:+x${PSTART2}}${DFO_TAG}"
CAMP="${CAMP:-$([ "$PDFO" = 1 ] && echo c1dfo || ([ -n "$PSTART" ] && echo c1fit || echo c1))}"
SIM_ROOT="${SIM_ROOT:-$REPO_DIR}"        # where sim_out_af_a<A>_<mode> live (override for a local dry run)
CELL_Z=3.905; LAM=0.0196877
BOXZ=$(awk "BEGIN{printf \"%.3f\", ${THIN}*${CELL_Z}+2*${ZVAC}}")      # full box thickness [A]
# aberrations.json + probe_initial_true.mat are linked for the WRITER (C3/C5 tableau; self-check
# against the truth when the trial C1 is the true one) -- never as probe_initial.mat.
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat aberrations.json probe_initial_true.mat)
# Tarball name carries alphas, NITER and the second: several submissions pasted together must not share
# (and overwrite) one tarball.
TS="$(date +%Y%m%d_%H%M%S)"; TAG="a$(echo ${ALPHAS} | tr ' ' '-')_n${NITER}${PS_TAG}"
PACK="${SHARE:+$SHARE/$USER}"; PACK="${PACK:-$REPO_DIR}/${CAMP}_results_${TAG}_${TS}.tgz"
DIRS_FILE="${REPO_DIR}/logs/${CAMP}_pack_${TAG}_${TS}.dirs"; RDIRS=()
[ -n "${DC1:-}" ] || [ -n "${C1:-}" ] || { echo "set DC1 (offsets from the TSV C1) or C1 (absolute values) [A]" >&2; exit 1; }
[ -n "${DC1:-}" ] && [ -n "${C1:-}" ] && { echo "set DC1 or C1, not both" >&2; exit 1; }
echo "C1 search [${CAMP}]: alphas ${ALPHAS}; modes ${MODES}; NITER ${NITER}; full box ${BOXZ} A; ${DC1:+dC1 = ${DC1}}${C1:+C1 = ${C1}}"
if [ -n "$PSTART" ]; then echo "  probe UPDATE${DFO_TAG:+ (DEFOCUS ONLY)}: presolve from iter ${PSTART}, full engine from ${PSTART2:-never (fixed)}; aperture constraint ${PSFFT}"
else echo "  probe FIXED at each trial C1"; fi

nl_full(){ awk "BEGIN{n=int(${BOXZ}*2*($1/1000)^2/${LAM}+0.5); if(n<1)n=1; print n}"; }
mem_for(){ case "$1" in 1) echo 175G;; 2) echo 96G;; *) echo 48G;; esac; }
grp_for(){ case "$1" in 1) echo 16;;  2) echo 32;; *) echo "";; esac; }
time_for(){ # $1 bin $2 niter -> HH:MM:SS : 15 min startup + NITER x (presolve+full s/iter) x 2.
            # s/iter from the packed logs: a70 (BIN 4) 0.7+2.2, a90 (BIN 2) 4.3+17.7.
    local s; case "$1" in 1) echo 24:00:00; return;; 2) s=22;; *) s=3;; esac
    [ -n "$PSTART" ] && s=$(( s * 3 / 2 ))          # probe update adds work per iteration
    local m=$(( 15 + ($2 * s * 2 + 59) / 60 )); printf '%02d:%02d:00' $((m/60)) $((m%60)); }

preflight_sim(){ # $1 aberrations.json $2 alpha $3 c3 $4 c1 -- the sim on disk must be the TSV's probe.
                 # (An older tarball's a70 sim carried C3 -4 um / C1 +2 A while the TSV said -5 um / -60.)
    local j="$1"; [ -f "$j" ] || { echo "PREFLIGHT FAIL: $j missing -- the sim predates aberrations.json; re-sim" >&2; exit 1; }
    local ja jc3 jc1
    ja=$(sed -n 's/.*"convergence_mrad": *\([-0-9.eE+]*\).*/\1/p' "$j" | head -1)
    jc3=$(sed -n 's/.*"C30": *\([-0-9.eE+]*\).*/\1/p' "$j" | head -1)
    jc1=$(sed -n 's/.*"defocus_A": *\([-0-9.eE+]*\).*/\1/p' "$j" | head -1)
    awk -v a="${ja:-nan}" -v b="$2" -v c="${jc3:-nan}" -v d="$3" -v e="${jc1:-nan}" -v f="$4" \
        'BEGIN{ if (a=="nan"||c=="nan"||e=="nan"||(a-b)^2>1e-6||(c-d)^2>1e-3||(e-f)^2>1e-6) exit 1 }' \
        || { echo "PREFLIGHT FAIL: $j is alpha=${ja} C3=${jc3} C1=${jc1} but ${TSV} says alpha=$2 C3=$3 C1=$4" >&2
             echo "                -> stale sim; re-run campaign/run_thin_atomfind.sh for it (OVERWRITE=1)" >&2; exit 1; }
}

recon_job(){ # $1 name $2 datadir $3 bin $4 nl $5 c1 $6 c3 $7 c5 -> jobid  (trial probe written in-job, then FIXED)
    local name="$1" datadir="$2" bin="$3" nl="$4" c1="$5" c3="$6" c5="$7"
    local rdir="${REPO_DIR}/recon_${CAMP}_${name}_NL${nl}"
    local grp; grp="$(grp_for "$bin")"; local gx=""; [ -n "$grp" ] && gx=",GROUPING=${grp}"
    if [ -n "$PSTART" ]; then
        gx="${gx},PROBE_START=${PSTART},PROBE_SUPPORT_FFT=${PSFFT},PROBE_DEFOCUS_ONLY=${PDFO}"; [ -n "$PSTART2" ] && gx="${gx},PROBE_START2=${PSTART2}"
    fi
    local cmd=(sbatch --parsable --job-name="${CAMP}_${name}" --time="$(time_for "$bin" "$NITER")" --mem="$(mem_for "$bin")"
               --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err"
               --export=ALL,NLAYERS="${nl}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",SAVE_EVERY="${SAVE}",BETA_LSQ="${BETA_LSQ}",PROBE_C1="${c1}",PROBE_C3="${c3}",PROBE_C5="${c5}"${gx}
               run_recon_synthetic_ML.slurm)
    echo "recon_${CAMP}_${name}_NL${nl}" >>"${DIRS_FILE}.pending"
    if [ "$DRYRUN" = 1 ]; then echo "DRYRUN: ${cmd[*]}" >&2; echo "dry"; return; fi
    mkdir -p "${rdir}/01"
    # A previous run's output must not survive into this one (it would be packed and analysed as new).
    if [ -d "${rdir}/analysis" ]; then mv "${rdir}/analysis" "${rdir}/analysis.prev_$(date +%Y%m%d_%H%M%S)"; fi
    local f; for f in "${INPUTS[@]}"; do ln -sf "${datadir}/01/${f}" "${rdir}/01/${f}"; done
    # No probe may pre-exist here: the in-job writer's output must be the ONLY probe_initial.mat
    # (it refuses to overwrite a symlink, so a stale link would kill the job -- correctly).
    rm -f "${rdir}/01/probe_initial.mat" "${rdir}/01/probe_initial.json"
    "${cmd[@]}"
}

# Pass 1 validates EVERY (alpha, mode) before anything is submitted: a stale or missing sim found at the
# third leg must not leave the first two already queued.
PLAN=()
for a in $ALPHAS; do
    read -r c5t c3 c1t bin < <(awk -F'\t' -v A="$a" '$1!~/^#/ && $2==A {print $3"\t"$4"\t"$5"\t"$7}' "$TSV")
    [ -n "${bin:-}" ] || { echo "  a${a}: not in $TSV" >&2; exit 1; }
    for m in $MODES; do
        D="${SIM_ROOT}/sim_out_af_a${a}_${m}"
        if [ ! -e "${D}/01/data_dp.hdf5" ]; then
            if [ "$DRYRUN" = 1 ]; then echo "  (dry run: ${D}/01/data_dp.hdf5 missing locally, continuing)" >&2
            else echo "  a${a} ${m}: ${D}/01/data_dp.hdf5 missing -- run campaign/run_thin_atomfind.sh first" >&2; exit 1; fi
        fi
        [ -e "${D}/01/probe_initial_true.mat" ] || echo "  a${a} ${m}: no probe_initial_true.mat -- the writer's self-check is skipped" >&2
        preflight_sim "${D}/01/aberrations.json" "$a" "$c3" "$c1t"
        PLAN+=("${a}|${m}|${D}|${c5t}|${c3}|${c1t}|${bin}")
    done
done

# Pass 2 submits.
RIDS=()
for row in "${PLAN[@]}"; do
    IFS='|' read -r a m D c5t c3 c1t bin <<<"$row"
    nl=$(nl_full "$a")
    list=""
    if [ -n "${C1:-}" ]; then list="$C1"
    else for d in $DC1; do list+=" $(awk "BEGIN{printf \"%g\", ${c1t}+(${d})}")"; done; fi
    done_list=""
    for c1 in $list; do
        c1s=$(printf '%g' "$c1")
        n=$(( $(printf '%s\n' ${done_list} | grep -cx -- "$c1s" || true) + 1 )); done_list+=" ${c1s}"
        rep=""; [ "$n" -gt 1 ] && rep="_r${n}"
        name="a${a}_${m}_df${c1s}${rep}${PS_TAG}_n${NITER}"
        R=$(recon_job "$name" "$D" "$bin" "$nl" "$c1s" "$c3" "$c5t")
        RIDS+=("$R")
        printf '  %-28s C1=%-7s (dC1 %+g) NL=%-2s bin=%s -> %s\n' "$name" "$c1s" \
            "$(awk "BEGIN{print ${c1s}-(${c1t})}")" "$nl" "$bin" "$R"
    done
done
[ ${#RIDS[@]} -gt 0 ] || { echo "nothing submitted" >&2; exit 1; }
if [ "$DRYRUN" = 1 ]; then
    echo "DRYRUN: ${#RIDS[@]} jobs would be submitted; pack -> ${PACK}; dirs:"; sed 's/^/    /' "${DIRS_FILE}.pending"
    rm -f "${DIRS_FILE}.pending"; exit 0
fi
mv "${DIRS_FILE}.pending" "${DIRS_FILE}"
DEP=$(IFS=:; echo "${RIDS[*]}")
PJ=$(sbatch --parsable --job-name="${CAMP}_pack" --time=00:30:00 --mem=8G --dependency="afterany:${DEP}" \
    --output="logs/${CAMP}_pack_%j.out" --error="logs/${CAMP}_pack_%j.err" \
    --export=ALL,PACK_H5="${PACK_H5}",PACK_DIRS_FILE="${DIRS_FILE}" \
    --wrap="bash '${REPO_DIR}/campaign/pack_results.sh' ${CAMP} '${PACK}' '${TSV}'")
echo; echo "${#RIDS[@]} recon jobs; pack ${PJ} -> ${PACK}  (PACK_H5=${PACK_H5})"
echo "scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK}' ~/Desktop/"
echo "now check the group root stayed clean:  ls /springbrook/share/physics/"
echo "then locally:  python analysis/c1_objective.py --root <extracted dir>"
