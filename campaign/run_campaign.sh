#!/usr/bin/env bash
# [campaign] Parameter-sweep driver for the thin-slab aberration experiments. For every row
# of a sweep .tsv it fans out the proven 3-leg comparison (perfect / ab_known / ab_fitprobe),
# then submits ONE afterany-dependent job that tars every result for a single scp. Built on
# the thin slab + fit-probe-HARD settings de-risked in run_thin_aberration.sh (fast, no NaN).
#
#   CAMPAIGN=round    bash campaign/run_campaign.sh     # alpha sweep  (campaign/round_sweep.tsv, planner-made)
#   CAMPAIGN=nonround bash campaign/run_campaign.sh     # non-round sweep (campaign/nonround_sweep.tsv)
#   RECON_ONLY=1 CAMPAIGN=round bash campaign/run_campaign.sh   # re-fit from existing sims (no re-sim)
#   OVERWRITE=1 ...                                     # let sims clobber an existing sweep
#
# Sweep .tsv columns (TAB-sep; '#'/header/blank skipped): first 8 are read, rest is diagnostics
#   label  alpha  c5  c3  c1  df_perf  bin  aber_json
# c3/c1/c5 = round knobs [Å]; df_perf = aberration-free 4 Å reference defocus [Å]; bin from the
# planner (probe-size -> window). aber_json '-' = round-only (use c3/c5); else a full non-round
# abTEM dict (overrides c3/c5, e.g. '{"C30":-4e4,"C50":1e7,"C56":6e5,"phi56":0}').
#
# The 3 legs per row (each thin, NL=6, reg off, 200 iters):
#   perfect     df=df_perf, no aberration            -> resolution reference at THIS alpha
#   ab_known    true balanced probe fixed            -> ceiling (probe known)
#   ab_fitprobe starts from perfect probe, UPDATEs   -> blind retrieval test (fit probe HARD)
#
# Clean up later: delete campaign/ + the [campaign] blocks in simulate_4dstem.py & run_sim.slurm.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"; mkdir -p logs

CAMPAIGN="${CAMPAIGN:-round}"
TSV="${TSV:-campaign/${CAMPAIGN}_sweep.tsv}"
[ -f "${TSV}" ] || { echo "ERROR: sweep file '${TSV}' not found." >&2; exit 1; }

# thin-slab + fit-probe-HARD defaults (all env-overridable)
# SLICE (sim slab) is FIXED FINE at 0.9 A — below the ~1.94 A atomic-plane spacing along the beam,
# so the ground truth resolves the interatomic depth structure at EVERY alpha and only the recon
# NL (Nyquist per alpha, from the sweep .tsv) decides what's recovered. 2 A merges the planes.
THIN="${THIN:-3}"; SLICE="${SLICE:-0.9}"; STEP="${STEP:-0.5}"
NL_OV="${NL:-}"; NITER="${NITER:-200}"; PSTART="${PSTART:-8}"; BETA="${BETA:-0.1}"; PMODES="${PMODES:-1}"
C5DEF="${C5:-1e7}"; SAVE="${SAVE_EVERY:-25}"
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)               # probe chosen per-leg
TS="$(date +%Y%m%d_%H%M)"; PACK="${SHARE:-$REPO_DIR}/${CAMPAIGN}_results_${TS}.tgz"

# resources per BIN (planner emits bin 4 light / 2 medium / 1 heavy). See HPC_COMMANDS.md.
mem_for(){   case "$1" in 1) echo 175G;;   2) echo 96G;;      *) echo 48G;;      esac; }
grp_for(){   case "$1" in 1) echo 16;;     2) echo 32;;       *) echo "";;       esac; }  # ""=recon default
rtime_for(){ case "$1" in 1) echo 20:00:00;; 2) echo 08:00:00;; *) echo 04:00:00;; esac; }
stime_for(){ case "$1" in 1) echo 12:00:00;; 2) echo 04:00:00;; *) echo 02:00:00;; esac; }

sim_job(){   # $1 dir $2 alpha $3 bin $4 defocus $5 mode(perfect|round|json) $6 aberarg  -> jobid
    local dir="$1" alpha="$2" bin="$3" df="$4" mode="$5" arg="$6"
    local exp="ALL,JOB_DIR=${dir},THIN_CELLS=${THIN},SLICE_THICKNESS=${SLICE},SCAN_STEP=${STEP}"
    exp="${exp},CONVERGENCE=${alpha},BIN_FACTOR=${bin},DEFOCUS=${df},OVERWRITE=${OVERWRITE:-0}"
    export ABERRATIONS_JSON=""                                      # cleared unless json mode (avoids stale carry-over)
    case "$mode" in
        perfect) : ;;                                              # aberration-free: defocus only
        round)   exp="${exp},ABERRATED=1,PROBE_INITIAL=nominal,${arg}";;   # arg = CS=..,C5=..
        json)    export ABERRATIONS_JSON="${arg}"; exp="${exp},PROBE_INITIAL=nominal";;  # JSON via ALL (has commas)
    esac
    sbatch --parsable --job-name="${CAMPAIGN}_sim" --time="$(stime_for "$bin")" \
        --output="logs/${CAMPAIGN}_sim_%j.out" --error="logs/${CAMPAIGN}_sim_%j.err" \
        --export="${exp}" sim/run_sim.slurm
}
recon_job(){ # $1 name $2 datadir $3 probe $4 pstart(""=fixed) $5 bin $6 nl $7 dep(""=now)  -> jobid
    local name="$1" datadir="$2" probe="$3" pstart="$4" bin="$5" nl="$6" dep="$7"
    local rdir="${REPO_DIR}/recon_${CAMPAIGN}_${name}_NL${nl}"
    mkdir -p "${rdir}/01"
    local f; for f in "${INPUTS[@]}"; do ln -sf "${datadir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${probe}" "${rdir}/01/probe_initial.mat"               # start probe (may be another leg's)
    # extra probe modes only help the UPDATING (blind-fit) leg; fixed-probe controls stay 1 mode
    local pm=1; local psx=""; [ -n "${pstart}" ] && { psx=",PROBE_START=${pstart},BETA_LSQ=${BETA}"; pm="${PMODES}"; }
    local grp; grp="$(grp_for "$bin")"; [ -n "$grp" ] && psx="${psx},GROUPING=${grp}"
    local dep_arg=(); [ -n "${dep}" ] && dep_arg=(--dependency="afterok:${dep}")
    sbatch --parsable --job-name="${CAMPAIGN}_rec_${name}" --time="$(rtime_for "$bin")" --mem="$(mem_for "$bin")" \
        ${dep_arg[@]+"${dep_arg[@]}"} \
        --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${nl}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES="${pm}",NITER="${NITER}",SAVE_EVERY="${SAVE}"${psx} \
        run_recon_synthetic_ML.slurm
}

echo "CAMPAIGN=${CAMPAIGN}  tsv=${TSV}  NL=${NL_OV:-per-alpha(Nyquist)} sim_slab=${SLICE}A NITER=${NITER} PSTART=${PSTART} BETA=${BETA} PMODES=${PMODES}"
echo "fit-probe legs update the probe from iter ${PSTART}; results pack to ${PACK}"; echo
RIDS=()                                                             # every recon jobid, for the tar dep
while IFS=$'\t' read -r label alpha c5 c3 c1 dfp bin nl aj _rest; do
    case "$label" in ''|'#'*|label) continue;; esac                # skip blank/comment/header
    bin="${bin:-4}"; c5="${c5:--}"; [ "$c5" = "-" ] && c5="$C5DEF"
    nl="${nl:-6}"; [ "$nl" = "-" ] && nl=6; [ -n "$NL_OV" ] && nl="$NL_OV"   # per-alpha Nyquist unless NL= overrides
    PDIR="${REPO_DIR}/sim_out_${CAMPAIGN}_${label}_perf"
    ADIR="${REPO_DIR}/sim_out_${CAMPAIGN}_${label}_aber"
    # aberrated sim aberration spec: JSON (non-round) overrides the round c3/c5 knobs
    if [ "$aj" != "-" ] && [ -n "$aj" ]; then amode="json"; aarg="$aj"
    else                                       amode="round"; aarg="CS=${c3},C5=${c5}"; fi

    if [ "${RECON_ONLY:-0}" = "1" ]; then
        for d in "$PDIR" "$ADIR"; do [ -e "${d}/01/data_dp.hdf5" ] || {
            echo "ERROR: ${d}/01/data_dp.hdf5 missing — run sims first (unset RECON_ONLY)." >&2; exit 1; }; done
        SP=""; SA=""; depP=""; depK=""; depF=""
    else
        SP=$(sim_job "$PDIR" "$alpha" "$bin" "$dfp"  perfect "")        # perfect: aberration-free @ df_perf
        SA=$(sim_job "$ADIR" "$alpha" "$bin" "$c1"   "$amode" "$aarg")  # aberrated (+probe_initial_true.mat)
        depP="$SP"; depK="$SA"; depF="${SA}:${SP}"
    fi
    R1=$(recon_job "${label}_perfect"  "$PDIR" "${PDIR}/01/probe_initial.mat"      ""         "$bin" "$nl" "$depP")
    R2=$(recon_job "${label}_known"    "$ADIR" "${ADIR}/01/probe_initial_true.mat" ""         "$bin" "$nl" "$depK")
    R3=$(recon_job "${label}_fitprobe" "$ADIR" "${PDIR}/01/probe_initial.mat"      "$PSTART"  "$bin" "$nl" "$depF")
    RIDS+=("$R1" "$R2" "$R3")
    dz=$(awk "BEGIN{printf \"%.2f\", 11.7/${nl}}")
    printf 'row %-9s alpha=%-3s bin=%s NL=%-2s(dz %sA)  sims[P=%s A=%s]  recon[perf=%s known=%s fit=%s]\n' \
        "$label" "$alpha" "$bin" "$nl" "$dz" "${SP:-reuse}" "${SA:-reuse}" "$R1" "$R2" "$R3"
done < "${TSV}"

# one tar job, after ALL recons (afterany: partial sweeps still pack)
DEP=$(IFS=:; echo "${RIDS[*]}")
PJ=$(sbatch --parsable --job-name="${CAMPAIGN}_pack" --time=00:20:00 --mem=8G \
        --dependency="afterany:${DEP}" \
        --output="logs/${CAMPAIGN}_pack_%j.out" --error="logs/${CAMPAIGN}_pack_%j.err" \
        --wrap="bash '${REPO_DIR}/campaign/pack_results.sh' '${CAMPAIGN}' '${PACK}' '${TSV}'")
echo; echo "pack job ${PJ} -> ${PACK} (runs after all ${#RIDS[@]} recons)"
echo "when done:  scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK}' ~/Desktop/"
