#!/usr/bin/env bash
# [atomfind] Full-box thin recons + MATCHED aberrated PSFs across the ROUND alpha sweep, so
# atomfind can show depth localisation IMPROVING with alpha (and breaking once the probe blows
# up, >~100 mrad). Per alpha, reading (c3,c1,bin) from campaign/round_sweep.tsv:
#   lab : a KNOWN-probe labyrinth recon, THIN=5 cells + Z_VACUUM=4 A + --recon-full-box, so the
#         entrance/exit surface artifacts dump into the vacuum bands and the atomic planes stay
#         clean (the interface fix). This is the atom-finding target.
#   psfPb/psfTi : a Pb / Ti GRID through the SAME aberrated probe + same full box -> each blob is
#         the MATCHED system PSF (with the aberrations), for extract_psf.py -> atomfind's kernels.
# NL = Nyquist over the FULL padded box (thick/(depthres/2), depthres=lambda/alpha^2). All legs use
# the true aberrated probe fixed (no blind fit here — atomfind wants the best recon, probe known).
#
#   bash campaign/run_thin_atomfind.sh                 # alphas 50 70 90 100 (feasible BIN<=2)
#   ALPHAS="50 70 90 100 110 120" bash campaign/run_thin_atomfind.sh   # + the heavy BIN=1 break
# Then (when done): extract each PSF and run atomfind (see the echo at the end).
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
TSV="${TSV:-campaign/round_sweep.tsv}"; [ -f "$TSV" ] || { echo "no $TSV" >&2; exit 1; }
ALPHAS="${ALPHAS:-50 70 90 100}"
THIN="${THIN:-5}"; ZVAC="${ZVAC:-4}"; C5="${C5:-1e7}"; STEP="${STEP:-0.5}"; SLICE="${SLICE:-0.9}"
GRIDSP="${GRIDSP:-4}"; WIN="${WIN:-14}"; NITER="${NITER:-200}"; SAVE="${SAVE_EVERY:-25}"
CELL_Z=3.905; LAM=0.0196877
BOXZ=$(awk "BEGIN{printf \"%.3f\", ${THIN}*${CELL_Z}+2*${ZVAC}}")      # full box thickness [Å]
ATOMZ=$(awk "BEGIN{printf \"%.3f\", ${BOXZ}/2}")                        # PSF atom at box centre
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)
TS="$(date +%Y%m%d_%H%M)"; PACK="${SHARE:-$REPO_DIR}/atomfind_results_${TS}.tgz"
echo "full box ${BOXZ} A (THIN=${THIN} cells + 2x${ZVAC} A vac); PSF atom z=${ATOMZ}; alphas: ${ALPHAS}"

nl_full(){ awk "BEGIN{n=int(${BOXZ}*2*($1/1000)^2/${LAM}+0.5); if(n<1)n=1; print n}"; }
mem_for(){   case "$1" in 1) echo 175G;; 2) echo 96G;; *) echo 48G;; esac; }
grp_for(){   case "$1" in 1) echo 16;;  2) echo 32;; *) echo "";; esac; }
rtime_for(){ case "$1" in 1) echo 24:00:00;; 2) echo 10:00:00;; *) echo 05:00:00;; esac; }
stime_for(){ case "$1" in 1) echo 12:00:00;; 2) echo 05:00:00;; *) echo 03:00:00;; esac; }

sim_job(){   # $1 dir $2 alpha $3 bin $4 c3 $5 c1 $6 mode(lab|Pb|Ti) -> jobid
    local dir="$1" alpha="$2" bin="$3" c3="$4" c1="$5" mode="$6"
    local exp="ALL,JOB_DIR=${dir},SLICE_THICKNESS=${SLICE},SCAN_STEP=${STEP},CONVERGENCE=${alpha}"
    # PROBE_INITIAL=nominal so the sim emits BOTH probe_initial.mat (nominal) AND
    # probe_initial_true.mat (the true aberrated probe). recon_job symlinks the recon's
    # probe_initial.mat -> probe_initial_true.mat = a KNOWN-probe recon. (PROBE_INITIAL=true
    # writes the true probe AS probe_initial.mat and SKIPS probe_initial_true.mat, so the
    # symlink dangles and every recon dies on "File corrupt: probe_initial.mat" -- the
    # 2026-09-09 batch failure. Matches run_campaign.sh's ab_known leg.)
    exp="${exp},BIN_FACTOR=${bin},RECON_FULL_BOX=1,Z_VACUUM=${ZVAC},ABERRATED=1,PROBE_INITIAL=nominal"
    exp="${exp},CS=${c3},C5=${C5},DEFOCUS=${c1},OVERWRITE=${OVERWRITE:-0}"   # OVERWRITE=1 to re-sim over existing
    case "$mode" in
        lab) exp="${exp},THIN_CELLS=${THIN}";;
        *)   exp="${exp},SINGLE_ATOM=${mode},ATOM_Z=${ATOMZ},GRID_SPACING=${GRIDSP},SCAN_WINDOW=${WIN},GRID_BOX_Z=${BOXZ}";;
    esac
    sbatch --parsable --job-name="af_sim_${mode}" --time="$(stime_for "$bin")" \
        --output="logs/af_sim_%j.out" --error="logs/af_sim_%j.err" --export="${exp}" sim/run_sim.slurm
}
recon_job(){ # $1 name $2 datadir $3 bin $4 nl $5 dep -> jobid  (true probe fixed)
    local name="$1" datadir="$2" bin="$3" nl="$4" dep="$5"
    local rdir="${REPO_DIR}/recon_af_${name}_NL${nl}"; mkdir -p "${rdir}/01"
    local f; for f in "${INPUTS[@]}"; do ln -sf "${datadir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${datadir}/01/probe_initial_true.mat" "${rdir}/01/probe_initial.mat"     # known aberrated probe
    local grp; grp="$(grp_for "$bin")"; local gx=""; [ -n "$grp" ] && gx=",GROUPING=${grp}"
    local dep_arg=(); [ -n "$dep" ] && dep_arg=(--dependency="afterok:${dep}")   # empty dep (RECON_ONLY) -> run now
    sbatch --parsable --job-name="af_rec_${name}" --time="$(rtime_for "$bin")" --mem="$(mem_for "$bin")" \
        ${dep_arg[@]+"${dep_arg[@]}"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${nl}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",SAVE_EVERY="${SAVE}"${gx} \
        run_recon_synthetic_ML.slurm
}

RIDS=()
for a in $ALPHAS; do
    read -r c3 c1 bin < <(awk -F'\t' -v A="$a" '$1!~/^#/ && $2==A {print $4"\t"$5"\t"$7}' "$TSV")
    [ -n "${bin:-}" ] || { echo "  a${a}: not in $TSV, skipping" >&2; continue; }
    nl=$(nl_full "$a")
    LD="${REPO_DIR}/sim_out_af_a${a}_lab"; PD="${REPO_DIR}/sim_out_af_a${a}_Pb"; TD="${REPO_DIR}/sim_out_af_a${a}_Ti"
    if [ "${RECON_ONLY:-0}" = "1" ]; then     # reuse existing sims (no re-sim); recons run immediately
        for d in "$LD" "$PD" "$TD"; do [ -e "${d}/01/data_dp.hdf5" ] || { echo "  a${a}: ${d}/01 missing, skip" >&2; continue 2; }; done
        SL=""; SP=""; ST=""
    else
        SL=$(sim_job "$LD" "$a" "$bin" "$c3" "$c1" lab)
        SP=$(sim_job "$PD" "$a" "$bin" "$c3" "$c1" Pb)
        ST=$(sim_job "$TD" "$a" "$bin" "$c3" "$c1" Ti)
    fi
    R1=$(recon_job "a${a}_lab" "$LD" "$bin" "$nl" "$SL")
    R2=$(recon_job "a${a}_Pb"  "$PD" "$bin" "$nl" "$SP")
    R3=$(recon_job "a${a}_Ti"  "$TD" "$bin" "$nl" "$ST")
    RIDS+=("$R1" "$R2" "$R3")
    printf 'a%-3s bin=%s NL=%-2s  lab=%s Pb=%s Ti=%s\n' "$a" "$bin" "$nl" "$R1" "$R2" "$R3"
done
DEP=$(IFS=:; echo "${RIDS[*]}")
PJ=$(sbatch --parsable --job-name="af_pack" --time=00:20:00 --mem=8G --dependency="afterany:${DEP}" \
    --output="logs/af_pack_%j.out" --error="logs/af_pack_%j.err" \
    --wrap="bash '${REPO_DIR}/campaign/pack_results.sh' af '${PACK}' '${TSV}'")
echo; echo "pack ${PJ} -> ${PACK}  (recon_af_* .h5 + logs)"
echo "scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK}' ~/Desktop/"
echo "then per alpha:  python analysis/atomfind/extract_psf.py recon_af_a<A>_Pb_NL<NL> Pb_a<A>   (and Ti)"
echo "  -> psf_{Pb,Ti}_a<A>_vol.npy ; point config 'thin' single_atom_vol/ti_kernel_vol at them;"
echo "  atomfind --preset thin --recon recon_af_a<A>_lab_NL<NL>/.../*_recons.h5 --dz \$(bc<<<${BOXZ}/<NL>)"
