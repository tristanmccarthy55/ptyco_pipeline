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
#   ALPHAS="70 90" DOSES="1e7 1e6 1e5 1e4" bash campaign/run_thin_atomfind.sh   # relaxation step 2: shot noise
#   ALPHAS="70 90" PHONONS=16 PER_SPECIES=1 bash campaign/run_thin_atomfind.sh   # step 4: frozen phonons (sim x16)
#   ALPHAS="70 90" THIN=18 CELL_Z=3.889 GROUPING=16 RTIME=20:00:00 bash campaign/run_thin_atomfind.sh   # step 6: 70 A slab
#   TSV=campaign/nonround_sweep.tsv LABELS="nr1_C56_0p6w nr3_C56_2p5w" bash campaign/run_thin_atomfind.sh   # non-round, 70 mrad
# PHONONS>0 re-simulates every leg over that many frozen-phonon configurations (PHONON_SIGMA scalar, or
# PER_SPECIES=1 for the tabulated room-temperature Pb/Sr/Ti/O values); dirs gain _ph<N>. LABELS picks
# rows of the TSV by label instead of alpha, so a non-round row's aber_json (C56, C34, ...) reaches the
# sim as ABERRATIONS_JSON; dirs are named by the label. CELL_Z is the unit cell the driver's box formula
# uses (5 cells: 3.905 reproduces the 27.525 A box; 18 cells: the structure's 3.889 gives 78.00 vs the
# sim's actual 77.94 A). GROUPING / RTIME / STIME override the per-BIN defaults (NL 64 needs the first two).
# DOSES (e/A^2) reuses the existing noiseless sims: per dose and leg a CPU job writes a Poisson copy
# (sim/add_poisson_noise.py -> sim_out_af_a<A>_<mode>_dose<D>) and the recon runs on it once that job
# succeeds -- lab AND Pb/Ti kernels at the same dose (NEXT_PHASE rule 3), independent noise per leg
# (seed = DOSE_SEED + 0/1/2 for lab/Pb/Ti). Recon dirs: recon_af_a<A>_<mode>_dose<D>_NL<NL>.
# Then (when done): extract each PSF and run atomfind (see the echo at the end).
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"; mkdir -p logs
TSV="${TSV:-campaign/round_sweep.tsv}"; [ -f "$TSV" ] || { echo "no $TSV" >&2; exit 1; }
ALPHAS="${ALPHAS:-50 70 90 100}"
THIN="${THIN:-5}"; ZVAC="${ZVAC:-4}"; C5="${C5:-1e7}"; STEP="${STEP:-0.5}"; SLICE="${SLICE:-0.9}"
# WIN defaults to 20 = the sim's own SCAN_WINDOW_A, i.e. the window the lab legs use. The PSF
# kernel is only the matched system response if it comes out of the SAME pipeline as the data;
# the old 14 was a gratuitous difference that also cost the sparse grid its positional diversity.
# GRIDSP=3 (was 4) gives the denser grid that converges at high alpha. WIN=20 + GRIDSP=3 are the
# validated S1 settings (2026-09-11, all 8 kernels clean) -- keep every alpha on them for uniformity.
GRIDSP="${GRIDSP:-3}"; WIN="${WIN:-20}"; NITER="${NITER:-200}"; SAVE="${SAVE_EVERY:-25}"
MODES="${MODES:-lab Pb Ti}"          # e.g. MODES="Pb Ti" to rebuild only the PSF kernels
BETA_LSQ="${BETA_LSQ:-0.05}"         # one LSQ step for every leg (lab and kernels must match)
CELL_Z="${CELL_Z:-3.905}"; LAM=0.0196877
BOXZ=$(awk "BEGIN{printf \"%.3f\", ${THIN}*${CELL_Z}+2*${ZVAC}}")      # full box thickness [Å]
ATOMZ=$(awk "BEGIN{printf \"%.3f\", ${BOXZ}/2}")                        # PSF atom at box centre
INPUTS=(data_dp.hdf5 data_position.hdf5 sim_meta.mat)
DOSES="${DOSES:-}"; DOSE_SEED="${DOSE_SEED:-0}"
PHONONS="${PHONONS:-0}"; PHONON_SIGMA="${PHONON_SIGMA:-0.08}"; PER_SPECIES="${PER_SPECIES:-0}"; PHONON_SEED="${PHONON_SEED:-1}"
LABELS="${LABELS:-}"                      # non-empty: select TSV rows by label (col 1), not alpha
SFX=""; [ "$PHONONS" != 0 ] && SFX="_ph${PHONONS}"; [ "$THIN" != 5 ] && SFX="${SFX}_thin${THIN}"   # sim + recon dir suffix
# (never `VAR=$([ test ] && echo x)`: a false test fails the substitution and set -e exits the script)
PYBIN="${CONDA_ENV:-${SHARE:-}/phucrh/envs/abtem}/bin/python"
# Tarball name to the second, tagged by alphas (+doses): submissions pasted together must not share one.
TS="$(date +%Y%m%d_%H%M%S)"
TAG="$([ -n "$LABELS" ] && echo "$(echo ${LABELS} | tr ' ' '-')" || echo "a$(echo ${ALPHAS} | tr ' ' '-')")${DOSES:+_dose$(echo ${DOSES} | tr ' ' '-')}${SFX}"
PACK="${SHARE:+$SHARE/$USER}"; PACK="${PACK:-$REPO_DIR}/atomfind_results_${TAG}_${TS}.tgz"   # own subdir, not the shared group dir
DIRS_FILE="${REPO_DIR}/logs/af_pack_${TAG}_${TS}.dirs"; : >"${DIRS_FILE}"   # this submission's recon dirs only
WHAT="alphas: ${ALPHAS}"; [ -n "$LABELS" ] && WHAT="labels: ${LABELS}"   # ${X:+a}${X:-b} prints BOTH when X is set
echo "full box ${BOXZ} A (THIN=${THIN} cells + 2x${ZVAC} A vac); PSF atom z=${ATOMZ}; ${WHAT}${SFX:+; suffix ${SFX}}"

if [ -n "$DOSES" ] && [ ! -x "${PYBIN}" ]; then echo "DOSES set but no abtem env python at ${PYBIN}" >&2; exit 1; fi

nl_full(){ awk "BEGIN{n=int(${BOXZ}*2*($1/1000)^2/${LAM}+0.5); if(n<1)n=1; print n}"; }
mem_for(){   case "$1" in 1) echo 175G;; 2) echo 96G;; *) echo 48G;; esac; }
grp_for(){   [ -n "${GROUPING:-}" ] && { echo "$GROUPING"; return; }; case "$1" in 1) echo 16;;  2) echo 32;; *) echo "";; esac; }
rtime_for(){ [ -n "${RTIME:-}" ] && { echo "$RTIME"; return; }; case "$1" in 1) echo 24:00:00;; 2) echo 10:00:00;; *) echo 05:00:00;; esac; }
stime_for(){ [ -n "${STIME:-}" ] && { echo "$STIME"; return; }   # phonons multiply the sim time by PHONONS
             local h; case "$1" in 1) h=12;; 2) h=5;; *) h=3;; esac; [ "$PHONONS" != 0 ] && h=$(( h * 3 )); printf '%02d:00:00' $h; }

sim_job(){   # $1 dir $2 alpha $3 bin $4 c3 $5 c1 $6 mode(lab|Pb|Ti) [$7 aber_json] -> jobid
    local dir="$1" alpha="$2" bin="$3" c3="$4" c1="$5" mode="$6" aj="${7:-}"
    local exp="ALL,JOB_DIR=${dir},SLICE_THICKNESS=${SLICE},SCAN_STEP=${STEP},CONVERGENCE=${alpha}"
    exp="${exp},PHONONS=${PHONONS},PHONON_SIGMA=${PHONON_SIGMA},PER_SPECIES_SIGMA=${PER_SPECIES},PHONON_SEED=${PHONON_SEED}"
    # a non-round row's JSON has commas, so it cannot ride in --export's list: run_sim.slurm reads it from
    # the environment (ALL) instead, as campaign/run_campaign.sh does for its json legs
    if [ -n "$aj" ] && [ "$aj" != "-" ]; then export ABERRATIONS_JSON="$aj"; else unset ABERRATIONS_JSON; fi
    # PROBE_INITIAL=nominal so the sim emits BOTH probe_initial.mat (nominal) AND
    # probe_initial_true.mat (the true aberrated probe). recon_job symlinks the recon's
    # probe_initial.mat -> probe_initial_true.mat = a KNOWN-probe recon. (PROBE_INITIAL=true
    # writes the true probe AS probe_initial.mat and SKIPS probe_initial_true.mat, so the
    # symlink dangles and every recon dies on "File corrupt: probe_initial.mat" -- the
    # 2026-09-09 batch failure. Matches run_campaign.sh's ab_known leg.)
    exp="${exp},BIN_FACTOR=${bin},RECON_FULL_BOX=1,Z_VACUUM=${ZVAC},ABERRATED=1,PROBE_INITIAL=nominal"
    exp="${exp},CS=${c3},C5=${C5},DEFOCUS=${c1},OVERWRITE=${OVERWRITE:-0}"   # OVERWRITE=1 to re-sim over existing
    # SCAN_WINDOW goes to EVERY leg (it used to reach only the grids, leaving the lab on the sim
    # default). The scan field must stay larger than the probe or ptychography loses its positional
    # diversity: scan/d90 is 5.0/5.0/3.0/1.8 at a50-a100 but 0.8 at a110, where the 24.5 A probe
    # exceeds the 20 A field -- and that alpha reconstructs as featureless speckle. Widening it is
    # capped by the 70 A box: with the scan centred at x=40, WIN <= ~35 keeps the d90 core inside.
    exp="${exp},SCAN_WINDOW=${WIN}"
    case "$mode" in
        lab) exp="${exp},THIN_CELLS=${THIN}";;
        *)   exp="${exp},SINGLE_ATOM=${mode},ATOM_Z=${ATOMZ},GRID_SPACING=${GRIDSP},GRID_BOX_Z=${BOXZ}";;
    esac
    sbatch --parsable --job-name="af_sim_${mode}" --time="$(stime_for "$bin")" \
        --output="logs/af_sim_%j.out" --error="logs/af_sim_%j.err" --export="${exp}" sim/run_sim.slurm
}
noise_job(){ # $1 noiseless sim dir $2 noisy out dir $3 dose $4 seed -> jobid  (CPU; streamed, ~minutes)
    local src="$1" out="$2" dose="$3" seed="$4"
    [ -x "${PYBIN}" ] || { echo "no abtem env python at ${PYBIN}" >&2; exit 1; }
    sbatch --parsable --job-name="af_noise" --time=01:00:00 --mem=32G --cpus-per-task=2 \
        --output="logs/af_noise_%j.out" --error="logs/af_noise_%j.err" \
        --wrap="'${PYBIN}' '${REPO_DIR}/sim/add_poisson_noise.py' --in-dir '${src}' --out-dir '${out}' --dose ${dose} --seed ${seed}"
}
recon_job(){ # $1 name $2 datadir $3 bin $4 nl $5 dep -> jobid  (true probe fixed)
    local name="$1" datadir="$2" bin="$3" nl="$4" dep="$5"
    local rdir="${REPO_DIR}/recon_af_${name}_NL${nl}"; mkdir -p "${rdir}/01"
    echo "recon_af_${name}_NL${nl}" >>"${DIRS_FILE}"
    # A previous run's output must not survive into this one: if the new run fails, the OLD
    # *_recons.h5 would be packed and analysed as if new. Moved aside, not deleted; pack_results.sh
    # only matches */analysis/*, so an analysis.prev_* dir is never shipped.
    if [ -d "${rdir}/analysis" ]; then mv "${rdir}/analysis" "${rdir}/analysis.prev_$(date +%Y%m%d_%H%M%S)"; fi
    local f; for f in "${INPUTS[@]}"; do ln -sf "${datadir}/01/${f}" "${rdir}/01/${f}"; done
    ln -sf "${datadir}/01/probe_initial_true.mat" "${rdir}/01/probe_initial.mat"     # known aberrated probe
    local grp; grp="$(grp_for "$bin")"; local gx=""; [ -n "$grp" ] && gx=",GROUPING=${grp}"
    # BETA_LSQ is ALWAYS forwarded, one value for every leg: it is an engine setting, so a lab leg
    # and its kernel must share it. (It is a step size, not a penalty, so a mismatch is far milder
    # than REGLAYER's -- but with a fixed NITER it still moves the result, and on 2026-09-11 the
    # kernels ran at 0.05 against labs at the .m default 0.1. Pinned here so that cannot recur.)
    # 0.05 because a100_Ti diverged to NaN at 0.1. REGLAYER stays 0 for EVERY leg: it symmetrizes
    # information between layers, i.e. low-passes the depth axis. See aberration_experiment/PSF_KERNELS.md.
    gx="${gx},BETA_LSQ=${BETA_LSQ}"
    local dep_arg=(); [ -n "$dep" ] && dep_arg=(--dependency="afterok:${dep}")   # empty dep (RECON_ONLY) -> run now
    sbatch --parsable --job-name="af_rec_${name}" --time="$(rtime_for "$bin")" --mem="$(mem_for "$bin")" \
        ${dep_arg[@]+"${dep_arg[@]}"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${nl}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",SAVE_EVERY="${SAVE}"${gx} \
        run_recon_synthetic_ML.slurm
}

RIDS=()
ROWS="$ALPHAS"; [ -n "$LABELS" ] && ROWS="$LABELS"
for a in $ROWS; do
    if [ -n "$LABELS" ]; then       # by label: alpha and the full aberration JSON come from the row
        read -r alpha c3 c1 bin aj < <(awk -F'\t' -v L="$a" '$1!~/^#/ && $1==L {print $2"\t"$4"\t"$5"\t"$7"\t"$9}' "$TSV")
        [ -n "${bin:-}" ] || { echo "  ${a}: not in $TSV, skipping" >&2; continue; }
        leg="$a"
    else
        read -r c3 c1 bin < <(awk -F'\t' -v A="$a" '$1!~/^#/ && $2==A {print $4"\t"$5"\t"$7}' "$TSV")
        [ -n "${bin:-}" ] || { echo "  a${a}: not in $TSV, skipping" >&2; continue; }
        alpha="$a"; aj="-"; leg="a${a}"
    fi
    nl=$(nl_full "$alpha")
    line="$(printf '%-14s bin=%s NL=%-2s ' "${leg}${SFX}" "$bin" "$nl")"
    for m in $MODES; do                       # MODES="Pb Ti" re-does only the PSF kernels
        D="${REPO_DIR}/sim_out_af_${leg}_${m}${SFX}"
        if [ -n "$DOSES" ]; then               # step 2: Poisson copies of the EXISTING noiseless sim
            [ -e "${D}/01/data_dp.hdf5" ] || { echo "  a${a} ${m}: ${D}/01/data_dp.hdf5 missing -- run the noiseless sim first" >&2; exit 1; }
            case "$m" in lab) so=0;; Pb) so=1;; *) so=2;; esac
            for dose in $DOSES; do
                DN="${D}_dose${dose}"
                N=$(noise_job "$D" "$DN" "$dose" $(( DOSE_SEED + so )))
                R=$(recon_job "${leg}_${m}${SFX}_dose${dose}" "$DN" "$bin" "$nl" "$N")
                RIDS+=("$R"); line+=" ${m}@${dose}=${N}>${R}"
            done
            continue
        fi
        if [ "${RECON_ONLY:-0}" = "1" ]; then  # reuse existing sims; recons run immediately
            [ -e "${D}/01/data_dp.hdf5" ] || { echo "  a${a} ${m}: ${D}/01 missing, skip" >&2; continue; }
            S=""
        else
            S=$(sim_job "$D" "$alpha" "$bin" "$c3" "$c1" "$m" "$aj")
        fi
        R=$(recon_job "${leg}_${m}${SFX}" "$D" "$bin" "$nl" "$S")
        RIDS+=("$R"); line+=" ${m}=${R}"
    done
    echo "$line"
done
[ ${#RIDS[@]} -gt 0 ] || { echo "nothing submitted" >&2; exit 1; }
DEP=$(IFS=:; echo "${RIDS[*]}")
PJ=$(sbatch --parsable --job-name="af_pack" --time=00:30:00 --mem=8G --dependency="afterany:${DEP}" \
    --output="logs/af_pack_%j.out" --error="logs/af_pack_%j.err" --export=ALL,PACK_DIRS_FILE="${DIRS_FILE}" \
    --wrap="bash '${REPO_DIR}/campaign/pack_results.sh' af '${PACK}' '${TSV}'")
echo; echo "pack ${PJ} -> ${PACK}  (recon_af_* .h5 + logs)"
echo "scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK}' ~/Desktop/"
echo "now check the group root stayed clean:  ls /springbrook/share/physics/"
echo "then per alpha:  python analysis/atomfind/extract_psf.py recon_af_a<A>_Pb_NL<NL> Pb_a<A>   (and Ti)"
echo "  -> psf_{Pb,Ti}_a<A>_vol.npy ; point config 'thin' single_atom_vol/ti_kernel_vol at them;"
echo "  atomfind --preset thin --recon recon_af_a<A>_lab_NL<NL>/.../*_recons.h5 --dz \$(bc<<<${BOXZ}/<NL>)"
