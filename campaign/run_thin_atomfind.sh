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
# KICK_LABELS (2026-10-01): those labels (must be in LABELS) also get a "kick" leg on their lab simulation -- the probe
# RELEASED from a start probe built from the true tableau with C1, C3 and A1 scaled by KICK (default 1.05; written
# in-job by sim/make_probe.py), KICK_MODES probe modes (3), probe update from KICK_PSTART in the presolve and KICK_PSTART2
# in the full engine (40 / 20), the TEM aperture constraint on, GROUPING KICK_GROUPING ("8;4": modes multiply the GPU
# memory), walltime KICK_RTIME (48:00:00). Dir recon_af_<label>_kick_NL<n>; analyse_sweep.py compares its probe with the
# true one (kick/) and runs atomfind on it with the label's known-probe kernels. Part of the same pack, so CLEANDATA
# waits for it.
# PARTIAL COHERENCE (2026-10-02): FOCAL_SPREAD (rms A) or ENERGY_SPREAD_EV + CC_MM, SOURCE_FWHM (A) and COHERENCE_SAMPLES
# reach the sim (sim/simulate_4dstem.py: sampled jointly with PHONONS); dirs gain the values, e.g. _fs3p5, _src0p4,
# _dE0p09_cc2p4 (the cold-FEG setting below keeps its old _coh; see COH_EXP). With DOSES the chain is built in
# ONE submission: sim (when its data are missing) -> Poisson copy -> recon. KICK_ALL=1 runs EVERY leg (lab, Pb, Ti) as a
# kick leg -- start probe kicked, probe released, KICK_MODES modes -- so the kernels come out of the same reconstruction
# operator as the lab leg (PSF_KERNELS.md: byte-identical except the object).
# OUTER SEARCH + REFINED KICK (2026-10-04, the experiment-like route; probe update ON): SEARCH_C1F / SEARCH_C3F (factors on
# the kicked C1 / C3) run a job array of fixed-probe trials on the lab data, campaign/select_trial.py picks the start, and
# every leg runs as a kick leg from it (_rk). FIXED_TWIN=1 adds true-probe-fixed legs (_fixed); RTAG tags recon dirs;
# KICK_MODES / KICK_NITER / KICK_GROUPING take lists (one kick leg per entry, _m<modes>). Details above search_job().
#   LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv PHONONS=16 PER_SPECIES=1 ENERGY_SPREAD_EV=0.3 CC_MM=1.4 SOURCE_FWHM=0.4 \
#     DOSES="1e7 1e8" FIXED_TWIN=1 SEARCH_C1F="0.91 0.94 0.97 1 1.03 1.06 1.09" SEARCH_C3F="<same>" KICK_MODES=4 \
#     KICK_NITER=500 GROUPING="32;16" bash campaign/run_thin_atomfind.sh
# A coherence scan point (noiseless, true probe fixed; dirs sim_out_af_armf_a080_<mode>_fs2p5):
#   LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 FOCAL_SPREAD=2.5 CLEANDATA=1 bash campaign/run_thin_atomfind.sh
# Its probe-mode leg from the sim's own MIXED STATE (KICK_START=mixed, details at KICK_START below), on the kept sim:
#   RECON_ONLY=1 LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv FOCAL_SPREAD=5 KICK_ALL=1 KICK_START=mixed KICK_MODES=7 \
#     KICK_NITER=250 KICK_GROUPING="4;2" RTAG=_mix PACK_H5=0 bash campaign/run_thin_atomfind.sh
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
GRIDSP="${GRIDSP:-3}"; WIN="${WIN:-20}"; NITER="${NITER:-200}"
# One save, at the end. The solver wrote a full object copy every 25 iterations by default: ~60 MB a
# time at NL 14, eight per engine, two engines per leg, and pack_results.sh has never packed one of
# them. Across the campaign that was 347 GB of checkpoints nobody reads, and on 2026-09-22 it helped
# fill a 3.9 TiB SHARED share and took down other people's jobs. SAVE_EVERY=<n> if a mid-run object
# is genuinely wanted.
SAVE="${SAVE_EVERY:-$NITER}"
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
# PHONON_T=<K>: a cooled specimen (simulate_4dstem --phonon-temperature; needs PER_SPECIES=1), dirs gain _T<K>
PHONON_T="${PHONON_T:-}"
if [ -n "$PHONON_T" ]; then
    [ "$PHONONS" != 0 ] && [ "$PER_SPECIES" = 1 ] || { echo "PHONON_T needs PHONONS>0 and PER_SPECIES=1" >&2; exit 1; }
    SFX="${SFX}_T${PHONON_T}"
fi
# [coherence] partial-coherence settings for the sim. The dir suffix carries their VALUES (2026-10-12), so a scan over
# focal spread or source size gets one sim dir per point -- the bare _coh it replaced put every setting in ONE dir:
#   _fs<rms A> (FOCAL_SPREAD)  _dE<eV>_cc<mm> (ENERGY_SPREAD_EV + CC_MM)  _src<FWHM A> (SOURCE_FWHM)  _cs<n> (COHERENCE_SAMPLES)
# in that order, each value as %g with '.' written 'p' (3.50 -> 3p5): the engine's MATLAB fileparts reads a '.' in a
# path as an extension. The one setting every run before 2026-10-12 used, ENERGY_SPREAD_EV=0.3 CC_MM=1.4 SOURCE_FWHM=0.4
# (the ARM200F's cold FEG: 7.29 A rms + a 0.4 A source), keeps its old name _coh, so those dirs, tarballs and page tags
# still resolve. Settings the sim would silently ignore are refused, so a name cannot claim physics the data lack.
COH_EXP=""
for v in FOCAL_SPREAD ENERGY_SPREAD_EV CC_MM SOURCE_FWHM COHERENCE_SAMPLES; do
    [ -n "${!v:-}" ] || continue
    [[ "${!v}" =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]] && awk -v x="${!v}" 'BEGIN{exit !(x > 0)}' \
        || { echo "${v}=${!v}: a positive number, or leave it unset for none" >&2; exit 1; }
    COH_EXP="${COH_EXP},${v}=${!v}"
done
coh_num(){ awk -v x="$1" 'BEGIN{s = sprintf("%g", x); gsub(/\./, "p", s); print s}'; }
if [ -n "$COH_EXP" ]; then
    # the sim forms the focal spread from ENERGY_SPREAD_EV only WITH CC_MM, and FOCAL_SPREAD overrides both
    { [ -n "${ENERGY_SPREAD_EV:-}" ] && [ -n "${CC_MM:-}" ]; } || { [ -z "${ENERGY_SPREAD_EV:-}" ] && [ -z "${CC_MM:-}" ]; } \
        || { echo "ENERGY_SPREAD_EV and CC_MM go together (one alone gives no focal spread)" >&2; exit 1; }
    [ -n "${FOCAL_SPREAD:-}" ] && [ -n "${ENERGY_SPREAD_EV:-}" ] \
        && { echo "FOCAL_SPREAD or ENERGY_SPREAD_EV + CC_MM, not both (the sim would take FOCAL_SPREAD)" >&2; exit 1; }
    if [ -n "${COHERENCE_SAMPLES:-}" ]; then   # the sim takes PHONONS configurations when phonons are on
        [ "$PHONONS" = 0 ] || { echo "COHERENCE_SAMPLES is ignored with PHONONS>0 (configurations = PHONONS)" >&2; exit 1; }
        [ -n "${FOCAL_SPREAD:-}${ENERGY_SPREAD_EV:-}${SOURCE_FWHM:-}" ] \
            || { echo "COHERENCE_SAMPLES alone does nothing: set a focal spread or a source size" >&2; exit 1; }
    fi
    if [ -z "${FOCAL_SPREAD:-}${COHERENCE_SAMPLES:-}" ] && [ -n "${ENERGY_SPREAD_EV:-}" ] && [ -n "${SOURCE_FWHM:-}" ] \
       && awk -v e="$ENERGY_SPREAD_EV" -v c="$CC_MM" -v s="$SOURCE_FWHM" 'BEGIN{exit !(e == 0.3 && c == 1.4 && s == 0.4)}'; then
        SFX="${SFX}_coh"                                  # the pre-2026-10-12 name, kept for its runs
    else
        [ -n "${FOCAL_SPREAD:-}" ] && SFX="${SFX}_fs$(coh_num "$FOCAL_SPREAD")"
        [ -n "${ENERGY_SPREAD_EV:-}" ] && SFX="${SFX}_dE$(coh_num "$ENERGY_SPREAD_EV")_cc$(coh_num "$CC_MM")"
        [ -n "${SOURCE_FWHM:-}" ] && SFX="${SFX}_src$(coh_num "$SOURCE_FWHM")"
        [ -n "${COHERENCE_SAMPLES:-}" ] && SFX="${SFX}_cs$(coh_num "$COHERENCE_SAMPLES")"
    fi
fi
KICK_ALL="${KICK_ALL:-0}"
# (never `VAR=$([ test ] && echo x)`: a false test fails the substitution and set -e exits the script)
PYBIN="${CONDA_ENV:-${SHARE:-}/phucrh/envs/abtem}/bin/python"
# DETECTOR (2026-10-23, aberration_experiment/detector256): record each leg's sim on a smaller camera, without re-simulating.
# DET_CROP_MRAD keeps +- that angle, DET_BIN sums DET_BIN x DET_BIN pixels (a CPU job, sim/rebin_detector.py, writes
# sim_out_af_<leg>_<mode><SFX>_det[c<crop>][b<bin>] beside the sim, which is never touched or cleaned); every leg (lab, Pb,
# Ti alike: the kernels see the same camera) then reconstructs that copy with the true probe fixed, on the grid the solver
# models (probe_initial_true_up<k>.mat). UPSAMPLING=k has the solver model the summed pixels (p.detector.upsampling,
# 2^k = DET_BIN; recon dirs gain _up<k>); without it a binned leg models only the detector's own (2^k smaller) window.
#   RECON_ONLY=1 LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 FOCAL_SPREAD=2.5 SOURCE_FWHM=0.1 \
#     DET_CROP_MRAD=96 DET_BIN=4 UPSAMPLING=2 bash campaign/run_thin_atomfind.sh     # -> ..._fs2p5_src0p1_detc96b4_up2
# Known-probe legs only for now: refused with DOSES, kicks or searches.
DET_CROP_MRAD="${DET_CROP_MRAD:-}"; DET_BIN="${DET_BIN:-1}"; UPSAMPLING="${UPSAMPLING:-}"
DTAG=""; UTAG=""
if [ -n "$DET_CROP_MRAD" ] || [ "$DET_BIN" != 1 ] || [ -n "$UPSAMPLING" ]; then
    case "$DET_BIN" in 1|2|4|8) ;; *) echo "DET_BIN=${DET_BIN}: 1, 2, 4 or 8 (the engine upsamples by 2^k)" >&2; exit 1;; esac
    if [ -n "$DET_CROP_MRAD" ]; then
        [[ "$DET_CROP_MRAD" =~ ^[0-9]*\.?[0-9]+$ ]] && awk -v x="$DET_CROP_MRAD" 'BEGIN{exit !(x > 0)}' \
            || { echo "DET_CROP_MRAD=${DET_CROP_MRAD}: a positive angle in mrad" >&2; exit 1; }
    fi
    [ -n "$DET_CROP_MRAD" ] || [ "$DET_BIN" != 1 ] \
        || { echo "UPSAMPLING=${UPSAMPLING} needs a binned detector (DET_BIN)" >&2; exit 1; }
    case "${UPSAMPLING:-0}" in
        0) ;;
        1|2|3) [ $(( 1 << UPSAMPLING )) = "$DET_BIN" ] \
                   || { echo "UPSAMPLING=${UPSAMPLING} models $(( 1 << UPSAMPLING )) x $(( 1 << UPSAMPLING )) sums: needs DET_BIN=$(( 1 << UPSAMPLING ))" >&2; exit 1; }
               UTAG="_up${UPSAMPLING}";;
        *) echo "UPSAMPLING=${UPSAMPLING}: 0-3, with 2^UPSAMPLING = DET_BIN" >&2; exit 1;;
    esac
    [ -z "${DOSES}${KICK_LABELS:-}${SEARCH_C1F:-}${SEARCH_C3F:-}" ] && [ "$KICK_ALL" = 0 ] && [ "${PSEARCH:-0}" = 0 ] \
        || { echo "DET_CROP_MRAD / DET_BIN: known-probe legs only (not with DOSES, KICK_*, SEARCH_* or PSEARCH)" >&2; exit 1; }
    [ -x "${PYBIN}" ] || { echo "DET_CROP_MRAD / DET_BIN set but no abtem env python at ${PYBIN} (sim/rebin_detector.py)" >&2; exit 1; }
    DTAG="_det"; [ -n "$DET_CROP_MRAD" ] && DTAG="${DTAG}c$(coh_num "$DET_CROP_MRAD")"; [ "$DET_BIN" != 1 ] && DTAG="${DTAG}b${DET_BIN}"
fi
# Tarball name to the second, tagged by alphas (+doses): submissions pasted together must not share one.
TS="$(date +%Y%m%d_%H%M%S)"
TAG="$([ -n "$LABELS" ] && echo "$(echo ${LABELS} | tr ' ' '-')" || echo "a$(echo ${ALPHAS} | tr ' ' '-')")${DOSES:+_dose$(echo ${DOSES} | tr ' ' '-')}${SFX}${DTAG}${UTAG}${RTAG:-}"   # RTAG: a re-run's tarball is told apart too
PACK="${SHARE:+$SHARE/$USER}"; PACK="${PACK:-$REPO_DIR}/atomfind_results_${TAG}_${TS}.tgz"   # own subdir, not the shared group dir
DIRS_FILE="${REPO_DIR}/logs/af_pack_${TAG}_${TS}.dirs"; : >"${DIRS_FILE}"   # this submission's recon dirs only
WHAT="alphas: ${ALPHAS}"; [ -n "$LABELS" ] && WHAT="labels: ${LABELS}"   # ${X:+a}${X:-b} prints BOTH when X is set
echo "full box ${BOXZ} A (THIN=${THIN} cells + 2x${ZVAC} A vac); PSF atom z=${ATOMZ}; ${WHAT}${SFX:+; suffix ${SFX}}"
if [ -n "$DTAG" ]; then
    echo "detector: ${DET_CROP_MRAD:+crop +-${DET_CROP_MRAD} mrad, }${DET_BIN}x${DET_BIN} sums, upsampling ${UPSAMPLING:-0}; data ${DTAG}, recon dirs ${DTAG}${UTAG}"
fi

if [ -n "$DOSES" ] && [ ! -x "${PYBIN}" ]; then echo "DOSES set but no abtem env python at ${PYBIN}" >&2; exit 1; fi

# NL is Nyquist over the full box for that aperture. NL=<n> overrides it, for ONE purpose: testing
# whether a leg's residual is limited by the reconstruction's own slice thickness rather than by the
# physics under test. The sim slices at 0.9 A and the recon at 27.525/NL (1.97 A at NL 14), so the
# forward models differ by ~2x and the solver can never drive the residual to zero. A leg re-run at
# double NL that drops sharply was slice-limited; one that does not was not. Dirs carry the NL, so an
# override cannot overwrite the Nyquist run. (2026-09-22)
nl_full(){ [ -n "${NL:-}" ] && { echo "$NL"; return; }
           awk "BEGIN{n=int(${BOXZ}*2*($1/1000)^2/${LAM}+0.5); if(n<1)n=1; print n}"; }
# --mem (2026-10-08): a recon peaks at ~4.3x its raw patterns -- measured, 1305041 (4 modes, 1420 px, 1600 positions, 500 it):
# MaxRSS 55 GB on 12.9 GB of data -- and Blythe allocates CPUs in proportion to --mem (175G -> 24 CPUs against the gpu QOS's
# 160 CPUs per user: ~6 recons at once). Region rows therefore ask for 1.3 x 4.3 x their own data (at least 32G): ~72G at
# 1420 px, ~85G at 1536 px, about twice as many running at once. RMEM=<n>G overrides; rows without a region keep the table.
mem_for(){ [ -n "${RMEM:-}" ] && { echo "$RMEM"; return; }
           if [ -z "${ROW_SIDE:-}" ]; then case "$1" in 1) echo 175G;; 2) echo 96G;; *) echo 48G;; esac; return; fi
           awk -v s="$ROW_SIDE" -v b="$bin" -v t="${ROW_DETMAX:-200}" -v w="${ROW_WIN:-$WIN}" -v st="${ROW_STEP:-$STEP}" -v l="$LAM" \
               'BEGIN{n = 2*t/1000*(s/b)/l; np = int(w/st + 1e-9); np = np*np; g = int(1.3*4.3*np*n*n*4/1e9 + 0.999)
                      if (g < 32) g = 32; printf "%dG", g}'; }
grp_for(){   [ -n "${GROUPING:-}" ] && { echo "$GROUPING"; return; }; case "$1" in 1) echo 16;;  2) echo 32;; *) echo "";; esac; }
rtime_for(){ [ -n "${RTIME:-}" ] && { echo "$RTIME"; return; }; case "$1" in 1) echo 24:00:00;; 2) echo 10:00:00;; *) echo 05:00:00;; esac; }
# [region] rows carry their own box side / scan / recorded angle (cols 11-14 of the tsv: side win step detmax).
# In a region box BIN no longer says how big a pattern is, so resources follow N = 2 theta_max (side/BIN) / lambda:
# <= 400 px the BIN-4 class, <= 800 the BIN-2 class, else the BIN-1 class (175G, GROUPING 16, 24 h).
res_class(){ [ -z "${ROW_SIDE:-}" ] && { echo "$1"; return; }
             awk -v s="$ROW_SIDE" -v b="$1" -v t="${ROW_DETMAX:-200}" -v l="$LAM" \
                 'BEGIN{n=2*t/1000*(s/b)/l; print (n<=400)?4:((n<=800)?2:1)}'; }
stime_for(){ [ -n "${STIME:-}" ] && { echo "$STIME"; return; }   # phonons multiply the sim time by PHONONS
             local h; case "$1" in 1) h=12;; 2) h=5;; *) h=3;; esac; [ "$PHONONS" != 0 ] && h=$(( h * 3 )); printf '%02d:00:00' $h; }

sim_job(){   # $1 dir $2 alpha $3 bin $4 c3 $5 c1 $6 mode(lab|Pb|Ti) [$7 aber_json] -> jobid
    local dir="$1" alpha="$2" bin="$3" c3="$4" c1="$5" mode="$6" aj="${7:-}"
    # the job name carries the leg: sacct lists names, not directories, and with every sim in a
    # nine-label submission called "af_sim_lab" a failure cannot be attributed to a leg. That cost a
    # round trip on 2026-09-22 working out which simulation a failed recon had been waiting on.
    local jn="af_sim_$(basename "$dir" | sed 's/^sim_out_af_//')"
    local exp="ALL,JOB_DIR=${dir},SLICE_THICKNESS=${SLICE},SCAN_STEP=${ROW_STEP:-$STEP},CONVERGENCE=${alpha}"
    exp="${exp},PHONONS=${PHONONS},PHONON_SIGMA=${PHONON_SIGMA},PER_SPECIES_SIGMA=${PER_SPECIES},PHONON_SEED=${PHONON_SEED}${PHONON_T:+,PHONON_TEMPERATURE=${PHONON_T}}${COH_EXP}"
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
    # C5 from the row in LABELS mode (col 3): an ARM row (4 mm) and a CEOS row (1 mm) can share a submission. The JSON
    # carries C50 and wins in the simulator either way; this keeps the export (and the sim log) honest.
    exp="${exp},CS=${c3},C5=${ROW_C5:-$C5},DEFOCUS=${c1},OVERWRITE=${OVERWRITE:-0}"   # OVERWRITE=1 to re-sim over existing
    # SCAN_WINDOW goes to EVERY leg (it used to reach only the grids, leaving the lab on the sim
    # default). The scan field must stay larger than the probe or ptychography loses its positional
    # diversity: scan/d90 is 5.0/5.0/3.0/1.8 at a50-a100 but 0.8 at a110, where the 24.5 A probe
    # exceeds the 20 A field -- and that alpha reconstructs as featureless speckle. Widening it is
    # capped by the 70 A box: with the scan centred at x=40, WIN <= ~35 keeps the d90 core inside.
    exp="${exp},SCAN_WINDOW=${ROW_WIN:-$WIN}"
    [ -n "${ROW_SIDE:-}" ] && exp="${exp},REGION_SIDE=${ROW_SIDE}"
    [ -n "${ROW_DETMAX:-}" ] && exp="${exp},DETECTOR_MAX_ANGLE=${ROW_DETMAX}"
    case "$mode" in
        lab) exp="${exp},THIN_CELLS=${THIN}";;
        *)   exp="${exp},SINGLE_ATOM=${mode},ATOM_Z=${ATOMZ},GRID_SPACING=${GRIDSP},GRID_BOX_Z=${BOXZ}";;
    esac
    local sb="$bin"; [ -n "${ROW_SIDE:-}" ] && sb=1          # a region box is simulated on a large grid
    sbatch --parsable --job-name="${jn}" --time="$(stime_for "$sb")" \
        --output="logs/af_sim_%j.out" --error="logs/af_sim_%j.err" --export="${exp}" sim/run_sim.slurm
}
noise_job(){ # $1 noiseless sim dir $2 noisy out dir $3 dose $4 seed [$5 sim job to wait for] -> jobid  (CPU; ~minutes)
    local src="$1" out="$2" dose="$3" seed="$4" dep="${5:-}"
    [ -x "${PYBIN}" ] || { echo "no abtem env python at ${PYBIN}" >&2; exit 1; }
    local dep_arg=(); [ -n "$dep" ] && dep_arg=(--dependency="afterok:${dep}")
    sbatch --parsable --job-name="af_noise" --time=01:00:00 --mem=32G --cpus-per-task=2 ${dep_arg[@]+"${dep_arg[@]}"} \
        --output="logs/af_noise_%j.out" --error="logs/af_noise_%j.err" \
        --wrap="'${PYBIN}' '${REPO_DIR}/sim/add_poisson_noise.py' --in-dir '${src}' --out-dir '${out}' --dose ${dose} --seed ${seed}"
}
rebin_job(){ # $1 sim dir $2 out dir [$3 sim job to wait for] -> jobid  (CPU; minutes: sim/rebin_detector.py, DETECTOR above)
    local src="$1" out="$2" dep="${3:-}"
    local dep_arg=(); [ -n "$dep" ] && dep_arg=(--dependency="afterok:${dep}")
    sbatch --parsable --job-name="af_rebin_$(basename "$out" | sed 's/^sim_out_af_//')" --time=01:00:00 --mem=32G --cpus-per-task=2 \
        ${dep_arg[@]+"${dep_arg[@]}"} --output="logs/af_rebin_%j.out" --error="logs/af_rebin_%j.err" \
        --wrap="'${PYBIN}' '${REPO_DIR}/sim/rebin_detector.py' --in-dir '${src}' --out-dir '${out}'${DET_CROP_MRAD:+ --crop-mrad ${DET_CROP_MRAD}} --bin ${DET_BIN}"
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
    # known aberrated probe; on a smaller camera (DETECTOR above) the copy's true probe on the grid the solver models
    local pf="probe_initial_true.mat"; [ -n "$DTAG" ] && pf="probe_initial_true_up${UPSAMPLING:-0}.mat"
    ln -sf "${datadir}/01/${pf}" "${rdir}/01/probe_initial.mat"
    local cls; cls="$(res_class "$bin")"
    local grp; grp="$(grp_for "$cls")"; local gx=""; [ -n "$grp" ] && gx=",GROUPING=${grp}"
    # BETA_LSQ is ALWAYS forwarded, one value for every leg: it is an engine setting, so a lab leg
    # and its kernel must share it. (It is a step size, not a penalty, so a mismatch is far milder
    # than REGLAYER's -- but with a fixed NITER it still moves the result, and on 2026-09-11 the
    # kernels ran at 0.05 against labs at the .m default 0.1. Pinned here so that cannot recur.)
    # 0.05 because a100_Ti diverged to NaN at 0.1. REGLAYER stays 0 for EVERY leg: it symmetrizes
    # information between layers, i.e. low-passes the depth axis. See aberration_experiment/PSF_KERNELS.md.
    gx="${gx},BETA_LSQ=${BETA_LSQ}"
    [ -n "$UTAG" ] && gx="${gx},UPSAMPLING=${UPSAMPLING}"          # the solver models the summed pixels (DETECTOR above)
    local dep_arg=(); [ -n "$dep" ] && dep_arg=(--dependency="afterok:${dep}")   # empty dep (RECON_ONLY) -> run now
    sbatch --parsable --job-name="af_rec_${name}" --time="$(rtime_for "$cls")" --mem="$(mem_for "$cls")" \
        ${dep_arg[@]+"${dep_arg[@]}"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${nl}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES=1,NITER="${NITER}",SAVE_EVERY="${SAVE}"${gx} \
        run_recon_synthetic_ML.slurm
}

KICK_LABELS="${KICK_LABELS:-}"; KICK="${KICK:-1.05}"; KICK_MODES="${KICK_MODES:-3}"; KICK_PSTART="${KICK_PSTART:-40}"
KICK_PSTART2="${KICK_PSTART2:-20}"; KICK_GROUPING="${KICK_GROUPING:-8;4}"; KICK_RTIME="${KICK_RTIME:-48:00:00}"
KICK_NITER="${KICK_NITER:-}"
# KICK_START=mixed (2026-10-14): every kick leg starts from the MIXED STATE of its own sim instead of the tableau probe --
# the top KICK_MODES eigenmodes of the focal spread / source size in its sim_meta.mat, on the true probe
# (sim/make_mixed_probe.py, run in the job via PROBE_MIXED), so no mode starts at the 2 % seed. Needs a coherence setting;
# not with an outer search (that one chooses C1 / C3 of a coherent start). KICK is then unused.
KICK_START="${KICK_START:-tableau}"
case "$KICK_START" in tableau) ;; mixed) [ -n "$COH_EXP" ] \
    || { echo "KICK_START=mixed needs a coherence setting (FOCAL_SPREAD / ENERGY_SPREAD_EV+CC_MM / SOURCE_FWHM)" >&2; exit 1; };;
    *) echo "KICK_START=${KICK_START}: tableau or mixed" >&2; exit 1;; esac
# KICK_MODES / KICK_NITER / KICK_GROUPING may be LISTS (2026-10-04), e.g. KICK_MODES="3 6" KICK_NITER="500 250"
# KICK_GROUPING="8;4 4;2": one kick leg per KICK_MODES entry, tagged _m<modes> when there is more than one (a single entry
# keeps the old names); a missing KICK_NITER / KICK_GROUPING entry falls back to NITER / the first grouping. Time scales
# with the modes (80 mrad, 1420 px, NL 18: 11.2 h per 200 iterations at 3 modes, 14.7 h at 4).
read -r -a KM <<<"$KICK_MODES"; read -r -a KN <<<"$KICK_NITER"; read -r -a KG <<<"$KICK_GROUPING"
dep_arg_for(){ # a job id (-> afterok) or a whole dependency expression ("afterok:1,afterany:2") -> the --dependency option
    case "$1" in "") ;; *[!0-9]*) echo "--dependency=$1";; *) echo "--dependency=afterok:$1";; esac; }
kick_job(){ # $1 name $2 sim dir $3 bin $4 nl $5 dep $6 c1 $7 c3 $8 modes $9 niter $10 grouping [$11 start file] -> jobid
    # kicked start probe, probe RELEASED. With a start file (an outer search's search_best.txt, see SEARCH below) the start
    # probe's C1/C3 are read from it IN THE JOB (PROBE_FROM in run_recon_synthetic_ML.slurm) instead of c1/c3 x KICK
    local name="$1" datadir="$2" bin="$3" nl="$4" dep="$5" c1="$6" c3="$7" modes="$8" niter="$9" grp="${10}" from="${11:-}"
    local rdir="${REPO_DIR}/recon_af_${name}_NL${nl}"; mkdir -p "${rdir}/01"
    echo "recon_af_${name}_NL${nl}" >>"${DIRS_FILE}"
    if [ -d "${rdir}/analysis" ]; then mv "${rdir}/analysis" "${rdir}/analysis.prev_$(date +%Y%m%d_%H%M%S)"; fi
    # the writer needs aberrations.json (the tableau) and the true probe (its overlap record); the in-job probe must be
    # the only probe_initial.mat (make_probe.py refuses to write through a symlink)
    local f; for f in "${INPUTS[@]}" aberrations.json probe_initial_true.mat; do ln -sf "${datadir}/01/${f}" "${rdir}/01/${f}"; done
    rm -f "${rdir:?}/01/probe_initial.mat" "${rdir:?}/01/probe_initial.json"
    local start k1 k3
    if [ -n "$from" ]; then
        start="PROBE_FROM=${from}"     # with KICK_START=mixed the mixed state is built on that (searched) probe in the job
        if [ "$KICK_START" = mixed ]; then start="${start},PROBE_MIXED=${modes}"; fi
    elif [ "$KICK_START" = mixed ]; then start="PROBE_MIXED=${modes}"
    else k1=$(awk "BEGIN{printf \"%.4f\", ${c1}*${KICK}}"); k3=$(awk "BEGIN{printf \"%.2f\", ${c3}*${KICK}}"); start="PROBE_C1=${k1},PROBE_C3=${k3}"; fi
    local cls da; cls="$(res_class "$bin")"; da=$(dep_arg_for "$dep")
    # PROBE_SCALE rides in the environment (its value holds '='), as ABERRATIONS_JSON does for the sims
    PROBE_SCALE="C12=${KICK}" sbatch --parsable --job-name="af_kick_${name}" --time="${KICK_RTIME}" --mem="$(mem_for "$cls")" \
        ${da:+"$da"} --output="${rdir}/slurm_%j.out" --error="${rdir}/slurm_%j.err" \
        --export=ALL,NLAYERS="${nl}",SIM_BASE="${rdir}/",REGLAYER=0,PROBE_MODES="${modes}",NITER="${niter}",SAVE_EVERY="${SAVE_EVERY:-$niter}",BETA_LSQ="${BETA_LSQ}",GROUPING="${grp}",${start},PROBE_START="${KICK_PSTART}",PROBE_START2="${KICK_PSTART2}",PROBE_SUPPORT_FFT=1 \
        run_recon_synthetic_ML.slurm
}
kick_variants(){ # $1 name $2 sim dir $3 bin $4 nl $5 dep $6 c1 $7 c3 [$8 start file] -> job ids, space separated
    local i tag ids=""
    for i in "${!KM[@]}"; do
        tag=""; [ ${#KM[@]} -gt 1 ] && tag="_m${KM[$i]}"
        ids="${ids} $(kick_job "$1${tag}" "$2" "$3" "$4" "$5" "$6" "$7" "${KM[$i]}" "${KN[$i]:-$NITER}" "${KG[$i]:-${KG[0]}}" "${8:-}")"
    done
    echo $ids
}

# OUTER SEARCH (2026-10-04; the user: the experiment-like route is an outer search over the knobs into a probe refinement
# on the reconstruction -- probe update ON). SEARCH_C1F / SEARCH_C3F are factors on the KICKED start's C1 / C3 (c x KICK x f;
# one axis alone is allowed). Per label, and per dose with DOSES, the lab leg gets ONE job array of fixed-probe trials, one
# per (f_c1, f_c3) pair, SEARCH_NITER iterations each (50 ranked trials as 200 did, step 1 2026-09-18), on the very data the
# legs reconstruct; then campaign/select_trial.py (CPU, afterany on the array) writes search_af_<leg>/search_best.{json,txt}:
# the minimum of the trials' final error (a quadratic about the best grid point). Then EVERY leg in MODES (lab first; Pb and
# Ti take the lab leg's answer: one probe per instrument) runs as a kick leg from that start, tagged _rk, with KICK_MODES /
# KICK_NITER / KICK_GROUPING; A1 stays at its kicked tableau value. select always writes its files (with too few trials it
# falls back to the kicked start and says so) and a refined leg whose start file is missing fails at once, so a broken chain
# still reaches the pack. SEARCH_MEM: default the leg's own --mem (mem_for, ~72G at 1420 px).
# SEARCH_RTIME 04:00:00 (50 iterations ~1.2 h at 80 mrad).
# FIXED_TWIN=1: every leg ALSO runs with the true probe fixed, tagged _fixed -- the side control ("possible at all?").
# RTAG=<tag> is appended to every recon dir name (never the sims), so a re-run cannot land in an earlier run's dirs.
SEARCH_C1F="${SEARCH_C1F:-}"; SEARCH_C3F="${SEARCH_C3F:-}"; SEARCH_NITER="${SEARCH_NITER:-50}"
SEARCH_MEM="${SEARCH_MEM:-}"; SEARCH_RTIME="${SEARCH_RTIME:-04:00:00}"; FIXED_TWIN="${FIXED_TWIN:-0}"; RTAG="${RTAG:-}"
SEARCH=0
if [ -n "${SEARCH_C1F}${SEARCH_C3F}" ]; then
    SEARCH=1; SEARCH_C1F="${SEARCH_C1F:-1}"; SEARCH_C3F="${SEARCH_C3F:-1}"
    [ "${MODES%% *}" = lab ] || { echo "SEARCH needs MODES to start with lab (Pb and Ti take the lab leg's answer)" >&2; exit 1; }
    [ -x "${PYBIN}" ] || { echo "SEARCH set but no abtem env python at ${PYBIN} (campaign/select_trial.py)" >&2; exit 1; }
fi
search_job(){ # $1 name (the lab leg) $2 sim dir $3 bin $4 nl $5 dep $6 c1 $7 c3 -> "<array id> <select id> <start file>"
    local name="$1" datadir="$2" bin="$3" nl="$4" dep="$5" c1="$6" c3="$7"
    local sd="search_af_${name}"; mkdir -p "${REPO_DIR}/${sd}"; echo "${sd}" >>"${DIRS_FILE}"
    local man="${REPO_DIR}/${sd}/trials.tsv" i=0 f1 f3 k1 k3 rd f
    printf 'idx\tf_c1\tf_c3\tc1_A\tc3_A\trecon_dir\n' >"$man"
    for f1 in $SEARCH_C1F; do for f3 in $SEARCH_C3F; do
        k1=$(awk "BEGIN{printf \"%.4f\", ${c1}*${KICK}*${f1}}"); k3=$(awk "BEGIN{printf \"%.2f\", ${c3}*${KICK}*${f3}}")
        rd="recon_af_${name}_t${i}_NL${nl}"; mkdir -p "${REPO_DIR}/${rd}/01"; echo "${rd}" >>"${DIRS_FILE}"
        if [ -d "${REPO_DIR}/${rd}/analysis" ]; then mv "${REPO_DIR}/${rd}/analysis" "${REPO_DIR}/${rd}/analysis.prev_$(date +%Y%m%d_%H%M%S)"; fi
        for f in "${INPUTS[@]}" aberrations.json probe_initial_true.mat; do ln -sf "${datadir}/01/${f}" "${REPO_DIR}/${rd}/01/${f}"; done
        rm -f "${REPO_DIR:?}/${rd:?}/01/probe_initial.mat" "${REPO_DIR:?}/${rd:?}/01/probe_initial.json"
        printf '%d\t%s\t%s\t%s\t%s\t%s\n' "$i" "$f1" "$f3" "$k1" "$k3" "$rd" >>"$man"; i=$((i + 1))
    done; done
    local cls grp da A S; cls="$(res_class "$bin")"; grp="$(grp_for "$cls")"; da=$(dep_arg_for "$dep")
    A=$(PROBE_SCALE="C12=${KICK}" sbatch --parsable --array="0-$((i - 1))" --job-name="af_srch_${name}" --time="${SEARCH_RTIME}" \
        --mem="${SEARCH_MEM:-$(mem_for "$cls")}" ${da:+"$da"} \
        --output="${REPO_DIR}/recon_af_${name}_t%a_NL${nl}/slurm_%A_%a.out" --error="${REPO_DIR}/recon_af_${name}_t%a_NL${nl}/slurm_%A_%a.err" \
        --export=ALL,SEARCH_MANIFEST="${man}",NLAYERS="${nl}",REGLAYER=0,PROBE_MODES=1,NITER="${SEARCH_NITER}",SAVE_EVERY="${SEARCH_NITER}",BETA_LSQ="${BETA_LSQ}"${grp:+,GROUPING=${grp}} \
        campaign/search_trial.sh); A="${A%%;*}"
    S=$(sbatch --parsable --job-name="af_select_${name}" --time=00:30:00 --mem=8G --cpus-per-task=1 --dependency="afterany:${A}" \
        --output="${REPO_DIR}/${sd}/slurm_select_%j.out" --error="${REPO_DIR}/${sd}/slurm_select_%j.err" \
        --wrap="'${PYBIN}' '${REPO_DIR}/campaign/select_trial.py' --manifest '${man}' --out '${REPO_DIR}/${sd}'")
    echo "${A} ${S} ${REPO_DIR}/${sd}/search_best.txt"
}
# STAGED OUTER SEARCH FROM A DRAWN START (2026-10-16; the user: C1, C3 and A1, the start drawn from the corrector's quoted
# precision, the search bounded by it). PSEARCH=1: per lab leg (and dose) campaign/psearch_plan.py draws the start (seed
# PSEARCH_SEED; 95 % intervals PSEARCH_CI95="<C1> <C3> <A1>" in A, default the ARM run-3 tableau's 39.5 8040 59) and lays
# out stage 1 (C1 x C3, PSEARCH_N1 per axis, +- ci95 about the start, A1 at the start) and stage 2 (the A1 vector,
# PSEARCH_N2 per axis, +- ci95 per component, at stage 1's answer). Each stage is a job array of fixed-probe trials
# (campaign/search_trial.sh, SEARCH_NITER iterations) then a campaign/select_grid.py job; every leg in MODES then runs as a
# kick leg from stage 2's answer (search_best.txt "C1 C3 C12 phi12"), tagged _ps (Pb and Ti take the lab leg's answer); with
# KICK_START=mixed the mixed state is built on that probe. Higher orders stay at their true values. Not with SEARCH_C1F/C3F.
# PSEARCH_STAGES=3 (the default from 2026-10-20) adds stage 3, C1 x C3 again at stage 2's A1, and the legs start from
# stage 3's answer: with A1 still far off, stage 1's C1 x C3 minimum is biased (the first run: C1 32 A off), while A1 is
# robust to a C1 error. PSEARCH_RESUME=<an earlier search_af_<leg> dir with the same draw> takes the stages that dir
# finished (its stage<k>_best.txt) instead of re-running them; the draw must match (start.json) or nothing is submitted.
PSEARCH="${PSEARCH:-0}"; PSEARCH_SEED="${PSEARCH_SEED:-1}"; PSEARCH_CI95="${PSEARCH_CI95:-39.5 8040 59}"
PSEARCH_N1="${PSEARCH_N1:-7}"; PSEARCH_N2="${PSEARCH_N2:-5}"; PSEARCH_STAGES="${PSEARCH_STAGES:-3}"; PSEARCH_RESUME="${PSEARCH_RESUME:-}"
if [ "$PSEARCH" = 1 ]; then
    [ "$SEARCH" = 0 ] || { echo "PSEARCH or SEARCH_C1F / SEARCH_C3F, not both" >&2; exit 1; }
    [ "${MODES%% *}" = lab ] || { echo "PSEARCH needs MODES to start with lab (Pb and Ti take the lab leg's answer)" >&2; exit 1; }
    [ -x "${PYBIN}" ] || { echo "PSEARCH set but no abtem env python at ${PYBIN} (campaign/psearch_plan.py)" >&2; exit 1; }
    case "$PSEARCH_STAGES" in 2|3) ;; *) echo "PSEARCH_STAGES=${PSEARCH_STAGES}: 2 or 3" >&2; exit 1;; esac
    [ -z "$PSEARCH_RESUME" ] || [ -f "${PSEARCH_RESUME}/start.json" ] || { echo "PSEARCH_RESUME=${PSEARCH_RESUME}: no start.json there" >&2; exit 1; }
fi
psearch_job(){ # $1 name (the lab leg) $2 sim dir $3 bin $4 nl $5 dep $6 c1 $7 c3 -> "<job ids ...> <last select> <start file>"
    local name="$1" datadir="$2" bin="$3" nl="$4" dep="$5" c1="$6" c3="$7"
    local sd="search_af_${name}"; local sdd="${REPO_DIR}/${sd}"; mkdir -p "$sdd"; echo "${sd}" >>"${DIRS_FILE}"
    local c12 phi st rd f n cls grp da A S ids="" prev="" axes base fin
    read -r c12 phi < <(ABJ="${aj:--}" "${PYBIN}" -c 'import json, os; j = os.environ["ABJ"]; d = json.loads(j) if j not in ("", "-") else {}
print(d.get("C12", 0.0), d.get("phi12", 0.0))')
    "${PYBIN}" "${REPO_DIR}/campaign/psearch_plan.py" --name "$name" --nl "$nl" --c1="$c1" --c3="$c3" --c12="$c12" --phi12="$phi" \
        --ci95 "$PSEARCH_CI95" --seed "$PSEARCH_SEED" --n1 "$PSEARCH_N1" --n2 "$PSEARCH_N2" --stages "$PSEARCH_STAGES" \
        --out "$sdd" >&2 || return 1
    if [ -n "$PSEARCH_RESUME" ]; then          # the same draw, or nothing
        "${PYBIN}" -c 'import json, sys; a, b = (json.load(open(f + "/start.json"))["start"] for f in sys.argv[1:3])
ok = all(abs(a[k] - b[k]) <= 1e-6 * max(1.0, abs(b[k])) for k in ("c1_A", "c3_A", "c12_A", "phi12"))
sys.exit(0 if ok else "PSEARCH_RESUME: " + sys.argv[1] + " drew another start than this plan")' "$PSEARCH_RESUME" "$sdd" >&2 || return 1
    fi
    cls="$(res_class "$bin")"; grp="$(grp_for "$cls")"; da=$(dep_arg_for "$dep")
    local ex="NLAYERS=${nl},REGLAYER=0,PROBE_MODES=1,NITER=${SEARCH_NITER},SAVE_EVERY=${SEARCH_NITER},BETA_LSQ=${BETA_LSQ}${grp:+,GROUPING=${grp}}"
    for st in $(seq 1 "$PSEARCH_STAGES"); do
        if [ -n "$PSEARCH_RESUME" ] && [ -s "${PSEARCH_RESUME}/stage${st}_best.txt" ]; then    # done before: take its answer
            for f in "stage${st}_best.txt" "stage${st}_best.json" "stage${st}_trials.csv"; do
                [ -e "${PSEARCH_RESUME}/${f}" ] && cp "${PSEARCH_RESUME}/${f}" "${sdd}/${f}"
            done
            echo "  stage ${st}: resumed from ${PSEARCH_RESUME} ($(cat "${sdd}/stage${st}_best.txt"))" >&2
            prev=""; continue
        fi
        while IFS= read -r rd; do
            mkdir -p "${REPO_DIR}/${rd}/01"; echo "${rd}" >>"${DIRS_FILE}"
            if [ -d "${REPO_DIR}/${rd}/analysis" ]; then mv "${REPO_DIR}/${rd}/analysis" "${REPO_DIR}/${rd}/analysis.prev_$(date +%Y%m%d_%H%M%S)"; fi
            for f in "${INPUTS[@]}" aberrations.json probe_initial_true.mat; do ln -sf "${datadir}/01/${f}" "${REPO_DIR}/${rd}/01/${f}"; done
            rm -f "${REPO_DIR:?}/${rd:?}/01/probe_initial.mat" "${REPO_DIR:?}/${rd:?}/01/probe_initial.json"
        done < <(awk -F'\t' 'NR == 1 {for (j = 1; j <= NF; j++) if ($j == "recon_dir") c = j; next} {print $c}' "${sdd}/stage${st}.tsv")
        n=$(( $(wc -l <"${sdd}/stage${st}.tsv") - 1 ))
        axes=c1; [ "$st" = 2 ] && axes=a1
        base=""; [ "$st" -gt 1 ] && base="${sdd}/stage$((st - 1))_best.txt"
        fin=""; [ "$st" = "$PSEARCH_STAGES" ] && fin="--final"
        local dd_="$da"                                    # the sims (if any); after a stage run here, its select
        if [ -n "$prev" ]; then dd_="--dependency=afterany:${prev}"; fi
        A=$(sbatch --parsable --array="0-$((n - 1))" --job-name="af_ps${st}_${name}" --time="${SEARCH_RTIME}" \
            --mem="${SEARCH_MEM:-$(mem_for "$cls")}" ${dd_:+"$dd_"} \
            --output="${REPO_DIR}/recon_af_${name}_s${st}t%a_NL${nl}/slurm_%A_%a.out" --error="${REPO_DIR}/recon_af_${name}_s${st}t%a_NL${nl}/slurm_%A_%a.err" \
            --export=ALL,SEARCH_MANIFEST="${sdd}/stage${st}.tsv"${base:+,SEARCH_BASE_FILE=${base}},${ex} campaign/search_trial.sh); A="${A%%;*}"
        S=$(sbatch --parsable --job-name="af_psel${st}_${name}" --time=00:30:00 --mem=8G --cpus-per-task=1 --dependency="afterany:${A}" \
            --output="${sdd}/slurm_select${st}_%j.out" --error="${sdd}/slurm_select${st}_%j.err" \
            --wrap="'${PYBIN}' '${REPO_DIR}/campaign/select_grid.py' --manifest '${sdd}/stage${st}.tsv' --stage ${st} --axes ${axes}${base:+ --base '${base}'} ${fin} --out '${sdd}'")
        ids="${ids} ${A} ${S}"; prev="$S"
    done
    [ -n "$prev" ] || { echo "PSEARCH: every stage was resumed -- nothing to search, so nothing to start from" >&2; return 1; }
    echo "${ids# } ${sdd}/search_best.txt"
}
SELS=(); BESTS=()                  # per dose index: the lab leg's select job and its start file (Pb and Ti read them)
submit_legs(){ # $1 base name (<leg>_<mode><SFX>[_dose<D>]) $2 sim dir $3 dep $4 mode $5 dose index; appends to RIDS and line
    local base="$1" dd="$2" dep="$3" m="$4" di="$5" ids A S B
    if [ "$FIXED_TWIN" = 1 ]; then
        ids=$(recon_job "${base}${RTAG}_fixed" "$dd" "$bin" "$nl" "$dep"); RIDS+=($ids); line+=" ${m}:fixed=${ids}"
    fi
    if [ "$SEARCH" = 1 ]; then
        if [ "$m" = lab ]; then
            read -r A S B < <(search_job "${base}${RTAG}" "$dd" "$bin" "$nl" "$dep" "$c1" "$c3")
            SELS[$di]="$S"; BESTS[$di]="$B"; RIDS+=("$A" "$S"); line+=" ${m}:search=${A}>select=${S}"
        fi
        ids=$(kick_variants "${base}${RTAG}_rk" "$dd" "$bin" "$nl" "${dep:+afterok:${dep},}afterany:${SELS[$di]}" "$c1" "$c3" "${BESTS[$di]}")
        RIDS+=($ids); line+=" ${m}:rk=$(echo $ids | tr ' ' '+')"
    elif [ "$PSEARCH" = 1 ]; then
        if [ "$m" = lab ]; then
            local PS=() np
            read -r -a PS < <(psearch_job "${base}${RTAG}" "$dd" "$bin" "$nl" "$dep" "$c1" "$c3")
            np=${#PS[@]}
            [ "$np" -ge 3 ] || { echo "PSEARCH: planning the search for ${base}${RTAG} failed (above)" >&2; exit 1; }
            B="${PS[$((np - 1))]}"; SELS[$di]="${PS[$((np - 2))]}"; BESTS[$di]="$B"
            RIDS+=("${PS[@]:0:$((np - 1))}"); line+=" ${m}:psearch=$(echo "${PS[@]:0:$((np - 1))}" | tr ' ' '>')"
        fi
        ids=$(kick_variants "${base}${RTAG}_ps" "$dd" "$bin" "$nl" "${dep:+afterok:${dep},}afterany:${SELS[$di]}" "$c1" "$c3" "${BESTS[$di]}")
        RIDS+=($ids); line+=" ${m}:ps=$(echo $ids | tr ' ' '+')"
    elif [ "$KICK_ALL" = 1 ]; then
        ids=$(kick_variants "${base}${RTAG}" "$dd" "$bin" "$nl" "$dep" "$c1" "$c3"); RIDS+=($ids); line+=" ${m}=$(echo $ids | tr ' ' '+')"
    else
        ids=$(recon_job "${base}${RTAG}" "$dd" "$bin" "$nl" "$dep"); RIDS+=($ids); line+=" ${m}=${ids}"
    fi
}

RIDS=(); SIM_DIRS=()          # SIM_DIRS: this submission's own sim dirs, for CLEANDATA
ROWS="$ALPHAS"; [ -n "$LABELS" ] && ROWS="$LABELS"
for a in $ROWS; do
    if [ -n "$LABELS" ]; then       # by label: alpha and the full aberration JSON come from the row
        read -r alpha c3 c1 bin aj side win step detmax rnl rc5 < <(awk -F'\t' -v L="$a" '$1!~/^#/ && $1==L {
            for (i = 11; i <= 15; i++) if ($i == "") $i = "-"
            print $2"\t"$4"\t"$5"\t"$7"\t"$9"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$3}' "$TSV")
        ROW_C5="${rc5:-}"
        ROW_SIDE=""; ROW_WIN=""; ROW_STEP=""; ROW_DETMAX=""; ROW_NL=""   # [region] per-row geometry, "-" = driver default
        [ "${side:--}" != "-" ] && ROW_SIDE="$side"; [ "${win:--}" != "-" ] && ROW_WIN="$win"
        [ "${step:--}" != "-" ] && ROW_STEP="$step"; [ "${detmax:--}" != "-" ] && ROW_DETMAX="$detmax"
        [ "${rnl:--}" != "-" ] && ROW_NL="$rnl"
        [ -n "${bin:-}" ] || { echo "  ${a}: not in $TSV, skipping" >&2; continue; }
        leg="$a"
    else
        read -r c3 c1 bin < <(awk -F'\t' -v A="$a" '$1!~/^#/ && $2==A {print $4"\t"$5"\t"$7}' "$TSV")
        [ -n "${bin:-}" ] || { echo "  a${a}: not in $TSV, skipping" >&2; continue; }
        alpha="$a"; aj="-"; leg="a${a}"
    fi
    nl=$(nl_full "$alpha")
    # col 15 (nl_force) of a LABELS row forces its NL (a slice-count test that keeps the row's own label, sim and
    # kernels together); the NL env var still wins. Col 8 is not used: older tsvs carry stale values there.
    [ -z "${NL:-}" ] && [ -n "${ROW_NL:-}" ] && nl="$ROW_NL"
    line="$(printf '%-14s bin=%s NL=%-2s ' "${leg}${SFX}" "$bin" "$nl")"
    [ -n "${ROW_SIDE:-}" ] && line+="box=${ROW_SIDE} win=${ROW_WIN:-$WIN} step=${ROW_STEP:-$STEP} det=${ROW_DETMAX:-200} class=$(res_class "$bin") "
    for m in $MODES; do                       # MODES="Pb Ti" re-does only the PSF kernels
        D="${REPO_DIR}/sim_out_af_${leg}_${m}${SFX}"
        if [ -n "$DTAG" ]; then                # a smaller camera on the sim (DETECTOR above), as DOSES does with noise
            S=""
            if [ ! -e "${D}/01/data_dp.hdf5" ]; then   # not there: simulate it in this submission, then rebin it
                [ "${RECON_ONLY:-0}" = "1" ] && { echo "  ${leg} ${m}: ${D}/01/data_dp.hdf5 missing and RECON_ONLY=1" >&2; exit 1; }
                S=$(sim_job "$D" "$alpha" "$bin" "$c3" "$c1" "$m" "$aj"); SIM_DIRS+=("$D"); line+=" ${m}:sim=${S}"
            fi
            DN="${D}${DTAG}"
            if [ "${RECON_ONLY:-0}" = "1" ] && [ -e "${DN}/01/rebin_detector.json" ]; then   # written last: the copy is whole
                N=""; line+=" ${m}${DTAG}:reuse"      # rewriting it could pull the data from under a running leg
            else
                N=$(rebin_job "$D" "$DN" "$S"); line+=" ${m}${DTAG}:rebin=${N}"
            fi
            SIM_DIRS+=("$DN")                      # the copy only: a kept sim stays (CLEANDATA)
            submit_legs "${leg}_${m}${SFX}${DTAG}${UTAG}" "$DN" "$N" "$m" 0
            continue
        fi
        if [ -n "$DOSES" ]; then               # step 2: Poisson copies of the noiseless sim
            S=""
            if [ ! -e "${D}/01/data_dp.hdf5" ]; then   # not there yet: simulate it in this submission, then noise it
                [ "${RECON_ONLY:-0}" = "1" ] && { echo "  ${leg} ${m}: ${D}/01/data_dp.hdf5 missing and RECON_ONLY=1" >&2; exit 1; }
                S=$(sim_job "$D" "$alpha" "$bin" "$c3" "$c1" "$m" "$aj"); SIM_DIRS+=("$D"); line+=" ${m}:sim=${S}"
            fi
            case "$m" in lab) so=0;; Pb) so=1;; *) so=2;; esac
            di=0
            for dose in $DOSES; do
                DN="${D}_dose${dose}"
                if [ "${RECON_ONLY:-0}" = "1" ] && [ -e "${DN}/01/data_dp.hdf5" ]; then
                    N=""; line+=" ${m}@${dose}:reuse"   # RECON_ONLY: the Poisson copy is there (CLEANDATA=0 kept it) --
                    # re-noising would rewrite data_dp.hdf5 under any job still reading it
                else
                    N=$(noise_job "$D" "$DN" "$dose" $(( DOSE_SEED + so )) "$S"); line+=" ${m}@${dose}:noise=${N}"
                fi
                SIM_DIRS+=("$DN")   # $DN, not $D: the noiseless
                # original is the SOURCE for every other dose in this submission
                submit_legs "${leg}_${m}${SFX}_dose${dose}" "$DN" "$N" "$m" "$di"; di=$((di + 1))
            done
            continue
        fi
        if [ "${RECON_ONLY:-0}" = "1" ]; then  # reuse existing sims; recons run immediately
            [ -e "${D}/01/data_dp.hdf5" ] || { echo "  a${a} ${m}: ${D}/01 missing, skip" >&2; continue; }
            S=""
        else
            # run_sim.slurm refuses to overwrite a finished sim (exit 1, zero seconds), and every recon
            # waiting on it is then DependencyNeverSatisfied. On 2026-09-22 the whole 70 A retry died
            # this way after 90 min in the queue: the sims already existed, only BETA_LSQ had changed.
            # Catch it here, on the login node, instead of on the GPU.
            if [ -e "${D}/01/data_dp.hdf5" ] && [ "${OVERWRITE:-0}" != "1" ]; then
                echo "ERROR: ${D}/01/data_dp.hdf5 already exists -- the sim would be refused and every" >&2
                echo "       recon left stuck on an unsatisfiable dependency." >&2
                echo "       Re-running only the solver? RECON_ONLY=1.  Re-simulating? OVERWRITE=1." >&2
                exit 1
            fi
            S=$(sim_job "$D" "$alpha" "$bin" "$c3" "$c1" "$m" "$aj")
        fi
        SIM_DIRS+=("$D")
        submit_legs "${leg}_${m}${SFX}" "$D" "$S" "$m" 0
        if [ "$KICK_ALL" != 1 ] && [ "$SEARCH" != 1 ] && [ "$m" = lab ] && [ -n "$LABELS" ] && [[ " ${KICK_LABELS} " == *" ${leg} "* ]]; then
            K=$(kick_variants "${leg}_kick${SFX}${RTAG}" "$D" "$bin" "$nl" "$S" "$c1" "$c3")
            RIDS+=($K); line+=" kick=$(echo $K | tr ' ' '+')"
        fi
    done
    echo "$line"
done
[ ${#RIDS[@]} -gt 0 ] || { echo "nothing submitted" >&2; exit 1; }
DEP=$(IFS=:; echo "${RIDS[*]}")
# CLEANDATA=1: once this submission's results are tarred, delete ITS raw 4D data (data_dp/position
# .hdf5, 0.8 GB a leg at BIN 4 and 3.2 GB at BIN 2). The trade-off is that a later RECON_ONLY re-fit
# needs a re-simulation, which is minutes at BIN 4. Same flag and same mechanism as run_campaign.sh.
# Only this submission's own sim dirs are touched, never another run's.
CLEAN=""
if [ "${CLEANDATA:-0}" = "1" ] && [ "${RECON_ONLY:-0}" = "1" ]; then
    echo "WARNING: CLEANDATA=1 with RECON_ONLY=1 deletes raw data this submission did not create," >&2
    echo "         so a later re-fit of these legs needs a re-simulation. Ctrl-C now if unintended." >&2
    sleep 5
fi
if [ "${CLEANDATA:-0}" = "1" ]; then
    for leg_dir in "${SIM_DIRS[@]}"; do
        CLEAN="${CLEAN} && find '${leg_dir}' \\( -name data_dp.hdf5 -o -name data_position.hdf5 \\) -delete"
    done
fi
PJ=$(sbatch --parsable --job-name="af_pack" --time=00:30:00 --mem=8G --dependency="afterany:${DEP}" \
    --output="logs/af_pack_%j.out" --error="logs/af_pack_%j.err" --export=ALL,PACK_DIRS_FILE="${DIRS_FILE}" \
    --wrap="bash '${REPO_DIR}/campaign/pack_results.sh' af '${PACK}' '${TSV}'${CLEAN}")
echo; echo "pack ${PJ} -> ${PACK}  (recon_af_* .h5 + logs)$([ "${CLEANDATA:-0}" = 1 ] && echo " then deletes the raw 4D data of this submission" || true)"
echo "scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:${PACK}' ~/Desktop/"
echo "now check the group root stayed clean:  ls /springbrook/share/physics/"
echo "then analyse ON BLYTHE (kernels + atomfind + summary.csv, geometry read from each leg):"
echo "  LABELS=\"<labels>\" DEP=${PJ} GT_REGION=<box side, region rows> | GT=<gt dir> bash campaign/run_analysis.sh"
