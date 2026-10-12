# Phase 1 (2026-10-12): 80 mrad, ARM tableau, noiseless; per point 3 sims, 3 true-probe-fixed legs (_fixed, 1 mode, 200 it)
# and 3 true-start probe-update legs (KICK=1.0, modes from mode_spectrum.txt: 99 % of the real-gun state + 1, cap 8; 250 it).
# Run from the repo root on Blythe. Design, cost and what to watch: aberration_experiment/mono300/PHASE1.md.
AIDS=""
for X in "FOCAL_SPREAD=7.29 8 4;2" "FOCAL_SPREAD=5 7 4;2" "FOCAL_SPREAD=3.5 5 4;2" "FOCAL_SPREAD=2.5 4 8;4" \
         "FOCAL_SPREAD=1.5 4 8;4" "FOCAL_SPREAD=0.75 3 8;4" "SOURCE_FWHM=0.05 4 8;4" "SOURCE_FWHM=0.1 8 4;2" \
         "SOURCE_FWHM=0.2 8 4;2" "SOURCE_FWHM=0.4 8 4;2"; do
  read -r P M G <<<"$X"
  OUT=$(env "$P" FIXED_TWIN=1 KICK_ALL=1 KICK=1.0 KICK_MODES="$M" KICK_NITER=250 KICK_GROUPING="$G" \
        LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 CLEANDATA=1 \
        bash campaign/run_thin_atomfind.sh) || { echo "$OUT"; echo "STOPPED at $P"; break; }
  SFX=$(echo "$OUT" | sed -n 's/.*; suffix \(_[^ ;]*\)$/\1/p'); PJ=$(echo "$OUT" | awk '/^pack /{print $2}')
  [ -n "$SFX" ] && [ -n "$PJ" ] || { echo "$OUT"; echo "STOPPED: no suffix / pack id for $P"; break; }
  echo "$OUT" | sed -n '2p;/^pack /p'; echo "  ^ $P: $M modes, grouping $G"
  for S in "$SFX" "${SFX}_fixed"; do
    A=$(LABELS=armf_a080 SUFFIX="$S" DEP="$PJ" GT_REGION=210 bash campaign/run_analysis.sh | head -1)
    echo "$A"; AIDS="$AIDS $(echo "$A" | awk '{print $2}')"
  done
done
echo "ANALYSES:$AIDS"
