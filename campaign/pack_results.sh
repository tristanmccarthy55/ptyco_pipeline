#!/usr/bin/env bash
# [campaign] Gather ONE campaign's results into a single tgz so the whole sweep comes down
# in one scp. Submitted by run_campaign.sh as an afterany-dependent job (runs once every
# recon has finished, success or not, so partial sweeps still pack). Grabs the recon .h5
# (object + recovered probe + all params live inside), the true probes, sim_meta, the
# per-sim aberrations.json, and the sweep .tsv (provenance). Also the small sidecars written by
# ptycho/run_synthetic_recon_ML.m (*_error_trace.csv, *_layer_stats.csv) and, for trials whose
# probe was written in-job by sim/make_probe.py, that probe + its .json (real files only: the
# known-probe legs symlink theirs, and a symlink in a tarball is a dangling link locally).
# PACK_H5=0 skips the *_recons.h5 (a C1-search grid at BIN=2 is 25 x 240 MB; the objective lives in
# the sidecars).
# PACK_DIRS_FILE=<file> restricts the recon dirs to the ones listed (one per line, relative to the
# repo) -- the drivers write it so a submission packs ITS OWN recons. Without it every
# recon_<campaign>_* dir is packed, so two submissions running side by side would each ship the
# other's h5s as soon as they exist.
#   campaign/pack_results.sh <campaign> <out.tgz> [sweep.tsv]
set -uo pipefail
CAMP="${1:?campaign name}"; OUT="${2:?output tgz}"; TSV="${3:-}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"
list="$(mktemp)"
RDIRS=()
if [ -n "${PACK_DIRS_FILE:-}" ]; then
    [ -f "${PACK_DIRS_FILE}" ] || { echo "pack: PACK_DIRS_FILE ${PACK_DIRS_FILE} missing" >&2; exit 1; }
    while IFS= read -r d; do [ -n "$d" ] && [ -d "$d" ] && RDIRS+=("$d"); done <"${PACK_DIRS_FILE}"
    echo "pack: ${#RDIRS[@]} recon dirs from ${PACK_DIRS_FILE}"
else
    for d in "recon_${CAMP}_"*; do [ -d "$d" ] && RDIRS+=("$d"); done
fi
if [ ${#RDIRS[@]} -gt 0 ]; then
    if [ "${PACK_H5:-1}" = "1" ]; then
        find "${RDIRS[@]}" -path '*/analysis/*_recons.h5' 2>/dev/null >>"$list" || true
    else
        echo "pack: PACK_H5=0 -- *_recons.h5 NOT packed (sidecars, logs and probes only)"
    fi
    find "${RDIRS[@]}" -path '*/analysis/*_error_trace.csv' 2>/dev/null >>"$list" || true
    find "${RDIRS[@]}" -path '*/analysis/*_layer_stats.csv' 2>/dev/null >>"$list" || true
    find "${RDIRS[@]}" -path '*/01/probe_initial.json' -type f 2>/dev/null >>"$list" || true
    find "${RDIRS[@]}" -path '*/01/probe_initial.mat'  -type f 2>/dev/null >>"$list" || true
    find "${RDIRS[@]}" -name 'slurm_*.out'                    2>/dev/null >>"$list" || true
fi
# sim provenance (small). Real files only: a Poisson copy links these back to its noiseless sim,
# which is packed itself.
find "sim_out_${CAMP}_"*    -name 'probe_initial_true.mat' -type f 2>/dev/null >>"$list" || true
find "sim_out_${CAMP}_"*    -name 'sim_meta.mat'           -type f 2>/dev/null >>"$list" || true
find "sim_out_${CAMP}_"*    -name 'aberrations.json'       -type f 2>/dev/null >>"$list" || true
find "sim_out_${CAMP}_"*    -name 'poisson_noise.json'     -type f 2>/dev/null >>"$list" || true
[ -n "$TSV" ] && [ -f "$TSV" ] && echo "$TSV" >>"$list"
n=$(wc -l <"$list")
if [ "$n" -eq 0 ]; then echo "pack: nothing found for campaign '${CAMP}'" >&2; rm -f "$list"; exit 1; fi
tar czf "$OUT" -T "$list"
echo "packed ${n} files -> ${OUT}"; du -h "$OUT"
rm -f "$list"
