#!/usr/bin/env bash
# [campaign] Gather ONE campaign's results into a single tgz so the whole sweep comes down
# in one scp. Submitted by run_campaign.sh as an afterany-dependent job (runs once every
# recon has finished, success or not, so partial sweeps still pack). Grabs the recon .h5
# (object + recovered probe + all params live inside), the true probes, sim_meta, the
# per-sim aberrations.json, and the sweep .tsv (provenance). Also the small sidecars written by
# ptycho/run_synthetic_recon_ML.m (*_error_trace.csv, *_layer_stats.csv) and, for trials whose
# probe was written in-job by sim/make_probe.py, that probe + its .json (real files only: the
# known-probe legs symlink theirs, and a symlink in a tarball is a dangling link locally).
# PACK_H5=0 skips the *_recons.h5 (a C1-search grid at BIN=2 is 22 x 240 MB; the objective lives in
# the sidecars).
#   campaign/pack_results.sh <campaign> <out.tgz> [sweep.tsv]
set -uo pipefail
CAMP="${1:?campaign name}"; OUT="${2:?output tgz}"; TSV="${3:-}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"
list="$(mktemp)"
if [ "${PACK_H5:-1}" = "1" ]; then
    find "recon_${CAMP}_"*  -path '*/analysis/*_recons.h5' 2>/dev/null >>"$list" || true
else
    echo "pack: PACK_H5=0 -- *_recons.h5 NOT packed (sidecars, logs and probes only)"
fi
find "recon_${CAMP}_"*      -path '*/analysis/*_error_trace.csv' 2>/dev/null >>"$list" || true
find "recon_${CAMP}_"*      -path '*/analysis/*_layer_stats.csv' 2>/dev/null >>"$list" || true
find "recon_${CAMP}_"*      -path '*/01/probe_initial.json' -type f 2>/dev/null >>"$list" || true
find "recon_${CAMP}_"*      -path '*/01/probe_initial.mat'  -type f 2>/dev/null >>"$list" || true
find "recon_${CAMP}_"*      -name 'slurm_*.out'            2>/dev/null >>"$list" || true
find "sim_out_${CAMP}_"*    -name 'probe_initial_true.mat' 2>/dev/null >>"$list" || true
find "sim_out_${CAMP}_"*    -name 'sim_meta.mat'           2>/dev/null >>"$list" || true
find "sim_out_${CAMP}_"*    -name 'aberrations.json'       2>/dev/null >>"$list" || true
[ -n "$TSV" ] && [ -f "$TSV" ] && echo "$TSV" >>"$list"
n=$(wc -l <"$list")
if [ "$n" -eq 0 ]; then echo "pack: nothing found for campaign '${CAMP}'" >&2; rm -f "$list"; exit 1; fi
tar czf "$OUT" -T "$list"
echo "packed ${n} files -> ${OUT}"; du -h "$OUT"
rm -f "$list"
