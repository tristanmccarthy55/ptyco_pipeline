#!/usr/bin/env bash
# [campaign] Gather ONE campaign's results into a single tgz so the whole sweep comes down
# in one scp. Submitted by run_campaign.sh as an afterany-dependent job (runs once every
# recon has finished, success or not, so partial sweeps still pack). Grabs the recon .h5
# (object + recovered probe + all params live inside), the true probes, sim_meta, the
# per-sim aberrations.json, and the sweep .tsv (provenance).
#   campaign/pack_results.sh <campaign> <out.tgz> [sweep.tsv]
set -uo pipefail
CAMP="${1:?campaign name}"; OUT="${2:?output tgz}"; TSV="${3:-}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO_DIR}"
list="$(mktemp)"
find "recon_${CAMP}_"*      -path '*/analysis/*_recons.h5' 2>/dev/null >>"$list" || true
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
