#!/usr/bin/env bash
# [driver] Dry run: every launch line the ARM results page records (aberration_experiment/page/arm_results.json, its "cmd"
# and "physics_cmd" fields) through a REFERENCE driver (a git revision) and the working-tree driver, with a stub sbatch
# (an exported function: nothing is ever submitted); the sbatch calls, the driver's output and the links it makes are
# compared line for line. A driver change that keeps old runs' names and settings shows every line "same".
#   bash campaign/test_driver_legacy.sh            # against HEAD
#   bash campaign/test_driver_legacy.sh ac69bb2    # against an older revision
# RECON_ONLY lines are given the old names' data, so a renamed sim shows up as a difference, not a skip.
# ONE DELIBERATE RENAME is mapped onto the reference before the comparison: per-species phonon sims gained v2 on
# 2026-10-30 (_ph<N> -> _ph<N>v2; the abTEM sigma fix, sim/simulate_4dstem.py a8b9264), so for a PER_SPECIES=1 launch
# the reference's _ph<N> is read as _ph<N>v2 and RECON_ONLY roots hold both names' data. Anything else must be "same".
set -uo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO_DIR"
REF="${1:-HEAD}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
git show "${REF}:campaign/run_thin_atomfind.sh" >"$TMP/driver_ref.sh" || { echo "no driver at ${REF}" >&2; exit 1; }
mkdir -p "$TMP/env/bin"; printf '#!/bin/sh\nexit 0\n' >"$TMP/env/bin/python"; chmod +x "$TMP/env/bin/python"
sbatch(){ local n; n=$(( $(cat "$STUB_N" 2>/dev/null || echo 1000) + 1 )); echo "$n" >"$STUB_N"
          printf '%s | PROBE_SCALE=%s\n' "$*" "${PROBE_SCALE:-}" >>"$STUB_LOG"; echo "$n"; }
export -f sbatch
DVARS=(); while IFS= read -r v; do DVARS+=("$v"); done < <(grep -ohE '\$\{[A-Z][A-Z0-9_]*:[-+=]' campaign/run_thin_atomfind.sh \
    "$TMP/driver_ref.sh" | sed -E 's/^\$\{([A-Z0-9_]+).*/\1/' | sort -u)
PY="${PYTHON:-python3}"
"$PY" - "$REPO_DIR/aberration_experiment/page/arm_results.json" >"$TMP/cmds.txt" <<'EOF'
import json, sys
seen = []
def walk(o):
    if isinstance(o, dict):
        for k, v in o.items():
            if k in ("cmd", "physics_cmd") and isinstance(v, str) and "run_thin_atomfind" in v: yield v
            else: yield from walk(v)
    elif isinstance(o, list):
        for v in o: yield from walk(v)
for c in walk(json.load(open(sys.argv[1]))):
    if c not in seen: seen.append(c)
print("\n".join(seen))
EOF
[ -s "$TMP/cmds.txt" ] || { echo "no launch lines read from the manifest (python: ${PY})"; exit 1; }   # never pass on nothing
echo "reference ${REF} vs working tree, $(grep -c . "$TMP/cmds.txt") recorded launch lines"
i=0; bad=0
while IFS= read -r cmd; do
  i=$((i + 1))
  for which in ref new; do
    R="$TMP/$which/$i"; mkdir -p "$R/campaign"
    if [ $which = ref ]; then cp "$TMP/driver_ref.sh" "$R/campaign/run_thin_atomfind.sh"; else cp campaign/run_thin_atomfind.sh "$R/campaign/"; fi
    cp campaign/ceos_sweep.tsv "$R/campaign/"
    if [[ "$cmd" == *RECON_ONLY=1* ]]; then
      for s in sim_out_af_armf_a080_{lab,Pb,Ti}_ph16{,v2}_coh{,_dose1e7,_dose1e8}; do mkdir -p "$R/$s/01"; : >"$R/$s/01/data_dp.hdf5"; done
    fi
    c="${cmd//bash campaign\/run_thin_atomfind.sh/bash \"$R/campaign/run_thin_atomfind.sh\"}"
    u=(-u BASH_ENV -u ENV); for v in "${DVARS[@]}"; do u+=(-u "$v"); done
    env "${u[@]}" STUB_LOG="$R/sbatch.log" STUB_N="$R/n" CONDA_ENV="$TMP/env" SHARE= bash -c "$c" >"$R/out.txt" 2>&1
    echo "exit $?" >>"$R/out.txt"
    { cat "$R/out.txt" "$R/sbatch.log" 2>/dev/null
      (cd "$R" && find . -path ./campaign -prune -o -type l -print | sort | while read -r l; do echo "$l -> $(readlink "$l")"; done)
    } | sed -E "s#$R#ROOT#g; s/[0-9]{8}_[0-9]{6}/TS/g" >"$R/norm.txt"
    if [ $which = ref ] && [[ "$cmd" == *PER_SPECIES=1* ]] && [[ "$cmd" == *PHONONS=* ]]; then   # the v2 rename (header)
      sed -E 's/_ph([0-9]+)([^0-9v])/_ph\1v2\2/g; s/_ph([0-9]+)$/_ph\1v2/' "$R/norm.txt" >"$R/norm2.txt" && mv "$R/norm2.txt" "$R/norm.txt"
    fi
  done
  n=$(grep -c . "$TMP/new/$i/sbatch.log" 2>/dev/null || echo 0)
  if diff -q "$TMP/ref/$i/norm.txt" "$TMP/new/$i/norm.txt" >/dev/null; then
    echo "same  #$i  (${n} sbatch calls, $(tail -1 "$TMP/new/$i/out.txt"))  $(echo "$cmd" | cut -c1-80)"
  else
    echo "DIFF  #$i  $(echo "$cmd" | cut -c1-80)"; diff "$TMP/ref/$i/norm.txt" "$TMP/new/$i/norm.txt" | head -8; bad=1
  fi
done <"$TMP/cmds.txt"
[ "$bad" = 0 ] && echo "all launch lines unchanged" || { echo "SOME LAUNCH LINES CHANGED"; exit 1; }
