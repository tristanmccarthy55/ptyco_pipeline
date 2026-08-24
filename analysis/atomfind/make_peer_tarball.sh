#!/usr/bin/env bash
# @file make_peer_tarball.sh
# @brief Build the data tarball that accompanies the standalone atomfind repo.
#
# The code is small and lives in git; the inputs are ~95 MB and do not, so they ship as one
# checksummed tarball (a GitHub release asset). Contents are exactly what PEER.md lists.
#
#   bash analysis/atomfind/make_peer_tarball.sh [OUTDIR]     # default ~/
set -euo pipefail

PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$PKG/../.." && pwd)"
OUT="${1:-$HOME}"
STAGE="$(mktemp -d)"; trap 'rm -rf "$STAGE"' EXIT
NAME="atomfind_data_v1"

# name -> where to find it. Sizes/roles are documented in PEER.md.
copy () {                       # copy <basename> <source-path> <what-it-is>
  [ -f "$2" ] || { echo "MISSING: $3 -- expected at $2"; exit 1; }
  cp "$2" "$STAGE/$NAME/$1"
  printf '  %-42s %8s  %s\n' "$1" "$(du -h "$2" | cut -f1)" "$3"
}

mkdir -p "$STAGE/$NAME"
echo "Staging $NAME ..."
copy NL70_new_vol.npy        "$HOME/Desktop/NL70_new_vol.npy"       "reconstructed phase volume (70,404,404)"
copy psf_Pb_NL70_vol.npy     "$HOME/Desktop/psf_Pb_NL70_vol.npy"    "measured single-Pb kernel"
copy psf_Ti_NL70_vol.npy     "$HOME/Desktop/psf_Ti_NL70_vol.npy"    "measured single-Ti kernel"
copy gt_prepared.npz         "$PKG/data/gt_prepared.npz"            "reference structure, beam frame"
copy PTO6_STO6_18_18_labyrinthPoscar.vasp \
                             "$REPO/sim/PTO6_STO6_18_18_labyrinthPoscar.vasp" \
                                                                    "reference structure, raw"

( cd "$STAGE/$NAME" && shasum -a 256 * > SHA256SUMS )
cp "$PKG/PEER.md" "$STAGE/$NAME/PEER.md"

mkdir -p "$OUT"
tar czf "$OUT/$NAME.tar.gz" -C "$STAGE" "$NAME"
echo
echo "Wrote $OUT/$NAME.tar.gz  ($(du -h "$OUT/$NAME.tar.gz" | cut -f1))"
echo "Verify after transfer:  tar xzf $NAME.tar.gz && cd $NAME && shasum -a 256 -c SHA256SUMS"
