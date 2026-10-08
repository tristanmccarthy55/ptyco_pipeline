# Does the atom finder's score show depth? (2026-10-08)

Question (user): on the phonon + partial-coherence runs the finder "snipes atoms that aren't there". Everything here runs
OUTSIDE the atomfind package (read-only use of its functions); the package is unchanged.

## 1. Where the flattering numbers come from
- `validate.match_found_to_gt`: 0.6 A in plane, **2.0 A in depth, any species**. Along a Pb column (Pb every 3.9 A) any
  depth is within 1.95 A of a Pb, and depth rms is over matched atoms only (capped ~1.15 A).
- Guided (lattice-prior) detections are NOT the cause: 2-8 % of found atoms on every run, oxygen only.
- The finder already prints species confusion + health warnings (logs/atomfind_*.log); the page does not show them:
  known_a080 1.1 %, win_070 0.4 %, long_m6 7.3 %, fixed_1e7 29.0 %, fixed_1e8 21.6 %, hail_1e7 34.7 %.
- The published depth sections take their colour scale over ALL layers; on full-physics runs the vacuum layers hold
  bright junk, so the slab reads black and the rings look as if they sit on nothing (ring_support.png, left panels).

## 2. Strict score with a half-period control (local re-run reproduces the shipped atoms exactly)
strict = right element AND within 0.5 A in depth. control = the same atoms moved 1.95 A (half the Pb period).
Pb/Ti/O recall %:

| run | strict | strict, moved 1.95 A | slices doubled by interpolation, strict |
|---|---|---|---|
| known_a080 | 90/58/59 | 0/0/1 | 78/54/52 |
| win_070 | 95/56/55 | 0/0/2 | 80/49/46 |
| long_m6 | 66/47/66 | 0/0/1 | 64/45/49 |
| fixed_1e7 | 23/0/5 | 21/0/7 | 16/4/4 |
| fixed_1e8 | 23/3/5 | 20/0/6 | 11/8/2 |
| hail_1e7 | 4/0/2 | 5/0/2 | 3/0/3 |

Full-physics runs sit at chance (~26 % of Pb fall within 0.5 A of a Pb by luck: 1 A window / 3.9 A period).
Doubling the slices by interpolation (volume and kernel phase, same convention) adds detections (loose recall up,
precision down 0.97 -> 0.89 on known_a080) and lowers every strict score; no gain on any run.

## 3. Blind per-column depth evidence (scripts/colshift.py)
Per column: NNLS fit of the Pb kernel at the found depths vs the same atoms moved +-half their own spacing;
e = (chi2 shifted - chi2 found) / column depth-structure energy. Null: the column with the found atoms' depth
modulation replaced by its depth mean, put through the finder's own CLEAN + refine + quality cut, scored the same.
Median e (null median, null p95), Pb columns:

| run | Pb columns e | null | B-site e | null | O columns e | null |
|---|---|---|---|---|---|---|
| known_a080 | 0.95 | 0.07 (0.15) | 0.43 | 0.17 (0.58) | 0.88 | 0.34 (0.68) |
| win_070 | 0.94 | 0.06 (0.14) | 0.44 | 0.26 (0.63) | 0.78 | 0.44 (0.73) |
| long_m6 | 0.84 | 0.09 (0.20) | 0.30 | 0.20 (0.34) | 0.28 | 0.17 (0.43) |
| fixed_1e7 | 0.05 | 0.05 (0.58) | 0.13 | 0.07 | 0.08 | 0.07 |
| fixed_1e8 | 0.07 | 0.12 (0.24) | -0.11 | 0.08 | 0.10 | 0.06 |
| hail_1e7 | 0.04 | 0.06 (0.20) | - | - | 0.02 | 0.04 |

Pb columns separate cleanly with no ground truth; B-site columns are too weak for a per-atom mark.

Scripts are the working copies (quick, not tidy): run from the repo root with `S=<this folder>/scripts`, e.g.
`S=$PWD/aberration_experiment/results/2026-W41/atomfind_honesty/scripts FS=1,2 ~/hyperspy-bundle/bin/python $S/refind.py known_a080`.
