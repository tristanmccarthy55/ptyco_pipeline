# The Coherence Wall page — how to fill it in and republish

PUBLISHED: https://claude.ai/artifact/B7ep8BpAMVdKvrnXywACq3
- version 5, 2026-10-21, for the meeting: new section 06 "The Spectra Mono point" -- A (ceiling with modes: good at
  src 0.1, okay at 0.2), B (cooling: fixed probe off at 100 K and 20 K, modes okay with Pb 51 / 70 %, Ti at guessing),
  C (searched start: the first 2-stage search failed on its order -- C1 32.5 A off -- the 3-stage rerun pending), the
  search table (`search_table`), the phase-2 strip `fig_strip_p2.png`; the 7.29 A `_mix` leg shown as "earlier solver"
  until the `_mixF` reruns land (`score_fallbacks`); sections renumbered (instruments 07, discussion 08).
- version 4, 2026-10-16, for the meeting: text cut, a "Parameters and method" table (what the two knobs are), a steps table for what was done, no guessing lines on the wall chart, plain label for the real-source band.
- version 3, 2026-10-16, for the meeting: section 01 rewritten as what we did / what we found / what is next; no new runs.
- version 2, 2026-10-16: step 2 in. Joint points (the losses compound: 2.5 + 0.2 A off though each part passed alone),
  the 7.29 A probe-mode leg (off -> okay, state stays mixed), the 0.2 A source mode leg as a volume without atoms, the
  engine's mode-1 step rule as a caveat on every mode run. Pending: the fs 5 / src 0.2 / src 0.4 reruns (`_mix16`, all
  three legs at 16;8), provisional scores (`_mix_kern-..._fixed`), phase 2 (not designed).
- version 1, 2026-10-15: step 1 in full; step 2's joint points and probe-mode legs pending. Local copy of each published version: `PtycoShelves/artifacts/wcm_2026-10-12_coherence_wall.html`.

The 300 kV monochromated-microscope results (phase 1, `aberration_experiment/mono300/PHASE1.md`), built the same way as
the ARM page (`ARM_PAGE.md`): generated from the run outputs, never hand-edited. It imports `analysis/make_arm_page.py`
for the atom scores, verdicts, depth sections, tab viewer and its template's CSS, so the two pages stay one system.

| file | what it is |
|---|---|
| `aberration_experiment/page/mono_results.json` | the run list: tag, leg, group (section), focal spread / source, pending text; the instrument table with sources |
| `analysis/make_mono_page.py` | builds every figure into `aberration_experiment/figs/mono_page/` and writes the page |
| `aberration_experiment/page/mono_2026-10-12.html` | the template: prose, `{{...}}` slots, `FILL` blocks for step 2 |
| `aberration_experiment/page/built/mono_page_built.html` | the page to publish (gitignored) |

## When step 2 lands

1. **Pull** into `~/Desktop/mono_phase1` with the pull blocks in PHASE1.md; extract every tarball (analysis into
   `analysis/`, atomfind_results into `results/`).
2. **See what landed:** `~/hyperspy-bundle/bin/python analysis/make_mono_page.py --list`.
3. **Rebuild:** `~/hyperspy-bundle/bin/python analysis/make_mono_page.py` (about 15 s; prints the size, stay under 16 MB).
   Landed runs fill the joint grid (§05), the mode table, mode-power figure and strip (§06), the diamonds on the wall
   figure, and the status ledger, with no edits.
4. **Look first** (the user's rule): the strips and tabs of the new runs, then `figs/mono_page/fig_modes.png`. If the
   recovered mode powers collapsed onto mode 1, the reconstruction left the mixed state: say so before reading the atoms.
5. **Write** the `FILL` blocks (`grep -n FILL aberration_experiment/page/mono_2026-10-12.html`): headline (§01), joint (§05),
   modes (§06), discussion (§08); replace the matching `call pending` box. Short, plain, numbers from the tables.
6. **Rebuild and republish** with the Artifact tool, `url` = the PUBLISHED url above (a conversation that did not publish
   it reads it first), no `icon`. Copy the built page over the local copy in `PtycoShelves/artifacts/`.
7. **Commit** the template, run list, script and the `figs/mono_page/fig_*.png` print copies (not `built/`, not `row_*.png`).

## Adding a run

A mode point lists its tags in priority order with a matching `legs` list: the `_mix16` rerun (own kernels), then the
`_mix_kern-<control>` analysis (provisional: atoms read with the fixed-probe control's kernels, labelled so on the page),
then the original `_mix` leg (a reconstruction without atoms, or a crash). The page shows whichever landed first in that
order, and its state (in / provisional / volume only / crashed / pending).

Add a run to `runs` in the run list with its `group` (`focal`, `source`, `joint`, `modes`, or `ref`), `fs` and `src`; a
`modes` run also takes `control` (its fixed-probe twin) and `modes`. Phase 2 (phonons) will want a new group and section.

## Numbers the page computes

- **strongest mode**: the largest eigenvalue share of the simulator's own 16-configuration state (seed 7), from
  `mono300/mode_spectrum.py`'s own functions (its print loop is skipped).
- **vacuum / slab**: rms phase per slice in the two 4 Å vacuum bands over that in the slab (coherent reference 0.31).
- **instrument ranges**: `simulate_4dstem`'s focal-spread formula on the run list's energy spreads and Cc values.
- The wall's edges (last good, first off) and the energy spread for good are read from the verdicts, not typed.
