# The Coherence Wall page — how to fill it in and republish

PUBLISHED: https://claude.ai/artifact/B7ep8BpAMVdKvrnXywACq3
- version 14, 2026-10-21: installs sorted by gun (Mono, UltiMono, monochromated-unstated, X-CFEG/unstated; ARM for
  comparison; Titan dropped, the user: no better than the ARM). The UltiMono's 25-55 meV was a 30-60 kV figure: at 300 kV
  the quoted numbers are 110 meV (Stanford) and < 0.2 eV (Michigan), so its bar is now 0.11-0.2 eV and the page says the
  Mono and UltiMono are about the same at 300 kV; the target is "a monochromated Spectra 300 at ~0.1 eV". Europe: >= 5
  monochromated Spectra 300s (Cambridge Mono; Imperial, Barcelona, Jülich UltiMono; Stuttgart unstated).
- version 13, 2026-10-21: every outside number linked. mono_results.json gains "sources" (cite, url, what it
  supports; rendered as section C by `sources_list`, inline by `src_link`, per instrument via "refs") and "installs"
  (model, region uk/eu/row, site, url, note; `installs_table` in §06, "How common each is": lower bounds from public
  facility pages). Corrections from the sources: the X-FEG/Mono's specification is 0.2-0.3 eV (DS0297; < 0.2 eV is its
  tuning range), the X-FEG's own spread ~0.8 eV (Talos manual; MONO_CUT 0.8/0.09, so 0.1 A needs 4-8 pA, 0.2 A 15-32 pA),
  Leeds ~0.25 eV now linked, UltiMono 25 meV at 60 kV (no 300 kV spec found).
- version 12, 2026-10-21: simpler figures. Source figure: one bold line (the Mono) and one thin grey (no monochromator),
  labelled on the curves, dotted verdict lines, dots for the Mono's current at 0.1 / 0.2 Å, no legend or bars. fig_wall:
  the hatched "typical source size on real instruments" band removed (measured at 25-40 mrad with jitter; §06 has
  the context); its source axis says "at the sample".
- version 11, 2026-10-21: the user read the monochromator bars (Å of focus spread) as source sizes. The where-it-works
  figure now says so in its title ("a spread of FOCUS along the beam ... not the source size"), bars read "gives x Å of
  focus spread", rows read "source 0.1 Å at the sample (Mono: ≤ 4-9 pA at 80 mrad)", and its caption names the two Å.
- version 10, 2026-10-21: "focal spread" renamed "chromatic focal spread" everywhere (ΔE through Cc; ripple and
  lens-current noise named as other, unsimulated causes). New `fig_source_current.png` in §06 (`fig_source`): the gun's
  source image at the sample, d = (2/πα)√(I/B_r V), against probe current at 80 mrad; X-FEG B_r 0.7-1.5e8 (Thermo
  Fisher), Mono at 90 meV taken as ÷ 0.7/0.09 (an estimate); 0.1 Å needs ≤ 4-9 pA, 0.2 Å ≤ 17-37 pA; top axis the scan
  time for 1e7 e/Å². Where-it-works rows carry those currents. Probe modes defined exactly (true aberrations as the
  mixed state, no search, probe update from it 40 / 20, aperture constraint).
- version 9, 2026-10-21: `fig_where_it_works.png` at the top of §07 (`fig_where`, rows in `WHERE_ROWS`): every run with
  the source and vibration put back, at its focal spread, under the Spectra 300 monochromators only (Mono at 90 meV and
  at its 0.2 eV spec, UltiMono), on a linear 0-6 A axis. Hollow = probe fixed, filled = modes, grey = running. New runs
  appear on it once scored (add their keys to WHERE_ROWS).
- version 8, 2026-10-21: told as the story (the user): the ARM200F's realistic runs failed -> we scanned the coherence
  to find a microscope that could do it -> that turned up the Spectra 300 X-FEG Mono (90 meV in use, ~2.5 A focal spread;
  current and widely installed; S-CORR corrects A5) -> the physics put back at its setting. Sections reordered to match:
  §06 is now "Which microscope" (was §07, instruments), §07 "The physics put back" (was §06, A/B/C); every cross-
  reference moved with them. 0.1 / 0.2 A sources are not called realistic (real ones are quoted at 0.24-0.6 A).
- version 7, 2026-10-21: calmer formatting (the user: "a bit hard on my eyes"). One sans typeface for text; no all-caps
  monospace labels (eyebrow, nav, table heads, chips, box labels); call boxes as a left rule only; no figure shadows; the
  big headline panel and its seven tiles replaced by one "In short" sentence (`headline()` now returns only WALL_LINE);
  findings as a numbered list. Overrides sit in the template's own <style> (the ARM page's CSS is shared).
- version 6, 2026-10-21, rewritten for a reader who has not seen the runs (the user's PI): §01 opens with "What this
  page is" (question, test, score, aberrations, order), then the findings one limiting factor at a time, the main figure,
  then the run list in two parts and a terms table; every section opens with a question / answer box (`call qa`); table
  notes say what every column means (`COLS_NOTE`) and sit outside the scroll box (`p.tnote`, so wide tables never clip
  them); job ids off the page (they are in PHASE1.md); the Spectra Mono claim tightened: yes at the 90 meV quoted in
  use, not at its 0.2 eV specification. No plots removed.
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
