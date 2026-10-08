# The ARM200F results page — how to fill it in and republish

Built 2026-10-08 for the meeting of 2026-10-09. Published URL: see below.

PUBLISHED: https://claude.ai/artifact/MUkJ7yYrsek3nTLXcZgCrk (version 1, 2026-10-08 evening: everything to the search; refined, true-start and final runs pending.
Version 2, same evening: atoms scored as "right place" against a guessing level, the finder's own depth bars drawn, misplaced atoms highlighted -- see "The atom numbers" below.
Version 3: text trimmed ~40 % at the user's request; the 5 % kick section cut to one paragraph, one table, one figure. Keep new prose as short.
Version 4: appendix B "Reproduce" (part_reproduce): coherent vs full physics, every probe coefficient, the fixed settings and
each run family's launch line -- add a family to the manifest's "reproduce" list when a new kind of run lands.)

The page is generated, not hand-edited: figures, judged tables and tab viewers come from the run outputs; the prose
lives in the template. Runs that have not landed show as marked "pending" placeholders, so the page can be republished
at any point.

| file | what it is |
|---|---|
| `aberration_experiment/page/arm_results.json` | the manifest: every run (tag, leg, tab label, display name, pending note) and the Desktop folders searched |
| `analysis/make_arm_page.py` | builds every figure into `aberration_experiment/figs/arm_page/` and writes the page |
| `aberration_experiment/page/arm_2026-10-09.html` | the template: prose, figure slots (`{{FIG_*}}`, `{{FRAG_*}}`), and `FILL` blocks for pending results |
| `aberration_experiment/page/built/arm_page_built.html` | the page to publish (gitignored, rebuilt every time) |
| `aberration_experiment/results/arm_page_cache/` | region ground-truth caches, rebuilt on demand (gitignored) |

## When results land

1. **Pull** the tarballs into one of the folders in `search_dirs` (the pull blocks already name them: `~/Desktop/arm_true4`,
   `~/Desktop/arm_s2_1e7`, `~/Desktop/arm_block3_s2`, `~/Desktop/arm_final90`). Check the tar steps ran; extract any that
   did not (analysis tarballs into `<folder>/analysis/`, atomfind_results into `<folder>/results/`). A new folder name
   must be added to `search_dirs`: macOS blocks listing `~/Desktop` itself, so the builder only searches named folders.
2. **See what landed:** `~/hyperspy-bundle/bin/python analysis/make_arm_page.py --list` (IN / part / pending per run).
3. **Rebuild:** `~/hyperspy-bundle/bin/python analysis/make_arm_page.py`. It prints the page size; stay under 16 MB.
4. **Look before writing** (the user's rule): the phase images and `figs/arm_page/fig_refined.png`, `fig_final90.png`
   first, then the tables in the built page.
5. **Write the interpretation** into the template's `FILL` blocks (`grep -n FILL aberration_experiment/page/arm_2026-10-09.html`):
   `headline` (section 01 panel), `refined` (section 07), `final` (section 08), `discussion` (section 09). Replace the
   matching `call pending` box once its runs are in. Plain language, numbers from the tables, no new claims the figures
   do not show.
6. **Rebuild again and republish** with the Artifact tool: `action: publish`, `file_path` = the built page, `url` = the
   PUBLISHED url above (a conversation that did not publish it must `action: read` it first). No `icon` on a redeploy.
7. **Commit** the template, manifest and script changes, and the `figs/arm_page/fig_*.png` print copies when they change (never `built/`; the `tab_*` images are gitignored).

## Adding a run

Add an entry under `runs` in the manifest: `tag` (list; the first is the run with atoms, later ones are fallbacks that
supply only the volume), `leg` (the recon dir core name), `short` (tab label), `name`, and `pending` (job ids and the
folder it will be pulled to). Then use it in the relevant `part_*` function of `make_arm_page.py` (each section has one).

## Conventions kept

- One colour per species everywhere (Pb `#2a78d6`, Ti `#eb6834`, O `#1baf7a`); never for anything else.
- Real atoms are dots, found atoms small rings with the finder's own 95 % depth bar (`halfwidth95_z_A`), vacuum edges
  light-blue dashes; depth runs downward. Ring states (`atom_states()`): element colour = within 0.5 Å of a true atom of
  its element; faded = further, but its bar reaches one; white with a black edge = misplaced (its bar misses, or no atom
  of that element there). Depth sections take their colour scale from the slab's slices only (`xz_section()`).
- No number typed into a figure: everything comes from summary.csv, kick records, error traces, the sweep table, the
  plan JSONs or the measured tableau.
- Verdicts (`verdict()`), on Pb and Ti in the right place: good = Pb ≥ 80 % and Ti ≥ 45 %; okay = Pb ≥ 50 %; off
  otherwise. Oxygen is reported, not judged. Stated on the page under the setup table.
- Tab images and renders go in as JPEG at moderate width to keep the page small; figures stay PNG.

## The atom numbers (since version 2)

`atom_scores()` in make_arm_page.py, cached in `results/arm_page_cache/atom_scores.json` by found-atoms file and time:

- **right place**, per species over the bulk: share of real atoms with a found atom of the SAME element within 0.6 Å in
  plane and 0.5 Å in depth (the finder's greedy matcher, one element at a time).
- **guessing**: the same found atoms slid along their columns to 16 evenly spread depths over the 3.9 Å period, scored
  the same way and averaged. Pb by guessing is ~26 % of its in-plane hits.
- **wrong element**: the finder's own confusion figure (its 0.6 Å / 2 Å match, any element), as in its log warnings.
- **95 % depth bar**: median of the finder's calibrated per-atom half-widths.

Why: summary.csv's recall matches within 2 Å in depth and any element, so along a Pb column almost any depth passes;
the full-physics runs looked "okay" on it with no depth information (aberration_experiment/results/2026-W41/atomfind_honesty/).
Runs that shipped no volume (block 2's kick legs) are scored in the frame recovered from their own found_atoms.csv,
accepted only when it reproduces their summary.csv recall within 2 points.

## Known gaps

- Block 2's kick legs shipped no volume (the kick pass saves none): they appear in the kick table only (scored from their
  own export), with full renders in the appendix.
- The 1e7 no-search baseline has atoms only once analysis 1310612 is pulled (tag `..._rk_kern-..._dose1e7_fixed`); until
  then it shows its volume only.
- The kick-record overlap compares with the central coherent probe; under partial coherence even a perfect recovery
  scores below 1 (said on the page).
