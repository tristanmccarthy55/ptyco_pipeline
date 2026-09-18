# analysis/paper — manuscript figures

Four figures, one claim each, all sharing one style. Everything is read from the campaign's own
outputs (`round_sweep.tsv`, `report.json`, `found_atoms.npy`, the `*_recons.h5`), so a figure
cannot drift from what was actually simulated and reconstructed — regenerate after any rerun and
the numbers follow.

This directory is the **manuscript** layer. The scripts one level up
(`make_ronchigram_fig.py`, `make_mep_volumes_fig.py`, `collate_atomfind_depth.py`) stay as
**diagnostics** and keep every check needed while working; the figures here import them for the
physics and the loaders, and then draw only what the claim needs.

## The figures

| | claim | inputs |
|---|---|---|
| `fig1_probe.py` | past ~70 mrad the Ronchigram stops being flat and the probe grows — conventional imaging is finished in exactly the regime the rest of the paper works in | `campaign/round_sweep.tsv` + `make_ronchigram_fig.py` |
| `fig2_depth.py` | atom recovery and depth accuracy both keep improving *through* that regime | `out/atomfind_a*/report.json` |
| `fig3_volumes.py` | at low α a column is one unbroken streak; by 90 mrad it separates into atoms at the right depths, and the blind finder puts them there | `recon_af_a*_lab_NL*/*_recons.h5`, `found_atoms.npy`, GT |
| `fig4_baselines.py` | on heavy columns any detector works; the matched single-atom kernel wins on the light atoms and on precision — and deconvolving first makes oxygen *worse* | `out/atomfind_a90/report.json` |

Fig. 1 is a 2×3 image block — Ronchigram over the wrapped aperture phase χ that produces it, at
three α — beside two line panels: the same phase as a ray displacement ∂W/∂θ (rays landing inside
±2 Å make a 4 Å probe), and flatness + d90 against α. The images show *what* and the profile shows
*why*; α is coded by ordered greys there, so the Okabe–Ito colours keep meaning species everywhere.

## Running

```bash
P=~/hyperspy-bundle/bin/python          # needs abtem for fig 1, h5py for fig 3
$P analysis/paper/fig1_probe.py
$P analysis/paper/fig2_depth.py  --glob '~/Desktop/thin_ab_af_final/out/atomfind_a*'
$P analysis/paper/fig3_volumes.py --recons ~/Desktop/thin_ab_af_final \
      --atomfind ~/Desktop/thin_ab_af_final/out --gt ~/Desktop/thin_ab_af/gtdata
$P analysis/paper/fig4_baselines.py --run ~/Desktop/thin_ab_af_final/out/atomfind_a90
```

Each writes `<name>.pdf` + `<name>.png` into `aberration_experiment/figs/<ISO week>/paper/`
(`--out DIR` to override). The PDF is the manuscript copy; the 600 dpi PNG is for slides.

## Style (`pubstyle.py`)

Width constants are journal columns — `COL1`/`COL15`/`COL2` = 85/127/180 mm — and every figure is
drawn at final size, so **do not rescale in the manuscript**: 8 pt text stays 8 pt. Colours are
Okabe–Ito (colour-blind safe) and each species keeps one colour across all four figures
(`SPECIES`: Pb blue, Ti orange, O green). PDF text is TrueType (`fonttype 42`), so it stays
editable and searchable.

Helpers: `panel(ax, "a")` for bold panel letters outside the frame, `imshow_clean()` for
tick-free image panels, `scalebar()`, `save()`.

## Conventions worth keeping if you add a fifth figure

- **One claim per figure**, stated at the top of the file as `ONE CLAIM:`. If a panel does not
  serve it, it belongs in a diagnostic script.
- **No hand-typed numbers.** Read them from the run outputs, even when it costs a loader.
- **Say when a number cannot be trusted.** Fig 2 draws open markers where atomfind's own
  species-confusion check exceeds 5 % — at α = 50 mrad the per-species split is label swapping,
  not physics, and the figure has to admit that.
- **Let the data pick the claim.** Fig 4 was originally drawn to show the matched kernel beating
  peak-picking outright; it does not (peak-picking gets Pb 100 % / Ti 97 % at 90 mrad). The figure
  was rewritten around what is actually true — oxygen and precision.
