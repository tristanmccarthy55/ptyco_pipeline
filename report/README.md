# PX915 report

Self-contained manuscript pulling together both project phases: 3-D atom localisation with
calibrated UQ from multislice electron ptychography, the polarisation map it yields, and the
first-principles EELS study of the along-beam component it cannot measure.

## Build

```bash
latexmk -pdf report_v2.tex     # current draft: 9 pp main + 3 pp appendix + 2 pp refs
latexmk -pdf report.tex        # v1, kept for reference
```

## Files

| File | Role |
|---|---|
| `report_v2.tex` | **current draft** — condensed and restyled from v1 (23% shorter main text), technical detail moved to Appendices A-E |
| `report.tex` | v1, superseded; kept so the two can be diffed |
| `4th_Year_Physics_Project___Final_Report.pdf` | style reference (author's MPhys dissertation) |
| `refs.bib` | merged + deduplicated from `analysis/atomfind/paper/{atomfind,sim_recon}_refs.bib` and `eels/eels_refs.bib`, plus intro/context refs added for the report |
| `figs/` | figures (see provenance below) |

## Figure provenance

| Figure | Source |
|---|---|
| `sample.pdf`, `tornado.pdf` | `analysis/atomfind/paper/make_figs.py` |
| `fig_method.pdf`, `fig_uncertainty.pdf` | `analysis/atomfind/paper/make_paper_figs.py` |
| `fig_polarisation.pdf` | `analysis/atomfind/paper/make_pol_fig.py` (needs `polarisation.npz` from `atomfind/polarisation.py`) |
| `eels_dichroism.png`, `eels_decomp.png` | cropped from `~/Desktop/eels_figs/fig{2_result,3_scan}.png` (M4/M5 OptaDOS output; raw spectra live on Blythe) |

Unused-but-available in `analysis/atomfind/paper/figs/`: `fig_comparison` (Table II carries the
same numbers), `technique`, `ptycho_inverse`, `debye_waller` (methods schematics, cut for length).

## Conventions in v2

- **Appendix floats are numbered by appendix letter** (Table A1, C1, C2, D1) via
  `\renewcommand{\thetable}{\thesection\arabic{table}}` plus a counter reset in each appendix
  section, so they can never be confused with main-text Tables I and II. The one long-range
  reference (parameters table, cited in Methods) names its appendix explicitly.
- `\clearpage` before `\appendix` so the appendix starts on a fresh page.
- Float **declarations sit immediately before the paragraph that cites them**. LaTeX cannot place
  a float earlier than its declaration point, and `figure*` in two-column can only reach the top
  of a later page, so declaring after the citing text made every figure surface two pages late.
  Keep this ordering when editing.

## Numbers

All localisation figures come from the current run in `~/Desktop/atomfind_out/report.json`
(1834 atoms, precision 0.97), **not** the older values in `analysis/atomfind/RESULTS.md` §1,
which predate the guided-deduplication guard. Polarisation numbers are `RESULTS.md` §8.
Dose-series numbers are `RESULTS.md` §4 (pre-dedup), flagged as such in the Table III caption.

## Not yet in the report

- EELS M4 cross-checks: `tet_Px` rotational-invariance check, equatorial-O multiplicity
  weighting for the full O-K edge (both noted as further work in the text).
- The STEM-EELS forward simulator (`eels/simulate_stem_eels.py`) beyond the up/down-domain
  sign-blindness result.
