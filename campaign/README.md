# Aberration-retrieval campaigns (thin slab)

Two parameter-sweep campaigns that stress-test blind probe retrieval in electron ptychography,
built on the fast, stable thin-slab setup de-risked in `../run_thin_aberration.sh` (3-cell
PTO/STO slab, BIN=4, few layers, ~2 min/leg, no NaN). Each sweep fans out the 3-leg comparison
per point and packs every result into **one tarball for a single scp**.

## The three legs (per sweep point)
| leg | probe | what it tells you |
|-----|-------|-------------------|
| `perfect`   | aberration-free, defocused to 4 Å (`df_perf`) | resolution ceiling at this α |
| `ab_known`  | the TRUE aberrated probe, fixed               | is the data usable if the probe is known? |
| `ab_fitprobe` | starts from the perfect probe, **updates** | the blind-retrieval test |

`ab_fitprobe ≈ ab_known ≈ perfect` ⇒ ptychography retrieved the aberration.

## Campaign A — round α-sweep (`round_sweep.tsv`)
A microscope with a FIXED residual **C5 = 1 mm** (uncorrectable 5th-order spherical). Opening
the aperture α grows its effect as α⁶. Two knobs hold the probe at a constant **4 Å** so α is the
only variable: **C3 (Cs) coarse @1 µm + C1 (defocus) fine**. Built by the planner:

```
~/hyperspy-bundle/bin/python campaign/plan_probe.py      # writes round_sweep.tsv (run locally, commit it)
```

The 2-knob balance holds 4 Å through **70 mrad**; past that C5 overruns it (d90 6.6/11/24/47 Å at
90/100/110/120) and BIN steps 4→2→1 — this is the "sweep until it breaks" axis. α120 (d99 > 70 Å
box) is past the wall; keep it only to show the failure.

**Depth sampling.** Each point reconstructs at **NL layers = Nyquist of its depth resolution**
δz = λ/α² (slice = δz/2): NL = 1/3/6/10/12/14/17 for α = 30/50/70/90/100/110/120. So low α can't
see the ~2 Å interatomic depth spacing and high α can — that difference *is* the result. For it to
be real the **sim slab is fixed fine at 0.9 Å** (< the 1.94 Å plane spacing; a 2 Å slab merges the
planes, erasing the depth structure), so only α + NL decide what's recovered. NL is per-row in the
`.tsv`; `SLICE` (sim slab) is the global fine default.

## Campaign B — non-round sweep (`nonround_sweep.tsv`)
Fixed α = 70 mrad, round part balanced as in A, then **non-round terms added on top, worst first**:
`C56` (6-fold astigmatism — the hexapole-corrector residual, uncancellable by round C1/C3),
escalated 0.6 → 1.2 → 2.5 waves, then a combined C56+C34 (4-fold) level. All stay BIN=4.

## Run
```bash
cd <repo>
CAMPAIGN=round    bash campaign/run_campaign.sh      # 7 α × 3 legs + pack job
CAMPAIGN=nonround bash campaign/run_campaign.sh      # 5 levels × 3 legs + pack job
```
Each prints a pack-job id and the exact `scp -O ...` line for the tarball
(`$SHARE/<campaign>_results_<timestamp>.tgz`). The pack job is `afterany`, so a partial sweep
still comes down.

Useful env overrides (defaults = fit-probe HARD): `NITER=200 PSTART=8 BETA=0.1 THIN=3 STEP=0.5
SLICE=0.9`. `NL` per-row (Nyquist) unless you set `NL=` to force one value on all points. `PMODES=1`
by default (coherent data ⇒ recover a single aberrated probe; controls stay 1 mode regardless) —
set `PMODES=2` to give only the blind-fit leg extra incoherent freedom.
Re-fit without re-simulating: `RECON_ONLY=1 CAMPAIGN=round bash campaign/run_campaign.sh`.
Re-simulate over an existing sweep: add `OVERWRITE=1`.

## `.tsv` schema
TAB-separated; `#`/header/blank lines skipped; the driver reads the **first 9** columns and
ignores the rest (planner appends d50/d90/d99/note as diagnostics):
```
label  alpha  c5  c3  c1  df_perf  bin  nl  aber_json
```
`c3/c1/c5` round knobs [Å]; `df_perf` aberration-free 4 Å reference defocus [Å]; `bin` from the
planner (probe size → real-space window); `nl` recon layers (Nyquist of λ/α² depth res); `aber_json`
`-` = round-only (use c3/c5), else a full abTEM Cnm/phi dict that **overrides** c3/c5 (non-round).

## Results / analysis
Each `*_recons.h5` holds the object volume (`reconstruction/object`, **NL×1×753×753** — NL varies
per α), the final probe (`reconstruction/probes`), the start probe and all params
(`reconstruction/p/...`). The true probe per point is in that point's `ab_known` h5. Plot with `../`
analysis (see `plot_thin_ab.py` pattern: dx = `p/dx_spec`·1e10 Å ≈ 0.049; **dz = 11.715/NL Å**, so
read NL from the object's first axis per file). PtychoShelves TIFFs are OFF (`save.store_images=0`);
render from the h5.

## Cleanup (revert the whole thing)
Delete `campaign/` and the `[campaign]` blocks in `../sim/simulate_4dstem.py` (the
`--aberrations-json` arg + its handler) and `../sim/run_sim.slurm` (the `AJ_ARG` line). The
`[thin-ab]` machinery it reuses is documented in `../run_thin_aberration.sh`.
