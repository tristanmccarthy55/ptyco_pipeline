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

> **Results before 2026-09-23 are wrong; a known non-round probe costs nothing.** The engine did not transpose
> the probe with the data (fixed 2026-09-22, `load_from_p.m`), and the fix reached only the presolve because
> `save_to_p.m` handed the full engine a flipped probe it flipped back (fixed 2026-09-23, `40d28a3`). On the
> fixed engine six-fold at 0.45 waves reconstructs like round at 70 and 90 mrad (ladder step 4).

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
(`$SHARE/$USER/<campaign>_results_<timestamp>.tgz`, i.e. your own subdirectory, not the group root). The pack job is `afterany`, so a partial sweep
still comes down.

Useful env overrides (defaults tuned for blind-fit convergence): `NITER=200 PSTART=40 BETA=0.05
THIN=3 STEP=0.5 SLICE=0.9`. `NL` per-row (Nyquist) unless you set `NL=` to force one value. `PMODES=1`
by default (coherent data ⇒ single aberrated probe; controls stay 1 mode) — `PMODES=2` gives only
the blind-fit leg extra incoherent freedom.

**Blind-fit machinery (fixed after the first campaign failed to recover any probe).** The fit leg
now (a) starts from an aberration-free **4 Å nominal** probe (`sim --probe-defocus df_perf`, BIN-matched
to the data at every α), (b) releases the probe **late** (`PSTART=40`) with a small step (`BETA=0.05`),
and (c) constrains the probe to the aperture in Fourier space (`PSFFT=1` → recon `PROBE_SUPPORT_FFT`),
without which the update absorbs high-frequency aliasing into a grid-artifact junk probe. The **perfect**
leg runs at its own `PERF_BIN=4` (a compact 4 Å probe is under-constrained in the aberrated legs' BIN=1
window). Set `PSFFT=0` to reproduce the original (failed) unconstrained fit.
Re-fit without re-simulating: `RECON_ONLY=1 CAMPAIGN=round bash campaign/run_campaign.sh`.
Re-simulate over an existing sweep: add `OVERWRITE=1`.

**Disk.** The raw 4D data is big and uncompressed (noiseless DPs are dense — gzip only ~1.2×):
~0.8 GB/point at BIN=4, ~3 GB at BIN=2, **~13 GB at BIN=1**; the full round sweep is ~78 GB. Add
`CLEANDATA=1` so the pack job deletes each point's `data_dp/position.hdf5` once results are tarred
(the recon `.h5` keeps object + probe + params). Old raw data is regenerable — reclaim it anytime with
`find <repo> \( -name data_dp.hdf5 -o -name data_position.hdf5 \) -delete`.

## C1 defocus search — `run_c1_search.sh` (aberration_experiment step 1)
C3 and C5 fixed at the sim's corrector values; **C1 found by an outer loop over fixed probes**. Each trial
reconstructs an existing `run_thin_atomfind.sh` sim (`sim_out_af_a<A>_<mode>`, no re-simulation) with a
probe that differs only in C1, and the objective is the solver's own final Fourier error.

- **Probe**: written *in the job* by `../sim/make_probe.py` (the recon slurm runs it when `PROBE_C1` is
  set) through the sim's own `build_initial_probe`, geometry from `sim_meta.mat`, C3/C5 from
  `aberrations.json`. At the sim's own C1 it must reproduce `probe_initial_true.mat` or the job stops
  before MATLAB (`SELF-CHECK PASS` in the log).
- **Objective**: `ptycho/run_synthetic_recon_ML.m` writes `<run>_error_trace.csv` (iteration, Fourier
  error of the full-resolution engine) and `<run>_layer_stats.csv` (per-layer phase mean/std over the
  illuminated field) beside the h5. Before 2026-09-17 neither the h5 nor the log carried the error.
- **Preflight**: every (α, mode) is checked before anything is submitted; a sim whose
  `aberrations.json` disagrees with `round_sweep.tsv` (an older a70 sim had C3 −4 µm, C1 +2 Å) aborts.
- **Dirs**: `recon_c1_a<A>_<mode>_df<C1>[_rN]_n<NITER>_NL<NL>`; a C1 listed twice runs twice (`_r2`,
  the noise floor). Engine constants pinned: `REGLAYER=0 BETA_LSQ=0.05 PROBE_MODES=1`.

```bash
ALPHAS=70 DC1="0" bash campaign/run_c1_search.sh                         # smoke test (true C1)
ALPHAS=70 DC1="-60 -50 -40 -30 -20 -10 -10 -8 -6 -4 -2 0 0 0 2 4 6 8 10 10 20 30 40 50 60" \
    NITER=50 bash campaign/run_c1_search.sh                                # objective grid (25 trials)
ALPHAS=70 C1="<fit>" MODES="lab Pb Ti" NITER=200 bash campaign/run_c1_search.sh   # final + matched kernels
DRYRUN=1 ... bash campaign/run_c1_search.sh                               # print the sbatch lines only
```
**Probe update (stage 2.5):** `PSTART=40` releases the probe in the presolve engine, `PSTART2=20` also in the
full engine, with the aperture constraint on; the trial probe is then only the starting point. Campaign
defaults to `c1fit`, dirs gain `_ps40[x20]`, and `analysis/probe_refit_check.py` summarises the refined probes.
**Defocus-only update (stage 2.5b):** add `PDFO=1` — the engine option `probe_defocus_only` moves the probe's
defocus alone (C3/C5 pinned, exact propagation per step); campaign `c1dfo`, dirs `_ps40x40dfo`, the engine's
own cumulative shift in the error-trace header.

Env: `ALPHAS`, `DC1` (offsets from the TSV C1) or `C1` (absolute), `MODES`, `NITER` (50), `PSTART`, `PSTART2`, `PSFFT`, `PACK_H5`
(1; 0 packs sidecars, probes and logs only), `DRYRUN`, `SAVE_EVERY` (= NITER), `BETA_LSQ`, `SIM_ROOT`.
Tarball: `$SHARE/$USER/c1_results_a<alphas>_n<NITER>_<ts>.tgz`, packing only that submission's recon dirs
(listed in `logs/c1_pack_*.dirs`), so several submissions can run side by side. Analyse locally:
```bash
~/hyperspy-bundle/bin/python analysis/c1_objective.py --root ~/Desktop/<fresh dir> --blind-start -40 -20 20 40
```

## Relaxation legs through the known-probe driver
`run_thin_atomfind.sh` also runs the later relaxation steps on the same pipeline (usage lines in its header):
`DOSES` (shot noise, Poisson copies), `PHONONS=16 PER_SPECIES=1` (frozen phonons, dirs `_ph16`),
`THIN=18 CELL_Z=3.889 GROUPING=16 RTIME=20:00:00` (the 70 Å slab, dirs `_thin18`), and
`TSV=campaign/nonround_sweep.tsv LABELS="nr1_C56_0p6w ..."` (non-round rows selected by label, dirs named by
the label). Each submission packs only its own recon dirs into a tarball named after them.

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
