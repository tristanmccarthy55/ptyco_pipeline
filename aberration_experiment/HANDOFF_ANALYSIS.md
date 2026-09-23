# Handoff — one job: a solid 90 mrad reconstruction with aberrations on the probe

Written 2026-09-23. **Read this whole file before running anything; it is short on purpose.**
Then `HANDOVER.md` for the experiment, `PSF_KERNELS.md` for the kernel rules.

---

## The job

Get a **good reconstruction at 90 mrad with six-fold astigmatism on the probe**, by reading the code and
fixing what is wrong with it. Nothing else. Do not extend the ladder, do not build tooling, do not write
pages. The sweep rows already exist (`nr90_C56_0p1w`, `nr90_C56_0p3w`, `nr90_C56_10um` in
`campaign/nonround_sweep.tsv`) and their simulations need regenerating after a storage clear-out.

## Why: the six-fold numbers are not believed, and should not be

**The physics.** The probe is exact, known, and held fixed. The solver is therefore fitting only the object.
A perfect probe is a perfect probe whatever its shape: the aberration is a known, deterministic part of the
forward model, and it should not limit the object the solver can reach. A 0.1-wave six-fold term costing
twenty points of lead recall is not a property of ptychography. It is a defect in this pipeline.

**What is on record, as a symptom rather than a result.** At 70 mrad, known probe, noiseless, static atoms,
after the probe-orientation fix:

| six-fold | lead | titanium | oxygen | residual |
|---|---|---|---|---|
| none | 98 % | 82 % | 75 % | 22.6 |
| 0.10 waves | 90 % | 71 % | 65 % | 37.8 |
| 0.20 waves | 76 % | 46 % | 15 % | 54.7 |
| 0.45 waves | 56 % | 36 % | 15 % | 74.0 |

`results/relaxation_ladder.csv` step 4. **Do not quote these as a tolerance.** They are what the pipeline
currently produces, and the pipeline is under suspicion.

**Two corroborating oddities.** The same idealised round-probe case gives residual 22.6 at 70 mrad but 5.8 at
90 mrad — a factor of four between two apertures that should both be near the model-mismatch floor. And
rotating a six-fold probe by 30° changed lead recall by up to 37 points while changing the residual by 4 %,
so the solve is violently sensitive to the probe in a way the data residual does not register.

## Leading suspect — the presolve solves a different problem

`ptycho/run_synthetic_recon_ML.m` runs two engines. The first has
`Np_presolve = 2*floor(Ndpx/4)`, half the detector width. `+engines/+GPU/+initialize/rescale_inputs.m`
then **crops the diffraction data** to that size and band-limits the probe to match.

A non-round aberration of order *n* scales as θ^(n+1). Halving the collection angle divides six-fold
astigmatism by **2⁶ = 64**. So the presolve solves an almost perfectly round problem, converges an object
under that assumption, and hands it to the full engine, which then has to reconcile it with a probe that has
a strong six-fold term. That is exactly the shape of the observed failure: worse with more aberration,
insensitive to the residual, and absent for round probes.

**Test it first, it is one job.** `PRESOLVE_NDP` (added 2026-09-23) forces the presolve's detector width;
set it to the full `Ndpx` and the downsampling disappears without touching the engine schedule. `Ndpx` is
356 at BIN 4 and 712 at BIN 2 — the recon log prints it as `sim_meta: Ndpx=...`.

```bash
cd /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline && git pull
# 0.45-wave leg, existing simulation, presolve at full resolution instead of half
PRESOLVE_NDP=356 TSV=campaign/nonround_sweep.tsv LABELS=nr0p45_C56_0p45w MODES=lab \
  RECON_ONLY=1 SAVE_EVERY=200 bash campaign/run_thin_atomfind.sh
```
The log must print `PRESOLVE_NDP: presolve detector width forced to 356`. Compare against the recorded 74.0,
and then — more importantly — compare the **object**, not the residual (see Traps).

## Other things worth checking in the code, in rough order

1. **Does 90 mrad show it at all?** BIN 2, 712 px, a 35 Å probe window instead of 17.5 Å. If six-fold costs
   nothing there, the problem is BIN 4 / 70 mrad conditioning, not the aberration. This is the job above.
2. **Probe window vs probe extent at BIN 4.** The reconstruction's real-space probe window is `box/BIN` =
   17.5 Å. A six-fold probe throws intensity into star arms; anything outside the window aliases back onto
   the object. Measured d99 is 10.5 Å at 0.45 waves so it nominally fits, but the *exit wave* after 27.5 Å of
   material is wider than the probe and nobody has checked that.
3. **Whether the full engine alone converges**, starting from a flat object rather than the presolve's.
4. **`grouping`** is `[64, 32]`; the LSQ update batches patterns, and a strongly structured probe may need
   smaller batches.

## What is solid — do not re-litigate

- **Focus.** An outer search over fixed trial probes recovers it and costs nothing (depth error 0.56 → 0.59 Å
  at 70 mrad). Fitting focus *inside* the solver fails and is closed.
- **Shot noise.** 10⁷ and 10⁶ e/Å² hold, 10⁵ works with the light atoms going, 10⁴ produces no measurable
  single-atom reference so the finder never runs.
- **The probe-orientation bug was real and the fix is right.** The engine transposes the diffraction data on
  load (`custom_data_flip = [0,0,1]`) and did not transpose the probe with it — invisible for a round probe,
  a 30° rotation for six-fold. `load_from_p.m` now flips the probe and its aperture mask with the data.
  Verified on objects: lead recall rose on all eight levels, mean +17.5 points.
- **The 70 Å sample fails**, at solver step 0.05 and 0.02. Phase saturates across every layer including the
  vacuum at 70 mrad; NaNs at 90 mrad. Separate problem, the user is picking it up later.
- **Aperture scaling** (`campaign/aberration_waves.py`) is analytic and depends on no reconstruction.

## How this project works

- **The user runs every cluster job.** No SSH. Give exact `sbatch` blocks; assume nothing until you see a log.
- **Everything lands under `$SHARE/phucrh`**, never the group root. Ship a matching `scp -O` line with every
  block, globbing the stable part of the tarball name, `tar` on its own line into a fresh directory.
- **The share is 3.9 TiB and shared with the department.** Submit with `CLEANDATA=1`. `SAVE_EVERY` now
  defaults to one save at the end. Check quota with `mmlsquota --block-size auto`, never `df`.
- Local Python is `~/hyperspy-bundle/bin/python`. Commit code, figures and results; leave the user's report
  and `analysis/atomfind/*.md` alone.
- The published record is one page, https://claude.ai/artifact/7ve93UM6yqCmiRcbfNuiJM, rebuilt from
  `aberration_experiment/page/logbook.html` by `page/build_page.py` and republished to that same URL.

## Traps

| symptom | cause | rule |
|---|---|---|
| **a probe change barely moves the residual, so it looks irrelevant** | a fixed-probe solve with a free object absorbs a wrong probe into the object: the fit hardly moves while the structure is wrong. Flipping the probe changed the residual 4 % and lead recall 37 points | **judge a probe change on the object, never on the residual.** This nearly caused a correct fix to be reverted |
| a leg passes triage and is still junk | it reported COMPLETED, wrote an h5 and logged no NaN, but the phase is saturated | `analysis/triage_recon.py` — a genuine leg wraps **zero** pixels; any wrapping at all is the flag |
| GPU jobs die anywhere with "Disk quota exceeded", naming a path that has no quota | `$HOME` is a 2 GB fileset and the CUDA JIT cache `~/.nv` fills it | both slurm scripts redirect every cache off `$HOME` and the recon aborts early if `$HOME` is unwritable |
| grepping a recon's slurm log reads the wrong run | recon dirs keep logs from every attempt; the driver moves only `analysis/` aside | check the job id, not the glob order |
| raising NL to chase a residual makes it worse | NL is Nyquist for the depth resolution; more layers with `REGLAYER=0` destabilises the solve (NL 28 gave 739 against 74 at NL 14) | leave NL alone; the `NL` override is for diagnosis |
| a comparison figure's panels differ for no reason | one RNG shared across panels | `make_ronchigram_fig.ronchigram` seeds per call |
| six jobs pending on `DependencyNeverSatisfied` | a simulation whose output exists exits 1 in zero seconds | `RECON_ONLY=1` when only the solver changed |
| `scp` finds nothing although jobs finished | the tarball name carries the submission timestamp | glob the stable part, never the date |

**Two rules that are physics, not plumbing:** `REGLAYER=0` on every leg (it low-passes the depth axis, which
is the measurement), and one `BETA_LSQ` for every leg in a comparison.

## State

`origin/main` plus whatever is unpushed locally — check `git log origin/main..HEAD`. Blythe needs `git pull`.
Raw simulation data was deleted in the storage clear-out, so any leg not currently on disk needs
re-simulating (about two hours of GPU at BIN 4 for a nine-label block).
