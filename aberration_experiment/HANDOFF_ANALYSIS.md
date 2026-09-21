# Handoff — analyse the 2026-09-22 runs, then finish the ladder

Written 2026-09-22 for whoever picks this up next. **Read `HANDOVER.md` first** for the experiment
itself, then `NEXT_PHASE.md` for the plan, then `PSF_KERNELS.md` for the kernel rules. This file is
only: what is waiting for you, how to analyse it, and what is left before the result is publishable.

Your job, in order:

1. **Analyse the six tarballs** from the 2026-09-22 submissions (below). They should be pulled to
   `~/Desktop/relax_0922/` by the time you start.
2. **Put the results in front of the group** — update the two pages listed under *Meeting pages*.
3. **Finish the ladder.** Two relaxations still need simulation code that does not exist, and the
   combined run has never been attempted. That is the road to publishable.

---

## How this project works

- **The user runs every cluster job.** You have no SSH. Give exact `sbatch` command blocks and the
  user pastes back the output. Assume nothing about a job until you see its log or its tarball.
- **Every Blythe path lands under `$SHARE/phucrh`** (`/springbrook/share/physics/phucrh/`), never
  the group root `$SHARE`, which other people use. After each submission ask the user to check
  `ls /springbrook/share/physics/` shows nothing new of ours.
- **Give a matching `scp -O` line with every command block you hand over** — one per tarball, into a
  fresh directory, with `tar xzf` on its own line. Tarball names carry a timestamp you cannot
  predict, so glob the stable part. The user asked for this explicitly.
- **Python is `~/hyperspy-bundle/bin/python`** for everything local (abTEM, h5py, scipy, matplotlib).
- **Commit code, figures and results tables to `origin/main`.** Leave the user's report and the
  `analysis/atomfind/*.md` files alone; the `.py` files there are fair game.
- **Verify against the code or the data before asserting anything.** Every wrong number in this
  campaign came from an assumption that looked right. Two examples are in *Traps* below.

---

## What is waiting for you

Six tarballs, from three experiments submitted 2026-09-22. Expected names (timestamps will differ):

| experiment | tarball glob | what it answers |
|---|---|---|
| **A** focus fitted with the object frozen | `c1dfo_results_a70_n200_ps40x1dfooinf_*` and `..._a90_...` | is the focus-fit runaway the focus/depth degeneracy? |
| **A2** focus fitted with a capped step | `c1dfo_results_a70_n200_ps40x40dfoc0.05_*` | control: confirms a step cap is *not* the fix |
| **B** the 70 Å sample, gentler solver step | `atomfind_results_a70_thin18_*`, `atomfind_results_a90_thin18_*` | does the thick sample reconstruct at `BETA_LSQ=0.02`? |
| **C** non-round threshold ladder | `atomfind_results_nr0p1_C56_0p1w-nr0p2_...-nr0p45_C56_0p45w_*` | **the key number**: how much six-fold astigmatism is tolerable? |

### A and A2 — the focus fitted inside the solver

Background: the outer search over trial probes recovers the focus to 2 Å and is the working route.
The tidy alternative, fitting focus as the solver's one free probe parameter, was built into the
engine (`p.probe_defocus_only`) and **failed**: at 90 mrad it moved the focus under 0.6 Å from starts
32 Å out; at 70 mrad, from starts at or above the truth, it walked the probe 20–33 Å the wrong way
into a reconstruction with mismatch ~260 against ~22 for the best fixed probe.

The first diagnosis (overshoot from oversized steps) was **wrong** and is disproved in the data:
median step 0.023 Å, only 3% above 0.5 Å, 54% of steps one-signed over 322 iterations. That is a slow
steady drift. The live hypothesis is the degeneracy the fixed-probe scan already measured: probe focus
and object depth trade one-for-one, so with both free the pair slides along a direction the data
cannot see, until the object runs out of box.

**A tests it**: `OSTART2=inf` freezes the object in the full-resolution engine, so the focus is fitted
against the object the presolve already converged.

- If the focus now **moves toward the truth** → the drift was the shared degeneracy. The fix is to
  alternate (converge object, then fit focus, then repeat) rather than co-optimise, and the in-solver
  route becomes viable. Say so plainly; it changes the method.
- If the focus **still does not move, or still drifts** → the one-parameter gradient is biased on this
  data and the outer search is the permanent answer. Also a clean result. Record it and move on;
  do not spend days on it.

The solver now logs, every tenth iteration, the cumulative focus shift **and** the sample's depth
centroid:

```
[GPU-1_MLs] : defocus-only probe update: cumulative C1 shift -6.56 A | sample depth centroid 17.41 A
```

If during a drift the centroid tracks the shift one-for-one, that is the degeneracy caught in the
act, and it is worth a figure on its own.

```bash
# read the trajectory out of any focus-fit run
grep -h "defocus-only probe update" <recon dir>/slurm_*.out | tail -20
```

### B — the 70 Å sample

The previous attempt diverged in five of six legs (at 90 mrad two thirds through the first pass, at
70 mrad on entry to the second) with the step size pinned at `BETA_LSQ=0.05`. This retry halves it.
If it reconstructs, **this is the biggest remaining step**, because thickness is the regime the method
exists for. Note it needs analysis work that does not exist yet — see *Still to build*.

### C — the non-round threshold

This is the most publication-relevant run. Six-fold astigmatism is what a hexapole corrector leaves
behind and round knobs cannot touch. Measured so far at 70 mrad: 0.6 waves takes lead recall from 98%
to 39%, 1.2 waves to 7%. So the tolerance is somewhere **below** 0.6, and C samples 0.1, 0.2, 0.3 and
0.45 waves. The output is a specification: *the residual must be under X waves for this to work.*
Quote it with the recall and depth error at each rung, and say plainly if even 0.1 waves hurts.

---

## How to analyse, step by step

### 1. Triage before anything else

```bash
cd ~/Desktop/relax_0922
for d in */; do echo "=== ${d%/}"; for r in $(find "$d" -maxdepth 2 -type d -name "recon_*" | sort); do
  h=$(find "$r" -name "*_recons.h5" | wc -l | tr -d ' '); L=$(ls "$r"/slurm_*.out 2>/dev/null | tail -1)
  printf "   %-54s h5=%s nan=%s\n" "$(basename $r)" "$h" "$(grep -c 'contains NaNs' $L)"; done; done
```

A job that shows COMPLETED in slurm but has **no `*_recons.h5` means MATLAB crashed**. Check the log
for `contains NaNs` and for `Number of probe positions: 1600`. Never analyse a tarball you have not
triaged; an older run's h5 at the same path will silently stand in for a leg that failed.

### 2. The focus-fit runs (A, A2)

```bash
cd <repo>
~/hyperspy-bundle/bin/python analysis/probe_refit_check.py --root ~/Desktop/relax_0922 --camp c1dfo
```

Prints, per run: whether the probe survived, its overlap with the true probe before and after, the C1
measured from the refined probe, the engine's own reported shift, and where the sample ended up. The
measured C1 and the engine's shift should agree; they did last time, which is how we know the
machinery is sound. Slow (a few minutes per a90 run) because it rebuilds probes with abTEM.

For the figure, `analysis/make_relaxation_figs.py` already has the layout (`figA_dfo.png`): focus out
against focus in, plus the same moves drawn as arrows on the fixed-probe landscape. Point its
`dfo_runs()` at the new root and add the frozen-object runs as a second series.

### 3. The reconstruction campaigns (B, C)

Each leg needs its own matched kernels, then the blind finder. The pattern is in
`/private/tmp/.../scratchpad/run_0921_atomfind.sh` from the previous session; rewrite it for the new
directory names. Per leg:

```bash
# kernels from that leg's own Pb and Ti grid recons (--zdrop = round(4 Å / dz): 2 at NL14, 3 at NL23)
~/hyperspy-bundle/bin/python analysis/atomfind/extract_psf.py <recon_..._Pb_...> Pb_<tag> --zdrop 2 --out <out>/psf
# then the finder on that leg's lab recon, with those kernels
~/hyperspy-bundle/bin/python analysis/atomfind/run_atomfind.py --preset thin --recon <lab h5> \
    --dz $(python -c "print(27.525/14)") --data-dir ~/Desktop/thin_ab_af/gtdata \
    --single-atom-vol <out>/psf/psf_Pb_<tag>_vol.npy --ti-kernel-vol <out>/psf/psf_Ti_<tag>_vol.npy \
    --out <out>/atomfind_<tag>
```

**Kernel health matters as much as the finder's output.** `extract_psf` prints peak, background and
peak/background; a degraded kernel means the *reconstruction* is failing, not just the finder. That
distinction is what made the non-round result interpretable (peak fell fourfold, background rose
fivefold), so report it.

Then one ladder row per leg:

```bash
~/hyperspy-bundle/bin/python analysis/relaxation_ladder.py --step 4 \
    --label "six-fold astigmatism, 0.3 waves" --alpha 70 \
    --atomfind <out>/atomfind_nr0p3_C56_0p3w --note "alpha 70 mrad, known probe, noiseless"
```

`--alpha` is required for legs whose directory name carries no aperture (the non-round labels). The
ladder is **one cumulative table**, `aberration_experiment/results/relaxation_ladder.csv`, 16 rows as
of this handoff. Rows are replaced, never duplicated, on `(step, label, alpha)`.

### 4. Figures

Two scripts, same house style (`make_relaxation_figs.py` imports its palette and helpers from
`make_meeting_figs.py`, so extend rather than restyle):

- `analysis/make_meeting_figs.py` — the six baseline figures: probe vs aperture, the depth section
  with found atoms, accuracy vs aperture, the focus search, the probe-update failure, the dose ladder.
- `analysis/make_relaxation_figs.py` — figures A–D: focus fitted in the solver, focus fitted vs known
  probe, phonons, non-round.

Both write to `aberration_experiment/figs/<ISO-week>/meeting/`. **Every value is read from the run
output; nothing is typed into the scripts.** Keep it that way — it is why the pages and the data
cannot drift apart.

The user's standing instruction on figures: *no useless confusing detail.* One quantity per axis,
pick one probe-size measure and call it "probe diameter", and put a small "what is in this figure"
explanation next to each one rather than cramming it into the axes.

---

## Meeting pages

Two published artifacts. **Update them, do not make new ones**: pass the URL as `url`, read it first,
and build your update on what comes back.

| page | URL | covers |
|---|---|---|
| Depth Past the Corrector | `https://claude.ai/artifact/NiYCfNo7uFAWFkyUZ5L3SF` | the baseline story: probe, depth vs aperture, focus search, dose |
| Removing the Six Assumptions | `https://claude.ai/artifact/3TfTwcipbbTxG9PGbeEtqi` | the 2026-09-21 results: focus end to end, in-solver failure, phonons, non-round, thick slab |

**Embed figures as base64 data URIs inside the HTML.** Publishing them as separate artifact files and
referencing them by relative path does not render, even though the files list correctly — this cost a
round trip on 2026-09-21. Budget: page under 16 MB, base64 inflates by a third.

House style of both pages, which the user approved: dense prose in a reading column, numbered
sections, wide figure blocks with a small "how to read this" note, a status ledger with state chips,
and tables with the baseline row tinted. Plain language, professor-level audience, no buzzwords, and
only claims the user could defend out loud.

---

## Where the ladder stands

From `aberration_experiment/results/relaxation_ladder.csv` (α 70 / α 90 where both exist):

| step | relaxation | outcome |
|---|---|---|
| 0 | known probe, no noise | depth error 0.56 / 0.37 Å; recall 98/82/75 and 95/94/95 % |
| 1 | focus unknown, corrector known | **no real cost**: 0.59 / 0.44 Å, recall 97/94/79 and 90/83/86 % |
| 1b | focus fitted inside the solver | **fails**; experiment A is the diagnosis |
| 2 | counting noise | works to 10⁵ e/Å²; at 10⁴ the matched kernel cannot be measured at all |
| 3 | thermal vibration | depth error roughly doubles; oxygen recall 10–17 %; species labels untrustworthy |
| 4 | non-round aberration | 0.6 waves → lead 39 %; 1.2 waves → 7 %. Experiment C finds the threshold |
| 5 | the 70 Å sample | diverged at the pinned step; experiment B retries gentler |
| 6 | partial coherence | **not built** |
| 7 | specimen tilt | **not built** |

---

## Still to build, in the order that matters

1. **Analysis for the thick sample.** If B reconstructs, atomfind cannot be pointed at it as-is:
   - a ground-truth cache for the 18-cell slab: `make_gt_cache --thin-cells 18 --z-vacuum 4`;
   - a `thick` preset in `analysis/atomfind/config.py` with `trim_z_A`, `zmax_show_A` and
     `clean_max_atoms` **derived from the box**, not typed (three bugs in this campaign came from
     constants typed for one geometry and silently wrong for the next);
   - `--dz` from `sim_meta.beam_thickness_A / NL`, not from the driver's box formula. The measured
     18-cell box is **77.94 Å** (atom span 69.94 Å), against 78.00 from the driver's formula.
2. **Partial coherence (step 6).** Needs new simulation code: an incoherent sum of diffraction
   patterns over a Gaussian spread of focus (cold FEG at 300 kV gives roughly 12–23 Å at 1σ, 5–7
   quadrature points), then over source position (σ ≈ 0.2–0.4 Å). The reconstruction side already has
   `PROBE_MODES` and `VARIABLE_PROBE` to absorb it.
3. **Specimen tilt (step 7).** A few degrees, so 1/2/3/5°. Needs the crystal rotated about an in-plane
   axis in the simulation (not the probe tilted), a ground truth built from the *same* rotated atoms,
   and a finder that handles columns which cross one another in projection — on the 70 Å slab the
   shear passes the 1.95 Å column spacing at about 2°. Validate on the thin slab first, where the
   shear stays below that at every tilt on the ladder.
4. **The combined run.** Every surviving relaxation switched on together, at 50/70/90/100 mrad, with
   matched kernels and the same figure scripts. A result that survives the combined configuration is
   publishable; one that survives only the idealised case is not.

Open questions worth answering along the way, cheaply:

- **Are the phonon and non-round recall collapses partly the finder's thresholds?** They were tuned on
  still atoms and clean wavefronts. Re-running the *existing* reconstructions with a relaxed quality
  cut costs nothing on the cluster and would separate physics from tuning. Do this before concluding
  the method fails at room temperature.
- **False positives rise with dose even at 10⁷ e/Å²**, almost all oxygen. Same suspicion, same test.
- **a90's focus interval is narrower than the search grid** (fit −158 ± 1 Å against a true −160). A
  0.5 Å scan over ±5 Å would pin it, about 6 GPU-hours, only worth it if a referee would ask.

---

## Traps

Each of these cost real time. They are fixed in the code; do not reintroduce them.

| symptom | cause | rule |
|---|---|---|
| figures missing from a published page | images published as separate files, referenced relatively | embed as base64 data URIs |
| tarball downloads but the directory is empty | `tar` chained after `scp` on one line gets skipped | `tar` on its own line, into a fresh directory |
| every a90 result labelled a70 | aperture parsed from a tarball name containing `a70-90` | `relaxation_ladder.py --alpha` |
| finder crashes with "array must not contain infs or NaNs" | a column detected at the crop edge averaged an empty window | fixed in `find.py`; verified byte-identical on the published noiseless case |
| aperture-constrained probe fit NaN'd at iteration 1 of engine 2 | the aperture mask was resized by keeping the array centre, which for an unshifted FFT is pure high frequency, so the presolve mask was all zeros and erased the probe | fixed in `load_from_p.m` / `init_solver.m` |
| driver exits silently with no message | `VAR=$([ test ] && echo x)` fails the substitution when the test is false, and `set -e` kills the script | use an `if`, never that idiom |
| kernel quality numbers disagree with the pipeline's | `extract_psf` saves `exp(i·phase)`; take the phase first | see `kernel_quality()` in `make_meeting_figs.py` |

Two rules that are physics, not plumbing, and must hold on every leg: **`REGLAYER=0`** (it low-passes
the depth axis, which is the quantity being measured) and **one `BETA_LSQ` for every leg in a
comparison** — which is exactly why experiment B changes it deliberately and only for the thick slab,
and why its result must be compared against a thin-slab leg at the same step size before anything is
concluded about thickness.
