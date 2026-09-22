# Handoff — finish the ladder, and turn two pages into one

Written 2026-09-22, evening, for whoever picks this up next. **Read `HANDOVER.md` first** for the
experiment itself, then `NEXT_PHASE.md` for the plan, then `PSF_KERNELS.md` for the kernel rules.
This file is only: what exists and what does not, how to analyse it, and what is left before the
result is publishable.

Your job, in order:

1. **Analyse what is in hand** — three experiments finished and are unanalysed (*What you have*).
2. **Build one page, not two.** The user asked for a single, better artifact in place of the two
   published now. See *The page*.
3. **Finish the ladder.** Two relaxations still need simulation code that does not exist, the thick
   sample has no analysis path, and the combined run has never been attempted. That is the road to
   publishable.

Nothing here is blocked on the user except cluster submissions and the git push noted at the end.

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

## What you have, and what you do not

### In hand, unanalysed

Four tarballs from the batch submitted 2026-09-21 at 22:56. **Note the stamp**: the jobs went in late
on the 21st, so the names read `20260921_2256*` and not the 22nd. A glob on the wrong date made them
look missing and cost a round trip. They are on Blythe under `/springbrook/share/physics/phucrh/`:

| experiment | tarball | size | what it answers |
|---|---|---|---|
| **A** focus fitted, object frozen, 70 mrad | `c1dfo_results_a70_n200_ps40x1dfooinf_20260921_225613.tgz` | 314 MB | is the focus-fit runaway the focus/depth degeneracy? |
| **A** the same at 90 mrad | `c1dfo_results_a90_n200_ps40x1dfooinf_20260921_225613.tgz` | 870 MB | as above, at the aperture that matters most |
| **A2** capped-step control | `c1dfo_results_a70_n200_ps40x40dfoc0.05_20260921_225613.tgz` | 128 MB | confirms a step cap is *not* the fix |
| **C** non-round threshold ladder | `atomfind_results_nr0p1_C56_0p1w-…-nr0p45_C56_0p45w_20260921_225614.tgz` | 855 MB | **the key number**: how much six-fold astigmatism is tolerable |

If `~/Desktop/relax_0922/` is empty, hand the user this:

```bash
mkdir -p ~/Desktop/relax_0922 && cd ~/Desktop/relax_0922
scp -O phucrh@blythe.scrtp.warwick.ac.uk:'/springbrook/share/physics/phucrh/*_20260921_2256*.tgz' .
for t in *.tgz; do d="${t%.tgz}"; mkdir -p "$d" && tar xzf "$t" -C "$d"; done
```

Already on the Mac, and needed for every comparison: `~/Desktop/relax_0921/` (the 2026-09-21 12:32
batch — phonons, the first non-round rungs at 0.6/1.2/2.5 waves, the first 70 Å attempt, and the
focus-fit runs with a free object) and `~/Desktop/thin_ab_af_final/` (step 0, the known-probe
baseline).

### Not in hand

- **Experiment B, the 70 Å sample at `BETA_LSQ=0.02`.** It did not run on 2026-09-22 and was
  resubmitted the same evening. Read *Experiment B* below before assuming anything about it, and
  check for `atomfind_results_a70-90_thin18_20260922_*.tgz` rather than trusting this file.
- **Steps 6 and 7**, partial coherence and specimen tilt. The simulation code does not exist.
- **An analysis path for the thick sample.** atomfind cannot be pointed at an 18-cell slab as it
  stands, so even a successful B is not immediately readable. See *Still to build*, item 1.
- **The combined run.** Never attempted, and it is what makes the result publishable.

What each of those experiments is for, and what a result would mean, follows.

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

### B — the 70 Å sample (did not run; resubmitted)

The attempt of 2026-09-22 failed before it started, and the failure is worth understanding because it
will recur otherwise. `sim/run_sim.slurm` refuses to overwrite a finished simulation and exits 1 in
zero seconds. Those 70 Å simulations were already on disk from the earlier batch, and the only thing
this run changed was the solver step size, so all six simulation jobs died instantly and the six
reconstructions queued behind them sat on `DependencyNeverSatisfied` for 90 minutes until they were
cancelled. The run never needed to simulate anything.

It was resubmitted as reconstruction-only, reusing those simulations:

```bash
cd /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline
RECON_ONLY=1 ALPHAS="70 90" THIN=18 CELL_Z=3.889 GROUPING=16 RTIME=20:00:00 \
  BETA_LSQ=0.02 bash campaign/run_thin_atomfind.sh
```

Six reconstructions, 39 layers at 70 mrad and 64 at 90, plus one pack job. Pull with:

```bash
mkdir -p ~/Desktop/thin18_0922 && cd ~/Desktop/thin18_0922
scp -O phucrh@blythe.scrtp.warwick.ac.uk:'/springbrook/share/physics/phucrh/atomfind_results_a70-90_thin18_20260922_*.tgz' .
```

The driver now refuses to submit a simulation whose output already exists, so it fails on the login
node with both escapes named instead of on the GPU queue. That guard is in an unpushed commit; see
*State of the repository*.

The previous attempt at `BETA_LSQ=0.05` diverged in five of six legs and is already on the Mac in
`~/Desktop/relax_0921/`. If this retry reconstructs, **it is the biggest remaining step**, because
thickness is the regime the method exists for. Remember the rule at the foot of this file: a thick
result at 0.02 must be compared against a thin leg at 0.02 before anything is concluded about
thickness.

### C — the non-round threshold

This is the most publication-relevant run. Six-fold astigmatism is what a hexapole corrector leaves
behind and round knobs cannot touch. Measured so far at 70 mrad: 0.6 waves takes lead recall from 98%
to 39%, 1.2 waves to 7%. So the tolerance is somewhere **below** 0.6, and C samples 0.1, 0.2, 0.3 and
0.45 waves. **The data is in hand and unanalysed** — this is the first thing to do. The output is a specification: *the residual must be under X waves for this to work.*
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

- `analysis/make_meeting_figs.py` — figures 1–7: probe against aperture, the depth section with
  found atoms, accuracy against aperture, the focus search, the probe-update failure, the dose
  ladder, and `fig7_ronchigram` (Ronchigram and wavefront as the aperture opens).
- `analysis/make_relaxation_figs.py` — figures A–E: focus fitted in the solver, focus fitted against
  known probe, phonons, non-round recall, and `figE_nonround_probe` (Ronchigram, wavefront and probe
  across the six-fold ladder).

Run one at a time with `--only`, e.g. `--only 7` or `--only E`. Both need
`~/hyperspy-bundle/bin/python`; the two Ronchigram figures build probes with abTEM and take about a
minute each.

Both write to `aberration_experiment/figs/<ISO-week>/meeting/`. **Every value is read from the run
output; nothing is typed into the scripts.** Keep it that way — it is why the pages and the data
cannot drift apart.

The user's standing instruction on figures: *no useless confusing detail.* One quantity per axis,
pick one probe-size measure and call it "probe diameter", and put a small "what is in this figure"
explanation next to each one rather than cramming it into the axes.

---

## The page — one artifact, not two

The user's instruction for this round is explicit: **one better page, replacing the two.** Today the
story is split across two links, which means the focus result is explained twice and the argument
never lands in one place.

| page | URL | covers today |
|---|---|---|
| Depth Past the Corrector | `https://claude.ai/artifact/NiYCfNo7uFAWFkyUZ5L3SF` | the baseline: what the aperture costs, depth against aperture, the focus search, dose |
| Removing the Six Assumptions | `https://claude.ai/artifact/3TfTwcipbbTxG9PGbeEtqi` | the 2026-09-21 results: focus end to end, the in-solver failure, phonons, non-round |

Build the single page **on the first URL**, so the link the user has already shared keeps working.
Read it with the Artifact tool first and publish your update to that same `url`. **Do not delete the
second page** without asking; leave it and tell the user plainly that it is superseded.

What the combined page must do that neither does now:

- **One spine.** What the method delivers, then each assumption removed in turn, ending on where it
  breaks and what the specification is. The ladder table is the backbone and every figure hangs off
  a rung of it.
- **Lead with the specification, not the chronology.** The publishable sentence has the shape *depth
  to X Å at Y mrad, provided the six-fold residual stays under Z waves and the dose above W*. X, Y
  and W are measured; **Z is what experiment C delivers**, so it is the headline of this round.
- **Say each thing once.** Both pages currently explain the focus result.
- **Every figure keeps its small "how to read this" note**, one quantity per axis, and no figure
  carries two scales.

Figures available, all regenerated from run output with nothing typed in:
`analysis/make_meeting_figs.py` gives 1–6 plus `fig7_ronchigram`; `analysis/make_relaxation_figs.py`
gives A–D plus `figE_nonround_probe`. The two Ronchigram figures are the strongest openers the study
has: one shows what opening the aperture costs, the other shows that six-fold astigmatism destroys
the recall while barely changing the probe diameter, which is the most quotable physical point in the
campaign.

**Embed figures as base64 data URIs inside the HTML.** Publishing them as separate artifact files and
referencing them by relative path does not render, even though the files list correctly — this cost a
round trip on 2026-09-21. Budget: page under 16 MB, base64 inflates by a third.

House style, which the user approved: dense prose in a reading column, numbered sections, wide figure
blocks with a small "how to read this" note, a status ledger with state chips, tables with the
baseline row tinted. Plain language, professor-level audience, no buzzwords, and only claims the user
could defend out loud.

## Where the ladder stands

From `aberration_experiment/results/relaxation_ladder.csv` (α 70 / α 90 where both exist):

| step | relaxation | outcome |
|---|---|---|
| 0 | known probe, no noise | depth error 0.56 / 0.37 Å; recall 98/82/75 and 95/94/95 % |
| 1 | focus unknown, corrector known | **no real cost**: 0.59 / 0.44 Å, recall 97/94/79 and 90/83/86 % |
| 1b | focus fitted inside the solver | **fails**; experiment A is the diagnosis, in hand and unanalysed |
| 2 | counting noise | works to 10⁵ e/Å²; at 10⁴ the matched kernel cannot be measured at all |
| 3 | thermal vibration | depth error roughly doubles; oxygen recall 10–17 %; species labels untrustworthy |
| 4 | non-round aberration | 0.6 waves → lead 39 %; 1.2 waves → 7 %. **Experiment C, 0.1–0.45 waves, is in hand and unanalysed — it is the headline number** |
| 5 | the 70 Å sample | diverged at the pinned step; the gentler retry is on the cluster, nothing down yet |
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
| six jobs pending forever on `DependencyNeverSatisfied` | a simulation whose output already existed exited 1 in zero seconds, so every reconstruction behind it could never start | `RECON_ONLY=1` when only the solver changed; the driver now refuses up front |
| `scp` finds nothing although the jobs finished | the tarball name carries the *submission* timestamp, so a batch sent at 22:56 is stamped the previous day | glob the stable part of the name, never the date |
| driver exits silently with no message | `VAR=$([ test ] && echo x)` fails the substitution when the test is false, and `set -e` kills the script | use an `if`, never that idiom |
| kernel quality numbers disagree with the pipeline's | `extract_psf` saves `exp(i·phase)`; take the phase first | see `kernel_quality()` in `make_meeting_figs.py` |

Two rules that are physics, not plumbing, and must hold on every leg: **`REGLAYER=0`** (it low-passes
the depth axis, which is the quantity being measured) and **one `BETA_LSQ` for every leg in a
comparison** — which is exactly why experiment B changes it deliberately and only for the thick slab,
and why its result must be compared against a thin-slab leg at the same step size before anything is
concluded about thickness.

---

## State of the repository

Two commits sit on `main` locally and are **not pushed**. The user's Mac needs `git push origin main`
(a previous attempt was blocked by a permission prompt), and Blythe then needs `git pull`:

| commit | why it matters |
|---|---|
| `e3f2572` campaign: refuse to submit a sim whose output already exists | stops the failure that killed experiment B |
| `d7bf927` figs: Ronchigram evolution for the round sweep and the six-fold ladder | `fig7_ronchigram` and `figE_nonround_probe`, both already published |

Ask for the push early. The experiment B rerun above works without it, because it uses a flag that
already exists on Blythe, but nothing else should be submitted until origin and Blythe agree.
