# Handoff — the non-round results were wrong; the engine is fixed and the campaign is re-running

Written 2026-09-22, after the round-two analysis. **Read `HANDOVER.md`** for the experiment itself,
then `NEXT_PHASE.md` for the plan and `PSF_KERNELS.md` for the kernel rules. The two published pages
are the fastest way in:

| page | URL | what it is |
|---|---|---|
| **Aperture Campaign Logbook** | `https://claude.ai/artifact/7ve93UM6yqCmiRcbfNuiJM` | **the one record**: the round-α campaign, every relaxation with its pictures, the aperture context, the rules. Point a new agent here first. |
| The Six-Fold Tolerance | `https://claude.ai/artifact/21fr4Y6JqCeaeB6KLMUKGz` | superseded — folded into the logbook's §5–6 on 2026-09-22; not maintained |

The logbook is rebuilt from `aberration_experiment/page/logbook.html` by `page/build_page.py` and republished
to the **same** URL. One page. Never start another for an update.

---

## How this project works

- **The user runs every cluster job.** No SSH from an agent session. Give exact `sbatch` blocks and
  wait for the log or the tarball; assume nothing about a job until you have seen it.
- **Every Blythe path lands under `$SHARE/phucrh`** (`/springbrook/share/physics/phucrh/`), never the
  group root. After each submission, `ls /springbrook/share/physics/` should show nothing new of ours.
- **Ship a matching `scp -O` line with every command block**, globbing the stable part of the tarball
  name (the stamp is the *submission* time, so a batch sent at 22:56 is stamped the previous day).
  `tar` on its own line, into a fresh directory.
- **Python is `~/hyperspy-bundle/bin/python`** for everything local.
- **Commit code, figures and results tables to `origin/main`.** Leave the user's report and the
  `analysis/atomfind/*.md` files alone; the `.py` files there are fair game.
- **Verify against the code or the data before asserting anything.** Every wrong number in this
  campaign came from an assumption that looked right — including one in the previous handoff, below.

---

## Where it stands

**Two results stand from the 2026-09-22 round, and one does not.**

**Valid — an unknown focus costs nothing.** An outer search over fixed trial probes, reading the solver's own
mismatch as the objective, puts its minimum on the true focus at both apertures (−2 ± 3 Å at 70 mrad, +2 ± 1 Å
at 90; the ±2 Å floor is the focus/depth degeneracy and is real). Reconstructing at the fitted focus costs
nothing measurable: depth error 0.56 → 0.59 Å at 70 mrad, 0.37 → 0.44 at 90. Fifty iterations rank the trials
as two hundred do, which is what makes it affordable. The tidier alternative — the solver fitting focus as its
one free probe parameter — is **closed as a failure**: at 90 mrad the focus never moves by 1 Å from starts 32 Å
out, at 70 mrad it moves the wrong way to about 21 Å below the truth whatever it started from, and neither
freezing the object nor capping the step changes it.

**Valid — shot noise behaves.** 10⁷ and 10⁶ e/Å² hold up, 10⁵ works with the light atoms going, and at 10⁴ no
single-atom reference can be measured from the data at all, so the finder never runs. Numbers in the ladder,
pictures on the page.

**VOID — everything non-round.** Every reconstruction with a non-round probe ever made in this campaign used a
probe rotated 30° from the one the data was made with (first row of *Traps*). That covers the August 2026
sweep, the 2026-09-21 escalation and the 2026-09-22 threshold ladder. A "six-fold tolerance below 0.1 waves"
was published on 2026-09-22 and withdrawn the same day. **The campaign currently has no measurement of how much
non-round aberration the method tolerates.** The engine is fixed, the simulations were never affected, and
every leg is re-reconstructing; see *Your job*.

**Still analytic, and unaffected by any of this**: how fast each residual grows with aperture. A term of order
n grows as α^(n+1), so with the residuals a real hexapole corrector leaves, every non-round term is over
0.1 waves by 70 mrad and the campaign's own six-fold is ten waves there. `campaign/aberration_waves.py` is the
one definition of "waves at the edge" and draws the figure. Note that the assumed instrument opened to 70 mrad
has a 25 Å probe against a 20 Å scan field — the geometry that ended the round sweep at 110 mrad — so it needs
no simulation to fail.

## Your job, in order

### 0/1. The campaign is on the cluster (submitted 2026-09-22, after the fix)

The queue was cleared and the whole non-round campaign resubmitted against the corrected engine, in four
blocks. The first is a single job with no dependency, so it reaches a GPU first and reports in about
fifteen minutes; **if it does not come back near 22, stop and cancel the rest** — the diagnosis is wrong
and nothing else is worth running.

```bash
cd /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline && git pull
grep -c "custom_data_flip applied to the PROBE" ptycho/+engines/+GPU/+initialize/load_from_p.m   # must print 1

# 1. confirmation: the 0.45-wave leg alone, existing sim, fixed engine. 71.2 before; ~22 = confirmed
TSV=campaign/nonround_sweep.tsv LABELS=nr0p45_C56_0p45w MODES=lab RECON_ONLY=1 bash campaign/run_thin_atomfind.sh

# 2. every already-simulated non-round leg, reconstructed again: 8 labels x 3 legs, ~10 min each
TSV=campaign/nonround_sweep.tsv RECON_ONLY=1 \
LABELS="nr0p1_C56_0p1w nr0p2_C56_0p2w nr0p3_C56_0p3w nr0p45_C56_0p45w nr1_C56_0p6w nr2_C56_1p2w nr3_C56_2p5w nr4_C56_C34" \
  bash campaign/run_thin_atomfind.sh

# 3. bracketing rungs + the other orders at 70 mrad: sim + recon, 9 labels x 3 legs.
#    OVERWRITE=1 because the cancelled batch may have left partial sims in those dirs.
TSV=campaign/nonround_sweep.tsv OVERWRITE=1 \
LABELS="nr0p02_C56_0p02w nr0p04_C56_0p04w nr0p07_C56_0p07w nrA1_C12_0p1w nrA1_C12_0p3w nrA2_C23_0p1w nrA2_C23_0p3w nrD4_C43_0p1w nrD4_C43_0p3w" \
  bash campaign/run_thin_atomfind.sh

# 4. six-fold at 90 mrad: sim + recon, 3 labels x 3 legs at BIN 2, ~1 h per recon
TSV=campaign/nonround_sweep.tsv OVERWRITE=1 LABELS="nr90_C56_0p1w nr90_C56_0p3w nr90_C56_10um" bash campaign/run_thin_atomfind.sh
```
```bash
mkdir -p ~/Desktop/nr_round3 && cd ~/Desktop/nr_round3
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_nr*_2026*.tgz' .
for t in *.tgz; do d="${t%.tgz}"; mkdir -p "$d"; tar xzf "$t" -C "$d"; done
tail -1 */recon_af_nr0p45_C56_0p45w_lab_NL14/analysis/*/*/*_error_trace.csv
```

About twenty GPU-hours over three cards. The 100 mrad rows (`nr100_C56_0p1w nr100_C56_0p3w nr100_C56_10um`)
are written but cost ~7.5 h per reconstruction; fire them only if 90 mrad is interesting.

**Analysing it.** Triage first (`analysis/triage_recon.py`), then `analysis/run_0922_atomfind.sh` is the
pattern — matched kernels from each leg's own Pb/Ti grids (`--zdrop 2` at NL14, `3` at NL23), then the blind
finder on its lab recon — then `analysis/relaxation_ladder.py --step 4 --alpha <A>` (it replaces rows on
`(step, label, alpha)`, so the void rows are overwritten in place), then `analysis/make_simple_figs.py`.
Delete the `VOID` prefix from a note only when that row has actually been remeasured.

**What the answers mean.** If 0.1 waves of two-fold costs far less than 0.1 waves of six-fold, the tolerance
depends on order and the specification must be written per term. If 90 mrad tolerates what 70 did not, it
loosens with aperture and the specification is per aperture. If the whole ladder now sits near the round
baseline, six-fold was never the problem and the earlier collapse was entirely the bug.

### 2. Read out experiment B when its tarball lands

It was running at the time of writing (`sacct` showed a70 lab and Pb complete, the rest queued).
Everything needed is built: `analysis/run_thin18_atomfind.sh <extracted dir>` triages, extracts the
matched kernels and runs the finder with the `thick` preset. Check first:

```bash
sacct -u phucrh --starttime 2026-09-22 --format=JobID,JobName%34,State,Elapsed | grep -E "thin18|af_pack"
```
```bash
mkdir -p ~/Desktop/thin18_0922 && cd ~/Desktop/thin18_0922
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_a70-90_thin18_20260922_*.tgz' .
for t in *.tgz; do d="${t%.tgz}"; mkdir -p "$d"; tar xzf "$t" -C "$d"; done
```
Note the `20260922` glob: the tarball already on disk stamped `20260921_123203` is the **diverged**
attempt, and it is 0 of 6 usable (see *Traps*).

**Its comparison must be against a thin slab at the same solver step**, or the result is about the
step size and not about thickness. That control has not been submitted yet:
`ALPHAS="70 90" RECON_ONLY=1 BETA_LSQ=0.02 bash campaign/run_thin_atomfind.sh` (it moves the existing
step-0 `analysis` dirs aside as `analysis.prev_<ts>`; the canonical 0.05 h5s are already on the Mac).

### 3. Then the remaining relaxations

Partial coherence (step 6) and specimen tilt (step 7) both need simulation code that does not exist;
`NEXT_PHASE.md` has the physics and the quantities. Then the combined run — every surviving
relaxation on together, across apertures — which is what makes the result publishable, and only then
the merge of the two older meeting pages onto the "Depth Past the Corrector" URL.

---

## Tools built this round

| what | where |
|---|---|
| triage that looks at the object, not just the log | `analysis/triage_recon.py` |
| waves ↔ C_nm for any order, any aperture | `campaign/aberration_waves.py` |
| depth constants derived from the box | `analysis/atomfind/config.derive_depth_constants` + the `thick` preset |
| dz and the depth constants from the simulation | `run_atomfind.py --sim-meta <sim_meta.mat>` |
| override any finder setting without editing a preset | `run_atomfind.py --set key=value` |
| the focus fit iteration by iteration | `make_relaxation_figs.py --only F` |
| 18-cell ground truth | `~/Desktop/thin18_gt/gt_prepared.npz` (19 139 atoms, z 4.000–73.935) |

---

## Traps

Each cost real time. They are fixed in the code; do not reintroduce them.

| symptom | cause | rule |
|---|---|---|
| **every non-round leg reconstructs badly, striped, stuck at a high mismatch, with the probe file exactly right** | `custom_data_flip = [0,0,1]` transposes the diffraction data on load and nothing transposed the probe with it. Both leave the sim in the same axis order. A round probe is its own transpose, so the orientation test (round probe) passed and every round leg was fine; a six-fold probe is rotated 30° (overlap with its transpose 0.95 / 0.82 / 0.69 / 0.57 at 0.1–0.45 waves), two-fold 90°, three-fold 30°; four-fold is the only term unaffected | fixed in `load_from_p.m`: the probe (and its aperture mask) now receive every flip the data does. **All non-round results before 2026-09-22 are void** and must be re-reconstructed `RECON_ONLY` — the simulations are fine |
| **an aberrated probe looks wrong in the h5** | `reconstruction/p/probe_initial` is stored **transposed** relative to `reconstruction/probes`, which is the one the solver used. A round probe is its own transpose so it never showed; a six-fold probe is not, and the untransposed comparison degrades with the aberration (0.95 / 0.82 / 0.69 / 0.57) | compare against `reconstruction/probes`, or transpose. Verified 2026-09-22: the probe the solver used matches the simulated aberrated probe with overlap **1.000000** on every leg |
| **a leg passes triage and is still junk** | it reported COMPLETED, wrote an h5 and logged no NaN, but its object is saturated: every layer at 1.5 rad phase std with 3 % of pixels wrapped. The 2026-09-21 thick batch was recorded as 5 of 6 failures; it was **6 of 6** | `analysis/triage_recon.py` — every genuine leg measured, including badly degraded ones, wraps **zero** pixels |
| a comparison figure's panels differ for no reason | one RNG shared across panels advanced on each call, so every Ronchigram panel drew a *different* amorphous film | `make_ronchigram_fig.ronchigram` seeds per call by default; pass a generator only for deliberately different films |
| figures missing from a published page | images published as separate files, referenced relatively | embed as base64 data URIs (`page/build_page.py`) |
| tarball downloads but the directory is empty | `tar` chained after `scp` on one line gets skipped | `tar` on its own line, into a fresh directory |
| every a90 result labelled a70 | aperture parsed from a tarball name containing `a70-90` | `relaxation_ladder.py --alpha` |
| six jobs pending forever on `DependencyNeverSatisfied` | a simulation whose output already existed exited 1 in zero seconds | `RECON_ONLY=1` when only the solver changed; the driver refuses up front |
| `scp` finds nothing although the jobs finished | the tarball name carries the *submission* timestamp | glob the stable part of the name, never the date |
| driver exits silently with no message | `VAR=$([ test ] && echo x)` fails the substitution when the test is false, and `set -e` kills the script | use an `if`, never that idiom |
| kernel quality numbers disagree with the pipeline's | `extract_psf` saves `exp(i·phase)`; take the phase first | `kernel_quality()` in `make_meeting_figs.py` |
| constants silently wrong for a new geometry | trims and caps typed for one slab thickness | `config.derive_depth_constants`; `--sim-meta` for dz |

Two rules that are physics, not plumbing, and must hold on every leg: **`REGLAYER=0`** (it low-passes
the depth axis, which is the quantity being measured) and **one `BETA_LSQ` for every leg in a
comparison** — which is why experiment B changes it deliberately and only for the thick slab, and why
its result must be compared against a thin-slab leg at the same step size.

---

## State of the repository

`origin/main` is current as of this handoff apart from the commits of this session, which are local
until pushed. Blythe pulled to `44f08f1` on 2026-09-22 and needs another `git pull` before the
submission in step 1 — the new sweep rows are in `campaign/nonround_sweep.tsv`.
