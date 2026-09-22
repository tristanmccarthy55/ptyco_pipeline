# Handoff — the tolerance is below the ladder; bracket it, then finish the relaxations

Written 2026-09-22, after the round-two analysis. **Read `HANDOVER.md`** for the experiment itself,
then `NEXT_PHASE.md` for the plan and `PSF_KERNELS.md` for the kernel rules. The two published pages
are the fastest way in:

| page | URL | what it is |
|---|---|---|
| **Aperture Campaign Logbook** | `https://claude.ai/artifact/7ve93UM6yqCmiRcbfNuiJM` | the permanent record: the round-α campaign, the rules, the conventions. **Point a new agent here first.** |
| **The Six-Fold Tolerance** | `https://claude.ai/artifact/21fr4Y6JqCeaeB6KLMUKGz` | round two (2026-09-22): the specification, the finder-threshold control, the close of the in-solver focus fit |

Both are rebuilt from `aberration_experiment/page/*.html` by `page/build_page.py` and republished to
the **same** URLs. Never start a new page for an update.

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

**The deliverable is a specification, and it is tighter than the ladder that measured it.** At
70 mrad, six-fold astigmatism must stay below **0.1 waves** at the aperture edge (C₅₆ < 10 µm) — and
that is an upper bound, not the tolerance, because 0.1 waves is the smallest residual tested and it
already fails: lead recall 98 → 80 %, oxygen 75 → 48 %, depth error 0.56 → 0.76 Å, species confusion
1.3 → 12.4 % against a 5 % health threshold. A real hexapole tableau leaves C₅₆ of order 1 mm, about
ten waves, so the requirement is ~100× better than as-left, and it tightens with aperture (0.1 waves
is C₅₆ < 2.2 µm at 90 mrad).

**Closed this round.** Fitting focus inside the solver: at 90 mrad the focus does not move by 1 Å from
starts up to 32 Å out, object free or frozen; at 70 mrad it moves but never toward the truth; a
0.05 Å/iteration cap does not stop the drift. The outer search over trial probes is the permanent
answer. Say so plainly rather than leaving it open.

**Answered this round.** Is the recall collapse the finder's thresholds? For non-round, partly — at
0.1 waves a relaxed quality cut brings lead back 80 → 89 % and oxygen 48 → 69 %, while the round-probe
control barely moves, so the recovery is real. For phonons, not at all: the identical change leaves
those numbers untouched to two decimal places. Tables in `results/2026-W39/`.

### Ladder — `aberration_experiment/results/relaxation_ladder.csv`, 20 rows

| step | relaxation | outcome |
|---|---|---|
| 0 | known probe, no noise | depth error 0.56 / 0.37 Å; recall 98/82/75 and 95/94/95 % |
| 1 | focus unknown, corrector known | **no real cost**, via an outer search over trial probes |
| 1b | focus fitted inside the solver | **fails, and now diagnosed** — closed |
| 2 | counting noise | works to 10⁵ e/Å²; at 10⁴ the matched kernel cannot be measured |
| 3 | thermal vibration | depth error roughly doubles; **not a finder artefact** |
| 4 | non-round (six-fold) | **tolerance below 0.1 waves** — the smallest tested, and it fails |
| 4b | non-round, other orders | rows written, **not yet run** — the next submission |
| 5 | the 70 Å sample | running on the cluster now; analysis path built and waiting |
| 6 | partial coherence | **not built** |
| 7 | specimen tilt | **not built** |

---

## Your job, in order

### 1. Submit the next sweep (rows are written and validated; nothing else is blocked on it)

Nine labels at 70 mrad, BIN 4, 27 reconstructions plus a pack job — about 10 minutes each on three
GPUs. Three rungs bracket six-fold from beneath (0.02, 0.04, 0.07 waves); six put one term per order
a hexapole actually leaves (two-fold C12, three-fold C23, three-lobe C43) at 0.1 and 0.3 waves, the
levels where six-fold is already measured, so the comparison is direct.

```bash
cd /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline && git pull
TSV=campaign/nonround_sweep.tsv \
LABELS="nr0p02_C56_0p02w nr0p04_C56_0p04w nr0p07_C56_0p07w nrA1_C12_0p1w nrA1_C12_0p3w nrA2_C23_0p1w nrA2_C23_0p3w nrD4_C43_0p1w nrD4_C43_0p3w" \
  bash campaign/run_thin_atomfind.sh
ls /springbrook/share/physics/        # nothing new of ours in the group root
```
```bash
mkdir -p ~/Desktop/nr_round3 && cd ~/Desktop/nr_round3
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_nr*_2026*.tgz' .
for t in *.tgz; do d="${t%.tgz}"; mkdir -p "$d"; tar xzf "$t" -C "$d"; done
```

Analyse exactly as round two did: `analysis/run_0922_atomfind.sh` is the pattern (kernels from each
rung's own grids, then the blind finder on its lab recon), then `analysis/relaxation_ladder.py
--step 4 --alpha 70`. The figure axis and the coefficients both come from
`campaign/aberration_waves.py`, so nothing is typed.

**What the answer looks like.** If 0.1 waves of two-fold costs far less than 0.1 waves of six-fold,
the tolerance depends on order and the specification must be written per term — expected, since for
the same edge phase a term of order n throws the edge ray by (n+1)λ/θ, 0.56 Å at n = 1 against 1.69 Å
at n = 5. If they all cost the same, one number in waves covers the whole tableau, which is the more
useful result and the easier one to publish.

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
