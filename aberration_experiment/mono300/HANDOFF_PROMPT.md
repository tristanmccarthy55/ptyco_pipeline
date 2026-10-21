# Prompt for the next session -- the 300 kV monochromator experiment (written 2026-10-21)

Paste the block below into a new session once the Blythe queue has emptied, followed by whatever the user has pasted
from Blythe since.

```text
You are continuing the 300 kV monochromated-microscope direction of the PX915 aberration campaign. Repo:
~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/ptychoshelves-clean (origin/main, push allowed).
Blythe checkout: /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline (updated by git pull). Results are pulled
to ~/Desktop/mono_phase1/{results,analysis}.

READ FIRST: memory "target-mono300"; aberration_experiment/mono300/PHASE1.md in full (the paper plan at the top, every
result table, every job id, the engine bug and the search-order lesson); the STATE block at the top of
aberration_experiment/HANDOFF_MONO300.md; aberration_experiment/page/MONO_PAGE.md (the results page, now this
direction's to keep up: version 8 published 2026-10-21, told as the story ARM fails -> coherence scan -> Spectra Mono -> physics put back, at https://claude.ai/artifact/B7ep8BpAMVdKvrnXywACq3).

THE PAPER (the user): short sections, one limiting factor each -- a bad (ARM) aberration set and opening the aperture;
partial coherence; cooling; then the chosen Spectra 300 X-FEG Mono (not UltiMono) and its outputs at doses x-y-z.

DONE AND SCORED: the coherence walls (focal spread good to 3.5 A, off from 5 A; source good to 0.1 A, okay 0.2, off 0.4),
the joint points (losses compound), A (ceiling fs 2.5 A + src 0.1 A with probe modes: good; + src 0.2: okay), B (A plus
phonons: fixed probe off at 100 K and 20 K; modes okay, Pb 51 / 70 %, Ti at guessing even at 20 K), C first try (failed
on the search ORDER: stage 1 searched C1 x C3 with A1 still ~40 A wrong, C1 ended 32.5 A off).

WHAT SHOULD HAVE FINISHED (set off 2026-10-20; confirm with sacct -j <ids> --format=JobID,JobName%30,State,Elapsed):
 - C rerun, 3-stage search (C1xC3 -> A1 -> C1xC3, stages 1-2 resumed): stage-3 array 1334455, select 1334456, refine
   legs 1334457/58/59, pack 1334460, analysis 1334461.
 - _mixF reruns of the four single wall points (phase1_step2_rerunF.sh, 8 modes at "4;2", 250 it, engine fix in):
   legs fs 7.29 1334466/67/68, fs 5 1334471/72/73, src 0.2 1334476/77/78, src 0.4 1334481/82/83; packs + analyses
   inferred as 1334469/70, 1334474/75, 1334479/80, 1334484/85 (confirm). At the 2026-10-21 check every leg seen was
   running and NaN-free; 1334478 was behind (64 iterations).
Pull lines for the user (run on the Mac):
 scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/analysis_armf_a080_fs2p5_src0p1_s3_ps_*.tgz' ~/Desktop/mono_phase1/analysis/
 scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_armf_a080_fs2p5_src0p1_s3_*.tgz' ~/Desktop/mono_phase1/results/
 scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/analysis_armf_a080_*_mixF_2026102*.tgz' ~/Desktop/mono_phase1/analysis/
 scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_armf_a080_*_mixF_2026102*.tgz' ~/Desktop/mono_phase1/results/
(A PSEARCH pack is named with RTAG but without _ps; the search folder search_af_armf_a080_lab_fs2p5_src0p1_s3 is inside
the atomfind_results pack.)

NEXT, IN ORDER:
 1. C: compare search_af_armf_a080_lab_fs2p5_src0p1_s3/search_best.json (and stage3's minimum) with start.json's
    truth record (C1 -958.6 A, C3 -406600 A, A1 12.4 A @ 1.468 rad); give each error in A and in waves at the aperture
    edge. Look at the depth sections of the _s3_ps legs before any score. Then score with
    aberration_experiment/mono300/phase1_look.py (add rows for fs2p5_src0p1_s3_ps).
 2. The _mixF single points: score them the same way (rows for the four _mixF points); the page shows the 7.29 A point
    from the earlier solver until they land (mono_results.json score_fallbacks).
 3. Update the page: fill the pending C and _mixF entries in aberration_experiment/page/mono_results.json and the §07 call
    and question/answer boxes in page/mono_2026-10-12.html (keep it readable cold: the PI reads it), rebuild with ~/hyperspy-bundle/bin/python analysis/make_mono_page.py, view the
    figures, Artifact read then publish with url= (no icon), copy the built page over
    PtycoShelves/artifacts/wcm_2026-10-12_coherence_wall.html, add a version entry to MONO_PAGE.md.
 4. If C works (search lands within the corrector's precision and Pb/Ti are at least okay), propose D -- the dose
    series 1e5-1e8 e/A^2 with the 3-stage search at the combined point (focal spread 2.5 A + source 0.1 A + phonons at
    20 K) -- with its cost (~250 GPU-h per dose at the current grids; smaller grids are an option) and build it only on
    the user's go. If C fails again, find out why from the stage tables before proposing anything.

HOW WE WORK: the user runs every cluster job -- give exact command blocks with their scp -O pull lines; never submit
anything yourself; dry-run every driver change (campaign/test_coh_names.sh, campaign/test_driver_legacy.sh, and a
stub-sbatch run of the submit script in a throwaway git copy); look at the phase images before judging; honest scores only
(make_arm_page.atom_scores: right element within 0.5 A in depth vs guessing; good = Pb >= 80 and Ti >= 45 %, okay = Pb
>= 50); plain language for a professor-level reader; stage files by name (other sessions have uncommitted work in
analysis/atomfind, eels, report); keep tests in the repo (macOS wipes /tmp).

TRAPS (each cost a round trip): the Mac has BSD tools and bash 3.2, Blythe GNU and bash 5 (GNU paste -s prints a blank
line on empty input; no associative arrays on the Mac); leak guards must test EXPORTED variables (printenv); the user
sometimes runs check lines from ~ (use absolute paths); mode legs need the all-modes illumination fix pulled; the object
step is still divided by the mode count, so mode legs can score a little below their fixed controls; 8 modes at "4;2"
~32 h (16;8 slower and NaNs more); search A1 before trusting C1; the seed-1 draw differs by ~1e-14 between the Mac and
Blythe (resume check uses 1e-6 relative); quota is 160 CPUs / 15 GPUs per user and --mem sets the CPUs.
```
