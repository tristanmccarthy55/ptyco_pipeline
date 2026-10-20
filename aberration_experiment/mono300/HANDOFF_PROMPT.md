# Prompt for the next session -- the 300 kV monochromator experiment (written 2026-10-20)

Paste the block below into a new session, followed by whatever the user has pasted from Blythe since.

```text
You are continuing the 300 kV monochromated-microscope direction of the PX915 aberration campaign. Repo:
~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/ptychoshelves-clean (origin/main, push allowed).
Blythe checkout: /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline (updated by git pull). Results are pulled
to ~/Desktop/mono_phase1/{results,analysis}.

READ FIRST: memory "target-mono300"; aberration_experiment/mono300/PHASE1.md in full (the paper plan at the top, every
result table, every job id, the engine bug and the search-order lesson); the STATE block at the top of
aberration_experiment/HANDOFF_MONO300.md. If the results page is to be updated, page/MONO_PAGE.md (another session owns it).

WHERE THINGS ARE: running on Blythe -- the C rerun (searched start: C1/C3/A1 drawn from the corrector's precision, a
3-stage outer search C1xC3 -> A1 -> C1xC3, refine with probe modes from the mixed state on the searched probe; ids in
PHASE1.md) and the _mixF reruns of the four single wall points (phase1_step2_rerunF.sh; record their ids from the user's
paste into PHASE1.md). Done and scored: the coherence walls, the joint points, A (ceiling with modes) and B (ceiling +
phonons at 100 K and 20 K: Pb okay, Ti at guessing).

NEXT: (1) when the C rerun lands, compare search_af_armf_a080_lab_fs2p5_src0p1_s3/search_best.json with start.json's
truth record, look at the depth sections, then score with aberration_experiment/mono300/phase1_look.py (add rows for
fs2p5_src0p1_s3_ps and the _mixF single points). (2) If C works, propose D -- the dose series 1e5-1e8 with the 3-stage
search at the combined point (focal spread 2.5 A + source 0.1 A + phonons at 20 K) -- with its cost (~250 GPU-h per dose
at the current grids; smaller grids are an option), and build it only on the user's go.

HOW WE WORK: the user runs every cluster job -- give exact command blocks with their scp -O pull lines; dry-run every
driver change (campaign/test_coh_names.sh, campaign/test_driver_legacy.sh, and a stub-sbatch run of the submit script in a
throwaway git copy); look at the phase images before judging; honest scores only (make_arm_page.atom_scores: right element
within 0.5 A in depth vs guessing; good = Pb >= 80 and Ti >= 45 %); plain language for a professor-level reader; stage files
by name (other sessions have uncommitted work); keep tests in the repo (macOS wipes /tmp).

TRAPS (each cost a round trip): the Mac has BSD tools and bash 3.2, Blythe GNU and bash 5 (GNU paste -s prints a blank
line on empty input; no associative arrays); leak guards must test EXPORTED variables (printenv); the user sometimes runs
check lines from ~ (use absolute paths); mode legs need the all-modes illumination fix pulled; 8 modes at "4;2" ~32 h
(16;8 slower and NaNs more); search A1 before trusting C1; a PSEARCH pack is named with RTAG but without _ps.
```
