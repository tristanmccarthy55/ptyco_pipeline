# What's left — revised 2026-08-24 (report due Fri 4 Sep, ~11 days)

Supersedes the 2026-08-15 list. The change since then is that **the report exists** (v2:
9 pp main + 3 pp appendix + 2 pp refs), which converts most of the old "would strengthen it"
items into "post-submission / viva" and promotes a handful of correctness items to blocking.

---

## Objectives scorecard (against the project brief)

| Brief deliverable | Status |
|---|---|
| **D1** Initial MEP reconstruction of the PTO/STO model | **met** — Sec. III A, Appendix A |
| **D2** Atom finding + UQ → 3-D coords + Ti–O offsets mapping the vortices | **met** — Sec. III C/E; `atomfind/polarisation.py` closed the last gap on 15 Aug |
| **D3** CASTEP: how unit-cell displacements affect EELS | **substantially met** — M5 displacement scan is the core result; **apical oxygen only** (see B1 below) |
| **Stretch** Apply EELS models to the full PTO/STO model | not attempted; correctly framed as further work. Leave it |
| **Peer-to-peer reproducible result** | **described but not delivered** — Appendix E specifies it; nothing is packaged and it would not run on a peer's machine (see A1) |
| Report: 5–10 pp, journal format, reproducible result + its uncertainty | **met** |
| RRI/TRI: collaborator data attribution | met (acknowledgments) |

So: **D1 and D2 are done. D3 needs one scoping fix (free) and ideally one HPC run.
The peer-to-peer exercise is the only outright unmet requirement.**

---

## A. Blocking — must be done before submission, none needs compute

### A1. **Package and test the peer-to-peer exercise** — the one unmet requirement (half a day)
Appendix E documents commands that **would fail on a peer's machine today**.
`atomfind/config.py` hardcodes `~/Desktop/NL70_new_vol.npy`, `~/Desktop/psf_Pb_NL70_vol.npy`,
`~/Desktop/atomfind_out` and assumes the `~/hyperspy-bundle` interpreter. A peer has none of
these. Minimum viable fix:
1. Make `recon_vol`, `single_atom_vol`, `ti_kernel_vol` and `out_dir` overridable by CLI flag or
   environment variable, defaulting to a path *relative to the package*.
2. Ship a tarball: volume + Pb kernel + Ti kernel + reference structure + the two scripts +
   `environment.yml` (numpy, scipy, ase, matplotlib — note abtem/skimage are only needed by the
   deconvolution baselines, so a peer running only the finder needs less).
3. **Run it from a clean environment on a different machine.** Untested instructions are the
   normal failure mode for this exercise.
4. Include a one-page protocol sheet: inputs, the two commands, expected output with tolerances
   (Appendix E already has the numbers).

### A2. **Scope the 78% dichroism to the apical oxygen** — correctness (10 min)
`eels/RESULTS.md` M4 is explicit that the number is the **apical** `O:exc` K-edge, and that the
full O~K edge needs multiplicity weighting of 1×apical + 2×equatorial. The report currently
says "the O~K edge is $78\%$ dichroic" in the abstract, Sec. III F and the conclusions, without
that qualifier. As written it overstates the scope. Either run B1 and report the weighted edge,
or add "apical" and one sentence noting the weighting is outstanding. **Do the scoping fix now
regardless** — it is free and makes the claim true; B1 then upgrades it.

### A3. **Resolve the finite-dose forward reference** (10 min, or fold into B4)
Sec. III B says the discriminating comparison against peak-picking "is at finite dose
(Sec. III D)", and Sec. III D does not make that comparison. Either run B4 or soften the
sentence to state it as expectation rather than as a result delivered later.

### A4. **Fix polarisation panel (d)** (30 min)
Still outstanding from the old list, and you flagged it independently. Replace the paired bars
with **signal-to-noise = true spread / propagated σ on a log axis, line at 1**: δx 25, δy 8.7,
δz 0.31 — one bar below the line, no arithmetic asked of the reader. `make_pol_fig.py`.

### A5. **DOIs and bib hygiene** (1 h)
~8 DOIs in `report/refs.bib` were inherited unverified (the `VERIFY` markers were stripped to
make BibTeX parse, so the reminder is gone). Author/title/journal/volume/year are correct
throughout. Also: the two entries I corrected exist **only** in `report/refs.bib` — the segment
files still carry the originals (16 `VERIFY` markers across
`analysis/atomfind/paper/*.bib` and `eels/eels_refs.bib`). Sync them: Bugnet is PRB **93**,
020102(R) (**2016**); Mizoguchi is JPCM **21**, 104204 (**2009**).

---

## B. High value and genuinely achievable in 11 days

### B1. Equatorial oxygen → the full O~K edge (1 HPC job + analysis) — **gates A2**
`tet_Pz_Oeq` cells are built and waiting. This is the difference between "the apical O~K edge
is 78% dichroic" and "the O~K edge is X% dichroic", which is what D3 really asks for. Queue it
first; it has the longest lead time of anything here and it can change a headline number.

### B2. Rotational-invariance cross-check (1 cheap HPC job)
`tet_Pz` q⊥ must equal `tet_Px` q∥. The methods already present this as a validation; it is
currently unrun. Cheap, and it is what makes the 78% credible.

### B3. The `real` cell — what you would actually see in this material (1 HPC job)
`tet_Pz` is the full-along-beam upper bound. The `real` cell (82° median tilt, 0.046 Å along
beam) gives the honest number for the labyrinth in this zone axis. One SCF + OptaDOS pass, and
it closes the argument the conclusions currently leave at "upper bound".

### B4. **Finite-dose head-to-head — runs LOCALLY, no HPC** (a few hours)
I under-rated this on 15 Aug. The dose volumes are on the Desktop and `dose_series.py` runs on
the Mac, so this needs no queue. Run `peaks3d_raw` / `peaks3d_rl` / `peaks3d_mem` at each dose
alongside v3 (all four already exist in `run_atomfind.py`; they simply aren't invoked per-dose).
Expected payoff: the oxygen gap should *widen* with dose, converting "competitive on noiseless
data" into "the margin is where it matters", and it fulfils A3. If it does not widen, that is
also worth knowing and the claim must come out.

### B5. Re-run the dose series **and** the noise sweep under the current config (same session as B4)
Table II (dose) and Table D1 (noise) are both pre-deduplication-guard runs, each carrying a
caption caveat. One re-run with the current config removes both caveats and makes every table in
the report describe the same code.

---

## C. Post-submission / viva material — do not start these now

- **Kernel-transfer test** (localise NL70 with the NL105 dose kernel) — the "inverse crime"
  question a viva panel will ask. Worth having an answer ready even if unrun.
- **Blind cage-completeness detector** for the 5% wrong-composition tail.
- **Coherent NL105 kernel**; **in-situ vacancy-difference kernels**.
- **Energy-resolved M6** (fold the computed Δ(E) through the β integration).
- **Vector EELS figures** — regenerate from `eels/runs/exc/*.exc.txt` on Blythe in the paper's
  serif style. Cosmetic; the cropped PNGs are legible.
- Stale σ scale factors for guided/exit atoms (`RESULTS.md` §7).
- A recon-vs-truth figure (`fig8_nl70_potential.py`) if a reviewer asks how good the
  reconstruction actually is.

---

## Suggested order

1. **Today:** queue B1, B2, B3 on Blythe — longest lead time, and B1 gates A2.
2. **Today, while they queue:** A2 scoping fix, A4 panel (d), A5 bib sync. All free.
3. **This week:** A1, the peer package. Budget a full day; it will surface path bugs.
4. **This week:** B4 + B5 locally in one session; then finalise A3 either way.
5. **~1 Sep:** fold whatever HPC results landed into Sec. III F; re-verify data integrity and
   float placement; rebuild.
6. **Leave** everything in C for the viva.

If HPC time runs out, the report is still complete and honest with A1–A5 alone — B1–B3 upgrade
D3 from "substantially met" to "met", but A2's scoping fix is what keeps it truthful without
them.
