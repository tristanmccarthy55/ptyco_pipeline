# fusion — handover (2026-09-14)

> **CLOSED 2026-09-15.** The sign is in the hollow data only by model test; recovering it from a
> reconstruction failed (§2). The follow-on experiment (single-cell monolayer, ptychography for
> δxy + DFT-assisted EELS for |δz|) lives in its own repo: `../../monolayer-fusion/HANDOVER.md`.

Pairs with `README.md` (design + physics) and the `hollow-fusion-poc` auto-memory. This file is the
state of play: what is established, what is still open, where everything lives, and the traps.

Repo `ptychoshelves-clean` == `origin/main` (github `tristanmccarthy55/ptyco_pipeline`). Blythe:
`$SHARE/phucrh/ptyco_baseline/ptyco_pipeline`. **Runs are sbatch-only — hand the user commands.**
Analysis env `~/hyperspy-bundle/bin/python`. `fusion/runs/` is gitignored; Desktop copies under
`~/Desktop/fusion_recons/round{2,3,4,5}/` are the durable record.

---

## 1. What is established

**The sign is in the data — shown by a model test, not measured.** `sign_test.py` picks sign(P_z) in
**4 of 4 domains**, 100 % of
individual patterns, using only detector pixels **outside** the 0.75α hole and **no depth
reconstruction**. LLR ±0.997 for A/B (|δz| = 0.311 Å) against ±0.071 for C/D (|δz| = 0.166 Å).
The candidates come from the same simulator as the patterns, with no noise and probe, positions and
in-plane δ known, so this bounds what the data contain; it is not a measurement. Recovering the sign
from a reconstruction failed (§2).

**It holds at every hole size tested** (jobs 1274400–1274403; 200 patterns per domain): 4/4 and
100 % of patterns at 50, 75, 90 and 95 mrad. `sign_test` normalises each pattern, so LLR/pattern is
evidence per *recorded* electron (≈ the KL divergence); multiply by the measured fraction outside
the hole for evidence per *incident* electron:

| hole | to EELS | LLR/pattern A (C) | per incident e⁻ A (C) | e⁻/pattern for 100:1, A (C) |
|---|---|---|---|---|
| 50 mrad | 25.0 % | 0.00323 (0.000255) | 0.00242 (0.000191) | 1.9e3 (2.4e4) |
| 75 mrad | 56.2 % | 0.00498 (0.000355) | 0.00218 (0.000156) | 2.1e3 (3.0e4) |
| 90 mrad | 80.8 % | 0.0103 (0.000669) | 0.00197 (0.000128) | 2.3e3 (3.6e4) |
| 95 mrad | 90.0 % | 0.0186 (0.00117) | 0.00187 (0.000118) | 2.5e3 (3.9e4) |

Sending 25 → 90 % of the beam to the spectrometer costs only **23 %** of the sign evidence per
incident electron for A/B (39 % for C/D): the information sits just outside the bright-field edge.
The last column is ln(100)/KL — an expected-LLR statement, no Poisson fluctuation about it.

Fusion is what makes this a two-hypothesis problem: EELS fixes |δz|, projection fixes δxy, and only
the stacking order along the beam is left to decide.

**The hole costs the spectroscopy nothing.** O-K contrast is flat within 8 % from 0.10α to 1.00α —
beyond α ≈ 5 mrad the convergence cone alone sets the momentum-transfer spread. Size the hole from
the ptychography side. A near-parallel beam is *worse* (5.3 % vs 6.8 % limiting contrast): past the
magic angle the anisotropy inverts and saturates rather than vanishing.

**Measured dose split** (on the simulated patterns, not estimated): 0.50α → 24.3 % EELS, **0.75α →
55.5 / 44.5**, 0.95α → 89.8 %. Simultaneous HAADF takes 1.0 % of the beam and cannot separate the
domains at all (0.99 % spread) — as it must not. Do not show that image: the probe is 20 Å
overfocused (the ptychography geometry) and there are no phonons, so columns image as rings.

**Degeneracies verified, not assumed.** A≡B and C≡D in EELS to < 1e-12 (a round aperture cannot tell
an x-chain from a y-chain). The sign-encoding null — δ purely in-plane, so flipping δz is a no-op —
gives exactly 0.000e+00.

**The blind reconstruction does measure δxy.** `analyze_fusion.inplane_offsets`: Pb peaks in the
projected phase of the blind 0.15 Å object follow each domain's in-plane polar shift — RMS offset
**0.170 → 0.020 Å** once the predicted shift is subtracted (known start 0.175 → 0.022 Å). Ti–O
columns (Ti + apical O, |w_Ti| half of |w_Pb|) track less well blind: 0.076 → 0.069 Å (known start
0.061 → 0.029). `fuse()` still takes δxy as given; this is what would feed it.

---

## 2. The depth-reconstruction investigation (the long thread)

Blind multislice returned no lattice comb at 24, 20 and 16 layers, 200–600 iterations, on noiseless
data with the exact probe and positions. Ruled out in order:

| candidate | test | verdict |
|---|---|---|
| data don't carry depth | `depth_constraint.py` (true vs smeared stacking, identical projection) | **no** — 3.3–6.4 % pattern difference, 3→20 cells |
| layers / iterations | NL16, NL20 at 600 iters | **no** — still no comb |
| binned window truncation | probe intensity outside 11.7 Å window | **no** — 0.2–0.5 % |
| engine can't model the data | `known_object.py`, object frozen | **no** — truth scores **73.1** vs blind converged **158.8** |
| probe / beta_LSQ / bin / thickness | the user's NL70 labyrinth run works at 70 Å with betaLSQ 0.1, NpbstInf, p1, bin 4, reg 0, dz 0.999 | **no** — those are our settings |
| **scan redundancy** | re-sim at 0.15 Å step (25 600 positions) | **partly** — removes the Nyquist mode; atoms still not placed |

**Scan sampling explains the Nyquist artefact — not the whole failure.** At 0.3 Å step the depth
power sits at the layer-grid Nyquist — the two-slice odd/even mode
`aberration_experiment/HANDOVER.md` documents for a110/a120 — with 37–47 % of in-band power in the
top two kz bins. At **0.15 Å** step (matching NL70) the peak moves to the real lattice:

| run | kz peak | period | Nyquist | top-2-bin power |
|---|---|---|---|---|
| NL70 labyrinth (reference, works) | 0.2574 Å⁻¹ | 3.885 Å | — | — |
| fusion NL24, 0.3 Å step | 0.4730 | 2.11 Å | 0.4730 | **47 %** |
| fusion NL16, 0.3 Å step | 0.3154 | 3.17 Å | 0.3155 | — |
| fusion NL26, 0.15 Å step, presolve (step01) | 0.2365 | 4.23 Å | 0.5125 | 12 % |
| **fusion NL26, 0.15 Å step, full res (step02)** | **0.2365** | **4.23 Å (3.08×)** | 0.5125 | **22.7 %** |
| fusion NL16 from the true object | 0.2365 | 4.23 Å | — | — |

**The lattice comb is not the atoms.** On the full-resolution 0.15 Å object (`round5/`):

| object | atom-depth RMS, best global z shift | corr(⟨Pb⟩,⟨TiO⟩) | sign readout | residual @200 (slope) |
|---|---|---|---|---|
| blind, 0.3 Å, NL24 | 1.38 Å | −0.49 | 33 % | 159.2 (−0.44 %/it) |
| blind, 0.15 Å, NL26 step02 | **1.43 Å** | −0.28 | **10/25 cells** (flipped z: 64 %) | 159.0 (−0.41 %/it) |
| known start, 0.3 Å, NL16 | **0.58 Å** | −0.72 | 100 % (flipped z: 0 %) | 35.2 (−1.7 %/it) |

A uniform guess inside the ±0.45c fit window gives 1.08 Å, so the blind fits carry no depth
information. At the geometric depth origin, Ti selected by scattering weight (the atom-fit figure, since replaced — the RMS column above used `zs[::2]`, which also caught the cap-plane apical O): blind
0.3 Å 1.43 Å and 9/25 cells, blind 0.15 Å **1.54 Å** and 10/25, known start **0.54 Å** and 25/25.
Same processing, k_z prominence / top-two-bin power: 5.6× / 47.1 %, 3.1× / 22.7 %, 6.9× / 13.5 %;
NL70 4.4× / 4.7 % (3.9× without the per-layer median removal `load_volume` applies). Flipping z sends the known start to 0 %, so the depth
convention is right. Both blind runs are still descending at 200 iterations. Next tests, in order: continue
the 0.15 Å run to ≥1000 iterations (`PROBE_START`/restart from step02), then a known start on the
0.15 Å data to see whether that basin is reachable at all from this sampling.

**atomfind confirms it** (`atomfind_sign.py`, 2026-09-15). The v3 finder with the NL70 Pb/Ti
kernels (same probe, bin and engine settings; dz 0.999 vs 0.976), then each located Ti against its
equatorial O ring. Of the three references tried this is the reliable one: the full O6 cage needs
all six O (71 of 180 Ti on the ideal volume, 93 % right) and the apical pair is biased, because Ti
and apical O share a column 2 Å apart and blur together (domain A wrong even on the ideal volume).

| volume | z RMS Pb / Ti / O | Ti usable | per-Ti sign | domains |
|---|---|---|---|---|
| ideal (true atoms x kernels), noiseless | 0.04 / 0.17 / 0.24 Å | 180/180 | 100 % | 4/4 |
| ideal, 0.1x recon background | 0.05 / 0.17 / 0.24 Å | 180/180 | 99 % | 4/4 |
| ideal, 1x recon background (white) | 0.04 / 0.20 / 1.55 Å | 0 of 56 found | — | — |
| blind 0.15 Å step02 | 0.91 / 0.60 / 1.58 Å | 44/141 | 64 % | A +0.02±0.08, B −0.13±0.19, C −0.06±0.24 (wrong), D −0.41±0.10 |

**With atomfind's calibrated 95 % intervals** (split-conformal on each volume; z coverage 97 % blind,
96 % reference; one-to-one species-agnostic match, so blind depth RMS reads Pb/Ti/O 0.87/0.48/0.76 Å)
no blind Ti interval clears zero (0 of 44 usable; intervals ±1.3–1.6 Å) and every inverse-variance
domain mean spans zero: A −0.01±0.35, B −0.04±0.47, C −0.11±0.66, D −0.38±0.41 Å. Reference at 0.1x
background: 62 of 180 clear zero, all on the right side; means +0.31, −0.34, +0.19, −0.18 Å (±0.05).
Figures 4 and 5 of `make_figures.py` draw exactly this.

Atom positions carry the sign when the reconstruction has them at their depths; the blind object
does not. The control uses the finder's own kernels, so it is optimistic, and the 1x background
includes the blind object's depth smearing, so that row is a pessimistic bound. atomfind's Pb depth
RMS per checkpoint is the convergence criterion for the long run (ideal 0.04 Å; target <~ 0.2 Å).
Cached finder output and JSON sit in the session scratchpad; re-run with
`--ideal-noise 0 0.1 1 --cache <prefix> --reuse`.

The metric is validated: run on `~/Desktop/NL70_new_vol.npy` it reproduces that run's published
0.257 Å⁻¹ / 3.9 Å / 4.0× as **0.2574 / 3.885 / 3.92×**.

Notching the top two kz bins does **not** rescue the 0.3 Å runs (33→33 %, 25→28 %), so the Nyquist
mode is the solver filling an unconstrained null space, not a veil over a good solution.

**Known-object diagnostics** (these are diagnostics, NOT measurements — started from the truth their
signs are not independent): frozen truth residual **73.1** vs blind converged **158.8**; released it
reaches **35.2** and keeps the comb (3.2–5.3 vs 0.4–0.7 blind); the matched-filter readout then gets
**100 % of cells in all four domains**. So the readout was never the problem.

---

## 3. Open

1. **Blind depth at 0.15 Å step does not place atoms** (§2): comb present, atoms at RMS 1.54 Å,
   sign readout 10/25 cells, residual still falling at 200 iterations. Run it longer; then a known start
   on the 0.15 Å data. **Beware `ls -1v`** on these folder names: `g32` (step02) sorts before `g64`
   (step01), so `| tail -1` yields the *presolve* file. Select `*step02*` explicitly.
2. **Leg A (marker specimen) was simulated but never reconstructed.** Both NL26 jobs loaded 25 600
   positions, i.e. both ran leg B — the recon folder tag did not include the data source, so the
   second launch re-pointed the first's symlinks. Fixed in `c25719b`; the sim
   (`fusion/runs/fusion_mark`, 288 Sr markers) is ready to reconstruct. Now likely unnecessary
   since scan redundancy explains the failure, but it is the clean control for the
   z-degeneracy hypothesis.
3. Frozen phonons and finite dose are both off. `--phonons` and `sim/add_poisson_noise.py` exist.
4. Per-column (rather than per-domain) sign, which needs the depth route working blind.
5. ~~HSA sweep of the sign test~~ — done, §1: 4/4 at every hole 50–95 mrad.
6. The sign test uses 200 patterns per domain at infinite dose. A Poisson draw at the
   e⁻/pattern of §1 would turn the expected-LLR table into an error rate.

---

## 4. Files

| file | role |
|---|---|
| `toy_sample.py` | the four-domain membrane + ground truth, with build gates; `--marker-z` breaks z-periodicity |
| `eels_forward.py` | per-column O-K spectra: chain-resolved tensor × convergent-aperture average × the M5 ladder |
| `simulate_fusion.py` | one abTEM 4D-STEM scan → PtychoShelves inputs, hollow masks, measured dose split, virtual HAADF/BF |
| `sign_test.py` | **the headline**: sign(P_z) from the hollow data by two-hypothesis likelihood ratio |
| `sign_encoding.py` | is the sign in the data, and how does it grow with thickness (noiseless) |
| `depth_constraint.py` | do the data constrain depth at all (true vs smeared stacking) |
| `known_object.py` | writes the TRUE object in the engine's frame; optimiser-vs-model diagnosis |
| `check_recon.py` | **run before reading physics off any reconstruction** — did it recover depth? |
| `atomfind_sign.py` | sign from located atoms: atomfind v3 on an ideal control and a reconstruction, Ti against its O |
| `analyze_fusion.py` | matched-filter depth readout, EELS magnitude inversion, fusion |
| `make_figures.py` | figures 1–5 (headline, EELS axis, model test vs hole, depth sections with atomfind picks, Ti–equatorial-O offsets) |
| `../ptycho/run_fusion_hollow.m` | MHP driver: `run_synthetic_recon_ML.m` + the hole, `KNOWN_OBJECT`, `OBJECT_START` |
| `test_fusion.py` | 19 tests, includes the labyrinth-code-untouched gate |

---

## 5. Running it

```bash
# specimen + local checks
~/hyperspy-bundle/bin/python toy_sample.py --preview
~/hyperspy-bundle/bin/python eels_forward.py --cache        # once; needs eels/runs/*_core_edge.dat
~/hyperspy-bundle/bin/python test_fusion.py

# Blythe: simulate (GPU, sbatch only)
SCAN_STEP=0.15 sbatch --output=fsim_%j.out --error=fsim_%j.err fusion/run_fusion_sim.slurm
# hollow reconstructions -- one sim serves every hole
bash fusion/run_fusion_recon.sh 0 0.50 0.75 0.95
# the blind sign test, and any other GPU experiment
bash fusion/run_gpu.sh sign_test
bash fusion/run_gpu.sh depth_constraint --thickness 3 5 10 20

# read out
~/hyperspy-bundle/bin/python check_recon.py --recon <...>/*step02*/Niter*.mat
~/hyperspy-bundle/bin/python make_figures.py              # figures 1-5, inputs under ~/Desktop/fusion_recons
```

`SCAN_STEP=0.15` removes the Nyquist mode; neither step yet gives atom depths blind (§2).
Figures: `~/hyperspy-bundle/bin/python make_figures.py` → `figs/fig1`–`fig5`.

---

## 6. Traps (all cost real time; most are now guarded by tests)

- `fourier_error_out` scores **every** detector pixel while the hollow mask is applied only inside
  `modulus_constraint.m`. It RISES on a hollow run and is **not** a valid convergence monitor there.
- `save_outputs` writes `dp = A.transpose(0,2,1)` for MATLAB; read back in Python it is the
  transpose of the simulated pattern. The radially symmetric mask hides it — symptom is a quiet 50 %.
- Object pixel size is `box/(Ndpx·bin)`, **not** window/Nx (1 % off = 0.3 Å drift across the field).
  The presolve engine has its own, coarser pixel size — read `outputs.pixel_size`.
- Comparing domains **within** one scan is dominated by intra-cell probe-position sampling, not
  structure. Compare uniform specimens at matched positions (`sign_encoding.py`).
- A GPU run once returned bit-identical patterns for two different structures (exact 0.0) where CPU
  gave 2.96 %. `assert_distinct()` now makes that fatal.
- The Blythe login node has no CUDA driver; `require_gpu_job()` refuses GPU work outside SLURM.
- `ls -1v` on the recon folder names picks the presolve engine — see §3.1.
- Recon folders must carry their data source in the tag, or two legs collide (fixed, `c25719b`).

**Hard constraint, tested:** `sim/simulate_4dstem.py` and `ptycho/run_synthetic_recon_ML.m` carry the
validated PTO/STO labyrinth geometry and are imported or copied, never edited.
`test_labyrinth_code_untouched` git-diffs them as a gate. Keep it that way.

---

## 7. Caveats for any write-up

- The EELS channel is a **kinematic** forward model: first-principles unit-cell spectra
  (CASTEP + OptaDOS core hole) combined by the validated uniaxial tensor form, averaged over the α
  and β cones. No channelling, no thickness dependence.
- All data are **noiseless** (1e10 e/pattern, no Poisson) and **coherent** (no phonons).
- The 4.7 % EELS contrast is a **domain-scale** number at realistic dose, not per probe position.
- Known-object runs are diagnostics; the blind result is the sign test.
