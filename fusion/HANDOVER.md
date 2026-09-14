# fusion — handover (2026-09-14)

Pairs with `README.md` (design + physics) and the `hollow-fusion-poc` auto-memory. This file is the
state of play: what is established, what is still open, where everything lives, and the traps.

Repo `ptychoshelves-clean` == `origin/main` (github `tristanmccarthy55/ptyco_pipeline`). Blythe:
`$SHARE/phucrh/ptyco_baseline/ptyco_pipeline`. **Runs are sbatch-only — hand the user commands.**
Analysis env `~/hyperspy-bundle/bin/python`. `fusion/runs/` is gitignored; Desktop copies under
`~/Desktop/fusion_recons/round{2,3,4}/` are the durable record.

---

## 1. What is established

**The proof of concept works.** `sign_test.py` recovers sign(P_z) in **4 of 4 domains**, 100 % of
individual patterns, using only detector pixels **outside** the 0.75α hole and **no depth
reconstruction**. LLR ±0.997 for A/B (|δz| = 0.311 Å) against ±0.071 for C/D (|δz| = 0.166 Å).
Fused δ matches ground truth exactly in all four domains.

Fusion is what makes this a two-hypothesis problem: EELS fixes |δz|, projection fixes δxy, and only
the stacking order along the beam is left to decide.

**The hole costs the spectroscopy nothing.** O-K contrast is flat within 8 % from 0.10α to 1.00α —
beyond α ≈ 5 mrad the convergence cone alone sets the momentum-transfer spread. Size the hole from
the ptychography side. A near-parallel beam is *worse* (5.3 % vs 6.8 % limiting contrast): past the
magic angle the anisotropy inverts and saturates rather than vanishing.

**Measured dose split** (on the simulated patterns, not estimated): 0.50α → 24.3 % EELS, **0.75α →
55.5 / 44.5**, 0.95α → 89.8 %. Simultaneous HAADF takes 1.0 % of the beam and cannot separate the
domains at all (0.99 % spread) — as it must not.

**Degeneracies verified, not assumed.** A≡B and C≡D in EELS to < 1e-12 (a round aperture cannot tell
an x-chain from a y-chain). The sign-encoding null — δ purely in-plane, so flipping δz is a no-op —
gives exactly 0.000e+00.

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
| **scan redundancy** | re-sim at 0.15 Å step (25 600 positions) | **YES** |

**The answer is scan sampling.** At 0.3 Å step the depth power sits at the layer-grid Nyquist — the
two-slice odd/even mode `aberration_experiment/HANDOVER.md` documents for a110/a120 — with 37–47 %
of in-band power in the top two kz bins. At **0.15 Å** step (matching NL70) the peak moves to the
real lattice:

| run | kz peak | period | Nyquist | top-2-bin power |
|---|---|---|---|---|
| NL70 labyrinth (reference, works) | 0.2574 Å⁻¹ | 3.885 Å | — | — |
| fusion NL24, 0.3 Å step | 0.4730 | 2.11 Å | 0.4730 | **47 %** |
| fusion NL16, 0.3 Å step | 0.3154 | 3.17 Å | 0.3155 | — |
| **fusion NL26, 0.15 Å step** | **0.2365** | **4.23 Å** | 0.5125 | **12 %** |
| fusion NL16 from the true object | 0.2365 | 4.23 Å | — | — |

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

1. **The 0.15 Å full-resolution object has not been analysed.** Only the *presolve* engine
   (`*_Ndp118_step01*`, 0.099 Å pixels) was pulled. On it the lattice is recovered but the sign
   readout gives 41 % with a common-mode positive bias on all four domains — expected, since the
   sign lives in sub-Å offsets between the Pb and Ti–O sublattices. Pull
   `recon_hsa0_NL26/01/*step02*/Niter200.mat` and re-run `check_recon.py` + the readout.
   **Beware `ls -1v`** on these folder names: `g32` (step02) sorts before `g64` (step01), so
   `| tail -1` yields the *presolve* file. Select `*step02*` explicitly.
2. **Leg A (marker specimen) was simulated but never reconstructed.** Both NL26 jobs loaded 25 600
   positions, i.e. both ran leg B — the recon folder tag did not include the data source, so the
   second launch re-pointed the first's symlinks. Fixed in `c25719b`; the sim
   (`fusion/runs/fusion_mark`, 288 Sr markers) is ready to reconstruct. Now likely unnecessary
   since scan redundancy explains the failure, but it is the clean control for the
   z-degeneracy hypothesis.
3. Frozen phonons and finite dose are both off. `--phonons` and `sim/add_poisson_noise.py` exist.
4. Per-column (rather than per-domain) sign, which needs the depth route working blind.
5. HSA sweep of the sign test itself (jobs 1274400–1274403) — results not yet read.

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
| `analyze_fusion.py` | matched-filter depth readout, EELS magnitude inversion, fusion |
| `make_figure.py` | the headline figure |
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
~/hyperspy-bundle/bin/python make_figure.py --budget runs/fusion/hollow_budget.json
```

`SCAN_STEP=0.15` is now the default worth using — 0.3 does not give depth.

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
