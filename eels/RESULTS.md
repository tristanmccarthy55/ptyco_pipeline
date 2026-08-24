# RESULTS — PbTiO₃ ELNES polarisation study

Milestone-by-milestone evidence. Append to this as HPC gates clear. (See README.md for the
full staircase and pass criteria.)

---

## M1 — structure build + validation — ✅ PASSED

`~/hyperspy-bundle/bin/python build_cells.py --corehole` (beam = +z). All 9 cells passed
every crystallographic gate; 16 base `.cell` + 35 core-hole `.cell` written.

| cell | space group | c/a | \|δ_Ti\| (Å) | along-beam / perp | \|δ_Pb\| (Å) |
|------|-------------|-----|-------------|-------------------|-------------|
| cubic | Pm-3m | 1.0000 | 0.000 | 0.000 / 0.000 | 0.000 |
| tet_Pz | P4mm | 1.0635 | 0.331 | **0.331 / 0.000** | 0.487 |
| tet_Px | P4mm | 1.0635 | 0.331 | 0.000 / 0.331 | 0.487 |
| real | P4mm | 1.0635 | 0.331 | 0.046 / 0.328 | 0.487 |
| scan_0.00 | P4/mmm | 1.0635 | 0.000 | — | 0.000 |
| scan_0.25→1.00 | P4mm | 1.0635 | 0.083→0.331 | along z | 0.122→0.487 |

- c/a = 1.0635 (target 1.0635 ✓); \|δ_Ti\| = 0.331 Å (target ~0.30, tol 0.08 ✓); \|δ_Pb\| =
  0.487 Å (lit. ~0.47 ✓).
- Orientation correct: `tet_Pz` puts **all** of δ_Ti along the beam, `tet_Px` **all**
  perpendicular → the two limiting cases for the M4 dichroism.
- Scan ladder along-beam δ_Ti = [0, 0.083, 0.165, 0.248, 0.331] Å, **linear r = 1.00000**,
  monotonic → clean M5 calibration axis.
- Core-hole cells: exactly one atom per file relabelled `X:exc`, chosen nearest the 2×2×2
  centre (max separation from its periodic images).

### Early scientific finding (from the real labyrinth, sim orientation)
Over **2560 / 3888** polar Ti (\|δ\|>0.15 Å) in `PTO6_STO6_18_18_labyrinthPoscar.vasp`:
median \|δ\| 0.227 Å (max 0.397); **angle-from-beam median 82°** (IQR 73–87°) → **median
along-beam fraction only ~0.14**. Viewed side-on (the pipeline's zone axis) the vortex
polarisation is *mostly in-plane*; the along-beam component we're chasing is intrinsically
small here. `real` cell built at the median 82° tilt captures this (0.046 Å along beam).
**Implication:** `tet_Pz` (full along-beam) is the upper bound; expect the real measurable
signal to be a small fraction of it. A different zone axis would raise the along-beam share.

---

## M6 — detectability geometry model — ✅ VALIDATED (local, synthetic)

`~/hyperspy-bundle/bin/python analyze_elnes.py --selftest`. The collection-aperture
averaging model is validated before any real spectrum:

- **θ_E(O K, 300 keV) = 1.09 mrad**; **magic angle = 4.32 mrad = 3.98 θ_E** (textbook ≈3.97 ✓).
- Surviving intrinsic dichroism vs collection semi-angle β:

  | β (mrad) | 1 | 2 | 5 | 10 | 25 | 100 |
  |---|---|---|---|---|---|---|
  | fraction surviving | 0.62 | 0.28 | −0.04 | −0.17 | −0.26 | **−0.33** |

- **Your 100 mrad ptychography probe → −0.33**: the anisotropy is averaged away / inverted.
  To measure along-beam polarisation, use a *small* EELS collection aperture, ideally
  ≲ 2 mrad (≪ the 4.3 mrad magic angle).
- SNR: a 2% fractional anisotropy needs ~4.5×10⁴ counts/channel at SNR 3 — sets the dose target.

---

## M6c — CONVERGENCE, not just collection — ⚠️ corrects the M6 recommendation (2026-08-24)

`analyze_elnes.py --convergence`. M6 above (`surviving_anisotropy`) is the textbook
**parallel-illumination** result: one incident direction, collect out to β. That is the only
geometry for which a magic angle is defined, and it was the whole basis for "use a small EELS
aperture". A focused STEM probe is not that geometry: **q_t = k(θ_f − θ_i)** is a *vector*
difference, so the convergence cone spreads **q** on its own however small β is.

Averaging the dipole weights over both cones (`surviving_anisotropy_convergent`, validated
against the analytic model to <1e-3 as α→0):

| surviving fraction | β=1 | β=2 | β=5 | β=10 | β=25 | β=100 mrad |
|---|---|---|---|---|---|---|
| **α = 0 (parallel)** | **0.62** | 0.28 | −0.04 | −0.17 | −0.26 | −0.33 |
| α = 1 | 0.52 | 0.29 | −0.04 | −0.17 | −0.26 | −0.34 |
| α = 2 | 0.29 | 0.25 | −0.03 | −0.16 | −0.26 | −0.34 |
| α = 5 | −0.04 | −0.03 | −0.01 | −0.16 | −0.26 | −0.34 |
| α = 20 | −0.24 | −0.25 | −0.24 | −0.24 | −0.24 | −0.32 |

**Closing the spectrometer aperture does not recover the signal.** Beyond α ≈ 5 mrad *no* β
does — the best achievable is −0.01. The convergence must itself come down to the order of
θ_E = 1.09 mrad, i.e. near-parallel illumination.

**This strengthens the instrumental conclusion.** The incompatibility between depth-resolved
imaging and polarisation spectroscopy is not a shared-aperture problem that a second aperture
would solve; it is an *illumination-mode* problem. Atomic-resolution depth sectioning needs a
strongly convergent probe; the dichroism needs a nearly parallel one.

**What was and was not accounted for before this.** `config.STEMEELS` (the forward simulator)
always separated α = 20 mrad from β = 100 mrad and says so explicitly. `config.Optics` (the M6
analysis) carries only `convergence_mrad = 100`, and `geometry_report` fed that number into the
β-model and reported it as the probe's surviving anisotropy. For the combined instrument the two
coincide (β = detector_max/2 = 100 mrad), so −0.33 was numerically right for that configuration,
but right by coincidence rather than because convergence was modelled. Both are now computed
properly and `geometry_report` states which angle it is using.

*Caveat:* kinematic dipole treatment with an **incoherent** average over incident directions. It
ignores coherence between incident directions and channelling, both of which the multislice route
(M6b) carries. Strictly more complete than the α = 0 model it replaces, not a final word.

---

## M6b — dynamical forward model (abtem multislice EELS) — ⚙️ API validated, EELS backend gated on gpaw

`simulate_eels.py`. Confirmed the exact abtem 1.0.5 core-loss call path works end-to-end
(`SubshellTransitions(Z,n,l).get_transition_potentials` → `Probe.transition_potential_scan`
→ detectors), on the same abtem engine as `sim/simulate_4dstem.py`.
- **Elastic multislice path validated without gpaw** (`--selftest` passes): build cells →
  potential → probe → scan produces finite channelling output.
- **Core-loss EELS needs `gpaw`**: abtem 1.0.5's atomic radial wavefunctions are gpaw-only
  (no hydrogenic fallback). gpaw is the *accurate* (all-electron) backend, so it fits the
  accuracy goal — but it's a real dependency to add to the bundle:
  `~/hyperspy-bundle/bin/pip install gpaw && ~/hyperspy-bundle/bin/gpaw install-data ...`
  (macOS may need libxc + a C toolchain). Once present, `--edge O_K --thickness-series`
  gives the Channel-B along-beam contrast (tet_Pz vs tet_Px) vs thickness.
- **Route comparison (the ask):** abtem (isolated-atom potentials) = Channel B only; OptaDOS =
  Channel A only → comparing them shows which dominates the measurable along-beam signal.
  BiP-PRISM/scatterem is an *approximation to* exact abtem multislice, so on unit cells exact
  abtem is the accuracy reference and BiP-PRISM adds nothing until scale.

## Pre-flight (CASTEP-independent) — ✅ ALL GREEN (Mac + Blythe HPC)

Also validated **on Blythe** (SLURM job 1240612, COMPLETED/0:0, in the `abtem` env, Python 3.11,
NumPy 2.x, no spglib): 10/10 tests, ALL PASSED build, and the abtem-1.0.9 elastic multislice
self-test passed unchanged — HPC env + code + `.cell` generation + the channelling engine all
confirmed. (`run_pyeels.slurm`; graceful `sg=n/a` without spglib; NumPy-2 `trapezoid` shim.)

`~/hyperspy-bundle/bin/python test_eels.py` → **10/10 passed**. Covers structure crystallography,
**CASTEP `.cell` writer validity (ASE round-trips the files, incl. the `O:exc` custom species)**,
the exactly-one-excited-atom invariant, the benchmark crystals, the M6 geometry model, the
OptaDOS parser (synthetic), and a submit-readiness check that all HPC templates exist. Re-run
this any time the code changes.

**M2(b) benchmark inputs READY:** `srtio3` (Pm-3m, a=3.905) and `tio2` (rutile P4₂/mnm,
a=4.593 c=2.959) build + validate; core-hole cells `srtio3_{Ti,O}` / `tio2_{Ti,O}` emitted.
So the first CASTEP calculation has its structures waiting.

## M0 — CASTEP + OptaDOS live on Blythe — ✅

CASTEP `24.1-foss-2023b` (`castep.mpi`); OptaDOS built from source (`$SHARE/phucrh/optados/
optados/optados.x`, big-endian). Smoke run: PbTiO₃ `tet_Pz` SCF `Final energy −8306.82 eV`,
and the reported Ti–O bonds (apical 1.745 Å / equatorial 1.980 Å) reproduce the ferroelectric
off-centering. Default OTFG pseudopotentials auto-generate (no SPECIES_POT for plain SCF).

## M2(a) — ferroelectric double-well — ✅ PASSED

Single-point SCF across `scan_0.00 → scan_1.00` (P4/mmm → P4mm, fixed exp. strain):

| s | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 |
|---|---|---|---|---|---|
| E − E(s=0), meV/f.u. | 0 | −30.5 | −106.3 | −185.4 | **−203.5** |

Monotonic; **polar 203 meV/f.u. below centrosymmetric**, and the descent flattens toward s=1
(last step −18 meV) → the **experimental displacement sits at the PBE minimum**. CASTEP+PBE
captures PbTiO₃ ferroelectricity. (Deeper than the ~50–100 meV cubic-referenced well because
the strained-non-polar reference is itself high-energy.)

## M2(b) — benchmark ELNES (rutile TiO₂ O-K) — ✅ PASSED, recipe LOCKED

`run_coreloss.slurm` on `tio2_O` (48-atom rutile, O:exc `{1s1}`, `charge:+1`): CASTEP wrote
`.elnes_bin` (big-endian), **OptaDOS read it cleanly** (endianness ✓), `_core_edge.dat` produced.
The physical spectrum is the `# O 1 K1 O:exc` block (OptaDOS writes one block per O atom; only
the excited one is the ELNES — **select, never sum**; `analyze_elnes.load_optados_core` does this).

Excited-atom O-K edge vs textbook rutile:
- **t₂g / e_g near-edge peaks split ~2.7 eV** (exp. ~2.5–3 eV) ✓
- **e_g taller than t₂g** (correct rutile intensity ordering) ✓
- broad O 2p–Ti 4sp band ~10–16 eV above onset ✓

→ the CASTEP 24.1 + OptaDOS core-hole recipe is **validated and locked**: `{1s1}`/`{2p5}`/`{3d9}`
OTFG core hole, `charge:+1`, `core_type:absorption`, O:exc-block selection.

## M3 — cubic null test — ✅ PASSED

Cubic Ti L₂,₃ (isotropic Ti site), q∥z vs q⊥z: dichroism **0.02%** (max|Δ|/max S and ∫|Δ|/∫S) —
the two spectra are identical, as symmetry demands. **The q-machinery is unbiased**; any
dichroism elsewhere is physical, not numerical. (OptaDOS writes the excited atom's full edge
family — Ti K/L1/L2,3/M1/M2,3; we take the L2,3 block.)

## M4 — tetragonal along-beam dichroism (tet_Pz, apical O-K) — ✅ LARGE SIGNAL

`tet_Pz` (c∥beam), apical O:exc K-edge, q∥c vs q⊥c: **max|Δ|/max S = 78%, ∫|Δ|/∫S = 47%** — a
clean σ/π anisotropy. q⊥c → sharp **π\*** peak ~8.5 eV (O p_{x,y}); q∥c → **σ\*** features
~11 & ~18.5 eV (O p_z along the Ti–O–Ti / polar axis). **Along-beam polar-axis orientation
imprints strongly on the O-K edge** — far above any detectability floor.

*Interpretation / caveat:* this apical-O dichroism is dominated by the **tetragonal backbone
σ/π anisotropy**, which the ferroelectric off-centering modulates. **M5** (dichroism/spectrum vs
displacement `s`) isolates the displacement-specific part (s=0 backbone baseline → s=1 polar);
also needed for the full O-K: multiplicity-weight 1×apical + 2×equatorial (`tet_Pz_Oeq`), and the
`tet_Pz` vs `tet_Px` rotational-invariance cross-check (pair in the crystal frame:
`tet_Px` has c∥x, so its `.qperp` output is the q∥c spectrum).

*Runtime note:* PbTiO₃ core-hole SCF ~40–45 min (32 cores); OptaDOS core-loss ~75–90 min/pass
(SERIAL) because it recomputes every atom/edge — the M5 bottleneck (see below).

## M5 — displacement scan (apical O-K, q∥c & q⊥c, s = 0 → 1) — ✅ DONE — the decomposition

Lattice held fixed (a, c); only the internal Ti/O displacement s varies, so the tetragonal
**strain** contribution is constant and the s-dependence is **purely the off-centering**.

| s | dichroism q∥c−q⊥ (peak) | along-beam (q∥c) spectrum change vs s=0 |
|---|---|---|
| 0.00 (zero-P) | **70.1%** | 0 |
| 0.50 | 70.5% | 4.3% |
| 1.00 (full P) | **78.2%** | **18.7%** |

**Key result:** ~70 of the 78% orientation dichroism is the **tetragonal backbone σ/π**, present
at *zero* polarization. The **ferroelectric off-centering's own signature** is a distinct,
**monotonic ~19%** change in the along-beam-probed O-K spectrum (and +8% to the dichroism),
cleanly isolated by the zero-P reference. Reproduces Bugnet et al. "O-K responds to Ti
off-centering" from first principles.

**Honest headline:** the O-K edge is 78% dichroic to the polar-**axis** orientation vs the beam
(what maps a vortex, since P and the tetragonal axis co-rotate); of that, ~19% is the
polarization's own displacement effect. Dipole limit: axis + magnitude, not sign.

PI-meeting figures (on ~/Desktop/eels_figs): fig1 validation, fig2 result + object + zero-P
overlay, fig3 decomposition/scan; SLIDE_NOTES.md.

## M4 cross-checks — ✅ DONE 2026-08-24 (B1 equatorial O, B2 rotational invariance)

Both jobs turned out to have completed in the M5 batch and never been analysed; the analysis
path itself had four bugs (see below), so none of this had ever been read.

### B2 — rotational invariance: ✅ PASS, essentially exact

```bash
~/hyperspy-bundle/bin/python analyze_elnes.py --compare tet_Pz_Oap tet_Px_Oap
```

| paired in the crystal frame | residual max\|Δ\|/max S | ∫\|Δ\|/∫S |
|---|---|---|
| `tet_Pz` q∥c vs `tet_Px` q∥c | **0.02 %** | 0.01 % |
| `tet_Pz` q⊥c vs `tet_Px` q⊥c | **0.03 %** | 0.01 % |

At the M3 cubic-null floor (0.02 %). **The 78 % apical dichroism is physical, not numerical.**

This is also the direct demonstration of the lab-vs-crystal q trap: pairing the files by NAME
instead of by crystal axis gives a residual of **78.2 %** — the entire signal — on a calculation
that is in fact exact to 0.0 %. Both HANDOFF.md and HANDOVER.md described the check that way.

### B1 — the multiplicity-weighted FULL O K edge

```bash
~/hyperspy-bundle/bin/python analyze_elnes.py --weighted tet_Pz
```

| site | max\|Δ\|/max S | ∫\|Δ\|/∫S |
|---|---|---|
| apical (1×) — *the previously published number* | **78.2 %** | 62.8 % |
| equatorial (2×) | **22.1 %** | 18.7 % |
| **weighted full edge (1 ap + 2 eq)** | **37.1 %** | 18.4 % |

**The headline nearly halves.** The apical oxygen sits on the Ti–O–Ti chain parallel to the
polar axis and is strongly dichroic; the equatorial oxygen, whose chain lies in-plane, is not.
Since the equatorial site carries twice the multiplicity, the edge of the cell as a whole is far
less dichroic than the apical site alone. **"The O K edge is 78 % dichroic" was never right; it
is the apical site that is 78 % dichroic, and the full edge is 37 %.**

The second equatorial orientation is now **measured, not estimated** (third OptaDOS pass,
`core_qdir 0 1 0`, `SKIP_CASTEP=1`). The 2c site holds two atoms, O_a at (½,0,z) and O_b at
(0,½,z), related by the C4 about z; `spectrum(O_b, x̂) = spectrum(O_a, ŷ)`. For **q∥c the two are
equivalent**, so counting one twice is exact there; for **q⊥c they are not**, and the measured
q=ŷ spectrum supplies the missing half.

**The earlier estimate was wrong, and not merely imprecise.** Counting O_a twice gave 46.4 % and
proxying O_b by the q=ẑ spectrum gave 39.1 %; those were quoted as a bracket on the true value.
The measured answer is **37.1 %, below both**. The proxy failed because q=ŷ is not a mild variant
of q=x̂ or q=ẑ: it differs from them by 64 % and 58 % respectively, and peaks at +19.0 eV rather
than +11.8 / +11.4 eV. Treat this as a caution about bracketing an unmeasured direction by
"similar" ones.

**Geometry self-check.** With Ti at (½,½,z), the 2c oxygen at (½,0,z) bonds along **y**, not x.
The spectra confirm it independently: q=ŷ (along its chain) gives a σ\*-like peak at **+19.0 eV**,
while q=x̂ and q=ẑ (both across it) give π\*-like peaks at **+11.8** and **+11.4 eV** — mirroring
the apical site, whose chain is along z and whose q=ẑ peak is at +18.7 eV.

**Remaining systematic: the chemical shift.** Still not fixed by these calculations (it needs the
Mizoguchi correction from a no-hole singlepoint). Sweeping ±1 eV on the weighted edge:

| shift (eV) | −1.0 | −0.5 | 0.0 | +0.5 | +1.0 |
|---|---|---|---|---|---|
| peak | 22.4 % | 31.5 % | **37.1 %** | 36.2 % | 34.3 % |
| integral | 15.6 % | 16.8 % | **18.4 %** | 20.3 % | 22.5 % |

So the weighted edge is **37 % assuming no shift, and 22–37 % across the sweep**. Quote it with
that systematic attached, or fix the shift.

### B3 — the `real` cell: the dichroism reports the polar-axis ORIENTATION, not its along-beam size

`real_Oap`, the apical oxygen of the labyrinth-derived cell whose polar axis is tilted 82° from
the beam (only ~14 % of |P| along it):

| cell | polar axis | max\|Δ\|/max S | ∫\|Δ\|/∫S |
|---|---|---|---|
| `tet_Pz` | fully along the beam | 78.2 % | 62.8 % |
| `real` | 82° from the beam (in-plane) | **76.3 %** | 60.4 % |

**Near-identical magnitude — which is not what "upper bound" framing predicts.** The reason is
that the character is *inverted*, and exactly so:

| cell | q along beam | q in-plane |
|---|---|---|
| `tet_Pz` (chain ∥ beam) | **σ\*** at +18.67 eV | π\* at +8.64 eV |
| `real` (chain in-plane) | **π\*** at +8.64 eV | σ\* at +18.67 eV |

The same two features at the same two energies, with their assignments swapped. In `real` the
whole tetragonal axis is rotated, not just the displacement: its short 1.745 Å Ti–O bond lies
along x and along the beam the titanium is centred (1.980 / 1.980 Å).

**Consequence.** The O K dichroism measures the **orientation of the polar axis relative to the
beam**, and *which* of q∥/q⊥ carries the σ\* says which way it points. It is not a magnitude that
shrinks toward zero as P rotates out of the beam. A measurement on the real labyrinth would
therefore still see a large dichroism, reporting "axis in-plane" rather than "no signal".

This does **not** rescue the along-beam component. M5 already showed only ~8 of the 78 points is
the off-centring's own signature, the rest being the tetragonal backbone; and the dipole limit
still forbids the sign of P. What changes is the framing: `tet_Pz` is not an upper bound that the
real geometry falls far below, it is one of two orientations that the edge distinguishes.

### Four bugs in the analysis path, all fixed

None of `analyze_seed` had ever run on real data (M4/M5 were read from hand-extracted
`.exc.txt`), and it contained: (1) a glob that could never match the slurm's output naming;
(2) `load_spectrum` instead of `load_optados_core`, i.e. summing all 162 per-atom blocks rather
than selecting `:exc`; (3) the lab-vs-crystal q-frame error above; and (4) **an analysis window
built from the ABSOLUTE edge energy (532 eV) when OptaDOS writes its axis RELATIVE to the edge**
(a real file spans −136.8 to +35.7 eV, onset +1.3 eV, near-edge structure +8 to +19 eV), so the
window selected zero points and the metric crashed.

**Metric conventions, now pinned by reproducing the published value.** `max|Δ|/max S` with
S = max over BOTH orientations reproduces 78.2 % exactly; normalising by the mean of the two
gives 105 %, which is not a usable fraction. The published **∫|Δ|/∫S = 47 % does NOT reproduce**
under any window or normaliser tried (closest 55 %, and 60–63 % on the natural conventions).
Recommend quoting the recomputed integral with its stated convention, or dropping it — the peak
metric is the one the paper leads with and it is exact.

Validation that the right block is being read: q⊥c gives a sharp π\* peak at **+8.64 eV** and
q∥c gives σ\* features at **+11 and +18.67 eV**, matching the M4 description (~8.5, ~11, ~18.5).

---

## M6 — detectability geometry — pending only the folding-through

tet_Px (rotational-invariance x-check) and tet_Pz_Oeq (full O-K) are **done** — see the M4
cross-checks section above. What remains is M6: folding the intrinsic dichroism through the
300 keV convergence/collection geometry (magic angle ≈ 4·θ_E ≈ 4.3 mrad; the 100 mrad
ptychography probe averages it away) for the measurable detectability number. The geometry model
itself is already validated (`analyze_elnes.py --selftest`); it now needs to be applied to the
**weighted** edge rather than the apical one, since the weighted edge is the physical observable.
