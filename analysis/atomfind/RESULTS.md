# atomfind — measured results & status

Everything here is measured on this machine, not estimated. Numbers are scored against
the ground-truth structure via the recon↔model map; the finder itself never sees GT.
**Status: portability pass complete** (2026-07-21) — the five portability bugs of the
previous §7 are fixed and verified on a second, independent reconstruction geometry.
Remaining open items in §7.

Reproduce: `~/hyperspy-bundle/bin/python atomfind/run_atomfind.py` (cwd
`ptychoshelves-clean/analysis`). Full machine-readable dump in `~/Desktop/atomfind_out/report.json`.

---

## 1. Headline — NL70 (0.15 Å, coherent, noiseless), gold Pb+Ti kernels

| metric | value | (pre-portability) |
|---|---|---|
| atoms found | 1878 (of which 53 lattice-guided, tagged) | 1761 (402 guided) |
| precision | 95 % | *97 %* |
| recall, all depths (Pb / Ti / O) | 88 % / 89 % / **90 %** | *88 / 91 / 83* |
| recall, **blind only** (Pb / Ti / O) | 88 % / 89 % / **86 %** | *88 / 89 / 50* |
| recall, bulk z 10–56 Å (Pb / Ti / O) | **96 % / 96 % / 96 %** | *96 / 97 / 92* |
| species confusion (off-diagonal) | **1.1 %** | *1.2 %* |
| in-plane accuracy (xy-RMS) | **0.032 Å** | *0.035 Å* |
| depth accuracy (z-RMS) | **0.37 Å** | *0.37 Å* |
| error-bar 1σ coverage (x / y / z) | 83 % / 92 % / 74 % (target ~68 %) | *87 / 92 / 71* |

The noise-relative CLEAN floor (§7) improved NL70 as a side effect: **bulk O 92 → 96 %**, and
**blind O recall 50 → 86 %** — oxygen is now mostly found directly by the model fit rather
than recovered by the lattice-guided pass, so the prior carries much less of the result.
The cost is precision 97 → 95 %: the extra ~35 detections are largely spurious *on this
noiseless volume* (bulk recall is unchanged-to-better), and they remain filterable via the
exported per-atom `quality`. That trade is what buys portability (§4); `clean_floor_k = 3.5`
restores 97 % precision here but costs dose1e10 Ti 74 → 63 %.

---

## 2. Detector comparison — including both deconvolution routes

All rows: same volume, same measured PSF, same validation. "pk" = 3-D local-maxima
peak-picking with sub-voxel refinement at its **best-effort threshold** (floor swept).

| method | Pb | Ti | **O** | precision | species labels | error bars |
|---|---|---|---|---|---|---|
| 3-D peak-pick, **raw** volume | 93 % | 88 % | **71 %** | 84 % | ✗ | ✗ |
| **RL deconv** + 3-D peak-pick | 93 % | 87 % | 58 % | 99 % | ✗ | ✗ |
| **MEM deconv** + 3-D peak-pick | 93 % | 88 % | 59 % | 86 % | ✗ | ✗ |
| v1: 1-D per-column NNLS spike deconv | 91 % | 65 % | 49 % | 84 % | ✗ | ✗ |
| **v3: 3-D model fit** (this work) | 88 % | 89 % | **90 %** | 95 % | ✓ | ✓ |
| v3, bulk only (z 10–56 Å) | **96 %** | **96 %** | **96 %** | | | |

### Key results
1. **Deconvolution-then-detect does not help, and RL actively hurts oxygen.** Both
   restoration routines *lose* O relative to peak-picking the raw volume (58–59 % vs 71 %):
   they concentrate dynamic range into the bright species so weak O falls below any
   relative floor, and they cannot restore the missing-cone frequencies that separate an
   O atom from the axial flank of its Ti neighbour. MEM ≈ RL here (Ishizuka et al. report
   MEM as the stronger routine; we reproduce no meaningful advantage on a *dense* lattice).
2. **On noiseless data a good peak-picker is competitive on the heavy atoms** — the honest
   comparison, and the reason our earlier "raw peaks = 24 % Ti" claim (a weak 1-D baseline)
   was retired.
3. **v3's edge is oxygen, labels, and uncertainty**: +19 pts O over the best baseline
   (90 %, and 96 % in bulk), 98.9 % correct species, calibrated per-atom σ. Its 95 %
   precision sits between the baselines' 84–99 %; the peak-pickers buy precision by
   simply not finding the hard atoms.

### Deconvolution implementations (both in `deconv.py`)
| routine | form | settings | role |
|---|---|---|---|
| Richardson–Lucy | multiplicative, non-negative | 15 iters (noiseless; dose-scaled otherwise), interior-trim, max-growth divergence guard | visualisation + baseline |
| Maximum entropy | Gull–Daniell exponential iteration | α = 0.5, 200 iters, stops on χ² stagnation | baseline only |

---

### Is the overlapped-oxygen margin a matching artefact? (checked 2026-08-24)

`cfg.match_tol_z_A` is **2.0 Å** and the axial overlap being scored is **1.95 Å**, so in
principle a detection sitting on the *titanium* can be credited as its apical oxygen, which
would inflate exactly the number the headline claim rests on. Re-scored the same detections at
tighter depth tolerances (no re-running of any finder — only the matcher changes):

| depth tolerance | v3 model fit | peak-pick raw | RL + peak-pick |
|---|---|---|---|
| 2.0 Å (shipped) | 0.821 | 0.258 | 0.092 |
| 1.0 Å | 0.796 | 0.258 | 0.098 |
| 0.6 Å | 0.725 | 0.246 | 0.087 |

**The margin is not an artefact.** The baselines are essentially *insensitive* to the tolerance
— they are not placing detections near the oxygen at all, at any tolerance — while v3 degrades
gracefully (0.82 → 0.80 → 0.73), which is the signature of detections that really do sit near
the oxygen sites and are progressively excluded as the gate closes below the localisation error.
The published 82% against 9–26% therefore stands, and would stand at a stricter gate.

**Caveat for the dose series (§4):** this reassurance is for the *coherent NL70* volume. On the
dose volumes the raw peak-picker scores a much higher overlapped-O recall than v3, while scoring
*lower* on in-plane-isolated oxygen — an ordering that is backwards from the physics (isolated
oxygen is the easier case) and therefore a warning that recall alone is not interpretable there.
Read it beside precision and the detection count; see §4.

---

## 3. Noise robustness — injected Gaussian noise, everything else fixed

Same geometry, kernel and calibration; only noise varies. Regenerated **2026-08-24** under the
current config with the seeded harness `noise_sweep.py` (there had been no committed harness —
the earlier table was produced ad hoc and predated the guided-deduplication guard, which is why
the report carried a caveat on it). Reproduce with:

```bash
python atomfind/noise_sweep.py --json noise.json      # seed 0, the ladder below
```

Measured intrinsic sub-median phase spread **σ ≈ 0.022** (the handover's estimate was 0.023);
an O peak is ≈ 0.06. The injected σ is on the same scale.

| injected σ | found | precision | Pb | Ti | O | **O overlapped** | **confusion** |
|---|---|---|---|---|---|---|---|
| 0 | 1834 | 0.97 | 88 % | 89 % | **89 %** | **82 %** | **1.1 %** |
| 0.01 | 1793 | 0.97 | 88 % | 91 % | **86 %** | 74 % | **1.5 %** |
| 0.02 (≈ noise floor) | 1421 | 0.98 | 87 % | 90 % | 57 % | 37 % | **3.7 %** |
| 0.04 |  795 | 0.96 | 87 % | 86 % |  5 % | 18 % | **58.7 %** |
| 0.08 |  408 | 0.95 | 87 % |  0 % |  0 % |  0 % | **74.6 %** |

Recalls are full-depth (not the bulk band). Versus the pre-dedup table the story is unchanged
and the numbers move only slightly (σ=0: 1878→1834 found, precision 0.95→**0.97**, the
deduplication guard removing duplicate detections rather than finding fewer atoms).

**The collapse at σ ≳ 0.02 is not a thresholding artifact**, it is oxygen falling below the
noise, and no floor rule fixes it. The new **O overlapped** column is the sharper statement:
axially overlapped oxygen (1.95 Å from a Ti along the beam) degrades 82 → 74 → 37 → 18 → 0 %,
i.e. it is the population that the noise takes first, exactly as it is the population the method
wins on when the data are clean.

**Precision is not a safety indicator.** At σ = 0.08 recall for Ti and O is literally zero and
three in four species labels are wrong, yet precision still reads 0.95 (the surviving bright Pb
are correctly placed). **Confusion rate and σ-coverage are the canaries**; `noise_sweep.py`
prints `validate.health_warnings` per rung, and they fire from σ = 0.04 upward. The failure mode
is at least conservative: the lattice-consistency amplitude gate makes guided-O *abstain* rather
than fabricate.

---

## 4. Dose series — real reconstruction noise (NL105, 0.1 Å, coherent)

Re-run **2026-08-24** under the current config (removes the pre-deduplication caveat), now with
the deconvolve-then-peak-pick baselines scored on the same volumes:

```bash
python atomfind/dose_series.py --json dose.json          # all four doses + baselines
```

### v3, bulk band (z = 10–58 Å)

| dose (e Å⁻²) | found | precision | Pb | Ti | O | confusion |
|---|---|---|---|---|---|---|
| 10¹⁰ | 1575 | 0.85 | 89 % | 74 % | 60 % | 12.4 % |
| 10⁸  | 1593 | 0.88 | 89 % | 77 % | 66 % | 11.3 % |
| 10⁶  |  521 | 0.96 | 76 % | 26 % |  5 % | 37.8 % |
| 10⁴  |    0 | —    | —    | —    | —    | —      |

Unchanged in substance from the pre-dedup run: usable to 10⁸, collapsing below 10⁶. At 10⁴ the
finder returns **nothing at all** rather than reporting noise, which is the intended behaviour.

### The head-to-head — and why it does NOT discriminate at dose

This was run to test the expectation that the oxygen margin *widens* with dose. **It does not,
and the reason is a measurement artefact in the baseline, not a loss for either method.**

| dose | detector | emitted | precision | O bulk | O overlapped | O isolated |
|---|---|---|---|---|---|---|
| 10¹⁰ | v3 | 1575 | **0.85** | 60 % | 46 % | 65 % |
|      | peak-pick raw | **2500 (cap)** | 0.56 | 59 % | 71 % | 49 % |
|      | RL + peak-pick | 2199 | 0.55 | 54 % | 35 % | 59 % |
|      | MEM + peak-pick | 1955 | 0.64 | 49 % | 35 % | 51 % |
| 10⁸  | v3 | 1593 | **0.88** | 66 % | 46 % | 72 % |
|      | peak-pick raw | **2500 (cap)** | 0.56 | 59 % | 77 % | 45 % |
|      | RL + peak-pick | 1781 | 0.67 | 56 % | 36 % | 58 % |
|      | MEM + peak-pick | 1802 | 0.68 | 48 % | 35 % | 49 % |
| 10⁶  | v3 | 521 | 0.96 | 5 % | 22 % | 0 % |
|      | MEM + peak-pick | 3165 | 0.45 | 69 % | 48 % | 72 % |
| 10⁴  | v3 | **0** | — | — | — | — |
|      | MEM + peak-pick | 3803 | 0.34 | 65 % | 49 % | 67 % |

Three facts make the baseline's apparent advantage uninterpretable:

1. **It is saturated.** `peaks3d` takes `max_atoms` (2500 for raw/RL, 4000 for MEM). Raw
   peak-picking emits **exactly 2500 at every one of the four doses** — the cap, not the data,
   is setting its output. At 10¹⁰ and 10⁸ it returns *identical* true-detection counts (1407)
   and *identical* false positives (1093), despite a hundredfold change in dose.
2. **Its recall ordering is physically impossible.** It scores axially overlapped oxygen (71 %,
   77 %) *higher* than in-plane-isolated oxygen (49 %, 45 %). Isolated oxygen occupies its own
   column and is the easy case for any detector. A method that finds the hard population more
   often than the easy one is not detecting oxygen; it is carpeting the volume and being
   credited for whatever it lands near.
3. **The true-detection counts are the same.** At 10⁸, v3 makes 1404 true detections from 1593
   emitted; raw peak-picking makes 1407 from 2500, buying the same true detections with
   **1093 false positives** (44 % of its output) against v3's ~190 (12 %).

Below 10⁶ this becomes stark: v3 abstains entirely at 10⁴ while MEM reports 49 % overlapped-O
recall at precision 0.34, i.e. it reports oxygen where there is no signal.

**Consequence for the paper.** The discriminating comparison is the *noiseless* one (§2), where
no detector is saturated — v3 1834, raw 1884, RL 1436, MEM 1670 emitted, at precision 0.84–0.99
— so recalls are commensurable. At finite dose the comparison stops being a comparison, because
the baseline's recall is set by its detection cap. Any claim that the margin widens with dose
must **not** be made; what dose actually demonstrates is the value of a method that abstains
(v3: 0 atoms at 10⁴) over one that does not (MEM: 3803 detections at precision 0.34).

A fair finite-dose comparison would need the baselines re-thresholded to emit a matched number
of detections, or scored on a precision-recall curve rather than at one operating point. That is
worth doing and is left as further work; it does not affect §2, which is the claim the paper
actually makes.

---

## 5. Kernel inventory (measured, not from the handover notes)

| kernel | shape | argmax | SNR | axial FWHM | verdict |
|---|---|---|---|---|---|
| `psf_Pb_NL70` | (46,60,60) | (24,30,30) centred | 81.5 | ~1.0 Å | **gold — NL70 default** |
| `psf_Ti_NL70` | (46,60,60) | (24,30,30) centred | 24.2 | — | clean, used for species shape |
| `psf_Pb_rev2_d1e10` | (81,60,60) | (43,30,30) centred | 23.6 | 2.0 Å | clean (NL105, 16-phonon) |
| `psf_Pb_rev2_d1e8` | (81,60,60) | (43,30,30) centred | 18.8 | — | clean (NL105, 16-phonon) |
| `psf_Ti_rev2_*` | (81,60,60) | corner (62,59,59) | 2.4 | — | **broken**; interior peak off-centre, SNR 4.7 |
| `psf_O_rev2_*` | (81,60,60) | corner (63,59,59) | 2.7 | — | **broken**; d1e10 interior peak on the *last layer* |
| data-derived (from dose1e10 Pb blobs) | — | — | — | 4.0 Å | **worse** (Pb 65 % vs 85 %) — in-crystal neighbours inflate the axial width |

Confirms the physics in `docs/history/PSF_SIM_RESPONSE.md`: isolated light atoms do not reconstruct
above the noise floor, and worse under TDS+dose. Use the **Pb shape for all species**.

---

## 6. Uncertainty quantification (rebuilt 2026-07-21)

Two explicitly separated stages (`find.py` model σ + `uncertainty.py` calibration). The old
tuned per-species floors are **retired**; they could only hit coverage on average and left
the overlapped oxygen optimistic and isolated atoms conservative (the failure below).

**Stage 1 — model σ (no ground truth):** σ² = σ²_stat + σ²_sys.
- σ_stat is the **JOINT Cramér–Rao bound**: the diagonal block of the inverse of the *full*
  Fisher matrix over all atoms in a tube, not the per-atom block `[JᵀJ]⁻¹`. The block form is
  the bound with every neighbour *known exactly*; the ratio to the joint bound is the variance
  inflation factor (VIF) of the overlapping design. **Measured on our kernel:**

  | separation | pair | VIF(β) | VIF(z) | σ_z understated |
  |---|---|---|---|---|
  | 3.90 Å | Pb–Pb, Ti–Ti | 1.04 | 1.04 | ×1.02 |
  | 1.95 Å | Ti–O (every apical O) | 2.28 | 2.02 | ×1.42 |
  | 0.99 Å | near-degenerate | 13.3 | 3.33 | ×1.83 |

  and **85 % of in-window atoms have a neighbour within 2.5 Å**, so the block form was
  optimistic almost everywhere. The joint inverse restores the inflation per atom, no tuned
  constant. It also gives σ_β, which drives the species posterior (below).
- σ_sys is the **kernel-mismatch** term (`find.kernel_mismatch_sigma`): refit a synthetic
  response made with one measured kernel using another, take the position spread. **Computable
  without GT, so it transfers.** On NL70 (Pb vs Ti kernels nearly identical) it is small
  (σ_z 0.02 Å); it grows on datasets whose kernels genuinely differ.

**Stage 2 — split-conformal calibration (Mondrian):** nonconformity s = |Δ|/σ per axis,
empirical (1−α) quantile *per stratum* (species × blind/guided × depth band), half-width =
q·σ. **Coverage holds per stratum by construction** — no Gaussian assumption, no hand-tuned
constant. Verified on NL70:

| target | overall coverage (x / y / z) | worst-stratum z |
|---|---|---|
| 68 % | 69 / 69 / 69 % | 68 % (O\|blind\|bulk) |
| 95 % | 96 / 96 / 96 % | 91 % (O\|guided\|entrance, n=11 — small-n conformal noise) |

The per-stratum q(z) ranges **7.3 → 14.3** — direct proof a single floor cannot represent the
uncertainty. Exports: `uq_conformal.json` (the q-table) and, per atom, the model σ, the **95 %
half-width as default** (a 1σ number is the least-conservative honest choice; the conservative
interval is what a downstream user picks up by accident), and a within-column amplitude
posterior `p_species` (relative confidence, not calibrated; atoms with p<0.9 are **~25× enriched
for mislabelling** over the 1 % base rate).

**Precision context:** model σ ≈ FWHM/10 per axis (in-plane σ_stat 0.001–0.002 Å vs kernel
FWHM 0.1 Å; axial 0.03 Å vs 1–3 Å) — the width/SNR scaling for a known template. On noiseless
simulation σ_stat → 0 and the *systematic* error dominates; conformal supplies that scale
(hence the large q), and σ_stat resumes its role on counting-limited data.

**Transfer — the decisive test.** Applying an NL70 q-table to a new volume assumes
exchangeability, which fails when registration quality differs. The intended workflow is to
re-run the conformal step on each volume's own matched atoms. On **dose1e10** (different
geometry NL105, ~half the phase amplitude, real recon noise, z-RMS 0.87 Å) this gives, with
the *same code and no tuned constants*:

| target | dose1e10 coverage (x / y / z) | vs old tuned floors |
|---|---|---|
| 68 % | **69 / 69 / 69 %** | 38 / 47 / 56 % (under-covered, warned) |
| 95 % | **96 / 96 / 96 %** | — |

q(z) spread on dose is even wider (10.9 → 23.3) — the errors are more heteroscedastic on the
harder volume, which is exactly the regime a single floor cannot serve and conformal handles
for free. Only σ_sys (kernel mismatch) transfers unconditionally.

---

## 7. Open items — read this first when resuming

**Portability bugs — ALL FIXED (2026-07-21), verified on the dose series (§4):**
1. ✅ **`clean_floor` made noise-relative.** Now `k · σ_MF`, where σ_MF is a per-tube MAD of
   the matched-filter response over quiet (sub-median) voxels — blind, no GT.
   `find.clean_floor_for`. `k = 2.0` chosen by sweeping **both** volumes (recall is monotone
   in k on each): it is the only setting meeting the dose targets with one config. Side
   effect on NL70: bulk O 92 → 96 %, blind O 50 → 86 %, precision 97 → 95 % (§1).
2. ✅ **Depth-offset refinement added.** `align.refine_with_atoms` now alternates depth and
   in-plane fits for 4 iterations — it *must* iterate, because the correspondences that
   diagnose the offset are themselves computed with the wrong offset. The z-gate (1.8 Å) is
   held below the 3.9 Å heavy-atom column spacing so a match cannot jump a lattice site. Also
   fits an optional depth-scale residual `mZ`. Two follow-on fixes proved necessary:
   **fiducials are selected by amplitude, not species label** (on a misregistered volume the
   labels are exactly what is wrong, so label-based selection feeds the error back into
   itself: dose1e10 OFF 0.43 → 0.29 Å, Ti 71 → 74 %), and OFF absorbs the shift so it stays
   interpretable. Measured: dose1e10 OFF **1.26 → 0.29 Å**, Ti **53 → 74 %**, confusion
   **19.0 → 13.7 %**.
3. ✅ **`clean_nms_z_layers` → `clean_nms_z_A`**, converted to layers per-volume at use.
4. ⚠️ **σ floors now PSF-scaled** (blind, kernel-FWHM relative) — dose1e10 z-coverage
   28 → 56 %. *Partially* fixed: the residual is misregistration, not blur, and is not
   knowable blind. See §6 — the coverage metric now warns when they don't hold.
5. ✅ **Health metrics printed at run time.** `validate.health_warnings` flags high species
   confusion, out-of-range σ-coverage, and collapsed per-species recall; the verdict prints
   them. dose1e6 is the case in point: precision 0.95, confusion 38 %.

**Still open:**
- **σ inflation factors for guided and exit-band atoms are now stale** (found 2026-07-21
  during doc verification, §6): the guided population collapsed 402 → 53 atoms as a *result*
  of the floor fix, and the surviving hard residue is under-covered (z 60 % at 1σ, 83 % at
  2σ); the exit band is likewise at 78 % at 2σ. Re-derive `guided_sigma_scale` and
  `exit_sigma_scale` on the current population. Small, contained job.
- **Dose-series absolute performance** (Ti 74 %, z-RMS 0.87 Å at 1e10) is well short of NL70
  and is dominated by reconstruction quality (corr_depth 0.44, 2.0 Å axial kernel), not by
  our constants. A **coherent NL105 kernel** (all current NL105 kernels are 16-phonon) is the
  most likely single improvement.
- **σ floors do not fully transfer** (item 4) — re-derive per dataset when coverage warns.
- **1e6 / 1e4 are not usable** and were not chased (deliberate).

**Pending data (owner: sims thread):** in-situ vacancy-difference kernels — the labyrinth
reconstructed minus exactly one atom (priority: O on a B–O column, then O on a pure-O
column, then Ti, then Pb). Spec in `docs/history/PSF_SIM_REQUEST.md` §"REQUEST 2". Would give a true
in-crystal O matched filter, which §5 shows is the one thing no existing kernel provides.
Also outstanding: a **coherent NL105 kernel** (all NL105 kernels are 16-phonon), and
re-extraction of the broken `psf_Pb_NL70_z10/z64` depth kernels with a smaller `--zdrop`.

**Known physics limitation:** per-species Debye–Waller factors differ (O B ≈ 0.80 vs
Ti 0.45 Å²), so under TDS oxygen's true kernel is genuinely broader than lead's — the
"one K for all species" assumption weakens exactly where it matters most.

---

## 8. Polarisation from the found atoms (added 2026-08-15)

`polarisation.py` — Ti-O6 off-centring (Eq. delta = r_Ti - centroid(6 nearest O)) computed from
the BLIND found atoms, uncertainty propagated by Monte Carlo from the exported 95% conformal
half-widths (no GT). Scored against the model. Run:
`~/hyperspy-bundle/bin/python atomfind/polarisation.py` (cwd `analysis/`).

NL70, bulk band z 10-56 A: **72%** of bulk Ti get a complete O6 cage from found O alone
(247 GT Ti in window -> 172 caged -> 171 matched to GT).

| component | median \|err\| | p90 | RMS | GT spread (sd) | r | MC sigma | 95% cov |
|---|---|---|---|---|---|---|---|
| delta_x | **0.005 A** | 0.013 | 0.080 | 0.272 | 0.96 | 0.011 | 95% |
| delta_y | **0.004 A** | 0.014 | 0.060 | 0.078 | 0.80 | 0.009 | 93% |
| delta_z | 0.122 A | 0.304 | 0.199 | **0.074** | 0.42 | **0.237** | 99% |

- **In-plane is quantitative**: vector error median 0.007 A / p90 0.017 A against a 0.274 A
  signal; DIRECTION error median 0.9 deg, p90 3.7 deg, 95% within 30 deg. Excluding the tail
  below, corr(recon,GT) delta_x = 1.000.
- **Along-beam is NOT measured**: propagated sigma_z (0.237 A) exceeds the entire GT spread of
  delta_z (0.074 A). This is knowable BLIND (compare propagated sigma to the recovered spread)
  — the UQ correctly declares the along-beam map uninformative. Motivates the EELS phase.
  Root cause: delta is a DIFFERENCE of coordinates; z-RMS 0.37 A per atom cannot resolve a
  0.07 A displacement, while xy-RMS 0.032 A resolves a 0.27 A one easily.
- **The 5% tail is a DETECTION failure, not a localisation one**: 8/171 Ti have a cage of the
  wrong composition (missing/misassigned O) -> in-plane error median 0.45 A. The propagated
  interval covers 99% of the good Ti but only **12%** of the tail — as expected, since the
  conformal interval is conditional on correct detection (§6). Report cage completeness +
  precision + confusion alongside any polarisation map; do NOT absorb this into the error bar.

---

## 9. Files
`METHODS.md` (method, plain-language) · `README.md` (usage) · `paper/atomfind_methods.tex`
+ `paper/atomfind_refs.bib` (manuscript segment, 6 refs, compiles clean) ·
`docs/history/PSF_SIM_REQUEST.md` / `docs/history/PSF_SIM_RESPONSE.md` (kernel provenance) ·
`dose_series.py` (portability harness: NL105 `.mat` loader + per-dose run) · this file.

Reproduce the portability benchmark:
`~/hyperspy-bundle/bin/python atomfind/dose_series.py [1e10 …]`
