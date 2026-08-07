# EELS on Blythe — operational handover

Everything a future agent (or a returning human) needs to run the PbTiO₃ core-loss EELS study
on the Warwick **Blythe** HPC. Companion docs: **README.md** (scientific milestone map),
**RESULTS.md** (evidence log). This file is the *operator's* guide: environment, sync flow,
the CASTEP+OptaDOS recipe, how to run each step, and the hard-won gotchas.

---

## 0. TL;DR — current state (as of 2026-08-07)

The question: **does the ferroelectric polarisation component ALONG the beam produce a
detectable EELS/ELNES signal** in PbTiO₃ — i.e. can EELS recover the 3-D vortex polarisation
that projected ptychography (in-plane only) misses? Method: **CASTEP** core-hole ELNES +
**OptaDOS** core-loss, with OptaDOS's `core_qdir` as the "beam direction" knob.

| Milestone | What | State |
|---|---|---|
| M0  | CASTEP 24.1 + OptaDOS on Blythe | ✅ both live (CASTEP smoke ran; OptaDOS built) |
| M1  | build + validate PTO/benchmark cells | ✅ 11/11 unit tests; HPC preflight green |
| M2a | ferroelectric double-well | ✅ **PASSED** — polar 203 meV/f.u. below centrosymmetric, min at exp. displacement |
| M2b | TiO₂/SrTiO₃ benchmark ELNES → lock recipe | ⏳ **in progress** — first TiO₂ O-K core-loss running |
| M3  | cubic null test (dichroism must vanish) | ☐ |
| M4  | tetragonal q∥c vs q⊥ dichroism (the result) | ☐ |
| M5  | dichroism vs \|P\| calibration | ☐ |
| M6  | detectability model (analytic magic-angle) | ✅ validated locally + on Blythe |
| M6b | abtem multislice EELS (Channel B) | ◐ elastic path works on Blythe; core-loss needs gpaw |

**Immediate next step:** finish M2(b) — read the TiO₂ O-K spectrum out of `runs/tio2_OK/`,
confirm endianness + shape vs the textbook rutile O-K edge; capture the **Pb OTFG string** to
complete the PbTiO₃ core-hole cells; then M3 (null) → M4 (dichroism).

---

## 1. The two machines and the sync flow (important)

Code is **edited on the Mac**, pushed to GitHub, and **pulled on Blythe**. There is no editing
on Blythe. Never assume the HPC has an uncommitted change.

```
Mac (dev)  ──git push──▶  github.com/tristanmccarthy55/ptyco_pipeline (branch: main)  ──git pull──▶  Blythe
```

- **Mac repo:** `/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/ptychoshelves-clean` (branch `main`). The `eels/` subproject is peer to `sim/` and `analysis/`.
- **Blythe repo:** `$SHARE/phucrh/ptyco_baseline/ptyco_pipeline` (branch `main`; the HPC was historically on `baseline-worked-example` — we moved it to `main`). `eels/` is at `.../ptyco_pipeline/eels`.
- **Sync a change:** on Mac `git add eels/ && git commit && git push origin main`; on Blythe `git pull origin main`. (Stage only `eels/` — the repo has unrelated atomfind-paper WIP.)
- `$SHARE` = `/springbrook/share/physics`, so `$SHARE/phucrh` = `/springbrook/share/physics/phucrh`. **`$HOME` on Blythe is tiny — keep everything on `$SHARE`.**

---

## 2. Blythe environment (paths, modules, binaries)

| Thing | Value |
|---|---|
| Login | `ssh phucrh@blythe.scrtp.warwick.ac.uk` |
| Shared FS | `$SHARE = /springbrook/share/physics`; project space `$SHARE/phucrh` |
| CPU partition | **`compute`** (default `*`; 21× `trent`, 2-day limit). Often 0 idle → jobs pend. |
| Quick-start partitions | `hmem` (humber) + `int` (cole) frequently have idle nodes — use for tiny test jobs |
| Python env | **`$SHARE/phucrh/envs/abtem/bin/python`** (Python 3.11.15) |
| CASTEP | module **`phys/CASTEP/24.1-foss-2023b`** → `castep.mpi` (MPI, gfortran 13.2, OpenMPI 4.1.6) |
| OptaDOS | **`$SHARE/phucrh/optados/optados/optados.x`** (built here; serial; v1.2.380) |
| CUDA (for cupy) | `module load CUDA/12.6.0` |

### The Python env (abtem)
Call it by **absolute path** — the conda shell hook is unreliable in batch jobs and silently
falls back to a wrong python. Override with the `CONDA_ENV` env var (the SLURM scripts honour it).

```
$SHARE/phucrh/envs/abtem/bin/python <script>
```
Has: `abtem 1.0.9`, `cupy-cuda12x`, `dask`, `h5py`, `ase`, `scipy`, `numpy` (2.x).
**Does NOT have:** `spglib`, `pymatgen`, `gpaw`. The code degrades gracefully:
- no spglib/pymatgen → `build_cells` prints `sg=n/a`, skips the space-group *string* check
  (lattice + Ti-displacement gates still run). `pip install spglib` re-enables it (optional).
- NumPy 2.x removed `np.trapz` → `analyze_elnes` uses a `np.trapezoid` shim.
- `gpaw` is needed only for the M6b abtem *core-loss* transition potentials (the elastic path
  runs without it).

### CASTEP (24.1)
- Binary path: `/springbrook/apps/software/CASTEP/24.1-foss-2023b/bin/castep.mpi`.
- **Default OTFG pseudopotentials** are generated on the fly — a plain SCF needs **no
  `SPECIES_POT` block**. (For core holes we DO specify it, see §4.)
- The `foss-2023b` toolchain = GCC 13.2 + OpenMPI. `module load foss/2023b` gives `gfortran`
  (used to build OptaDOS).

### OptaDOS (built from source)
Not packaged as a module — built once:
```bash
cd $SHARE/phucrh
git clone https://github.com/optados-developers/optados.git
cd optados/optados          # NB: source + Makefile are in the NESTED optados/ dir
module load foss/2023b
make SYSTEM=gfortran        # -> ../optados.x  (i.e. $SHARE/phucrh/optados/optados/optados.x)
```
It's a **debug build** (`-O0 -fcheck=all`, fine for light post-processing) and **big-endian**
(`-fconvert=big-endian`). CASTEP must WRITE its binaries big-endian for OptaDOS to read them
(the classic gotcha) — verified/checked on the first core-loss read (M2b). Rebuild optimised
later if it ever feels slow.

---

## 3. The CASTEP + OptaDOS core-loss recipe

Decoded from the OptaDOS examples that ship with the source (`optados/examples/Mizoguchi`
= NaGe core hole; `optados/examples/Si2_CORE`). This is the authoritative CASTEP-24.1 syntax.

**Workflow (two steps, one seed):** `castep.mpi <seed>` (core-hole SCF) → `optados <seed>`
(reads `<seed>.odi`, `task:core`). `run_coreloss.slurm` does both.

**Core-hole `.cell`** — excited atom gets a distinct label and a core-hole OTFG string:
```
%BLOCK POSITIONS_FRAC
  ...
  O:exc  0.40 0.60 0.25          # ONE excited atom
  ...
%ENDBLOCK POSITIONS_FRAC
KPOINT_MP_GRID 4 4 4
SPECTRAL_KPOINTS_MP_GRID 4 4 4   # plural form (matches Si2_CORE)
%BLOCK SPECIES_POT
O  2|1.1|17|20|23|20:21(qc=8)
Ti 3|1.8|9|10|11|30U:40:31:32(qc=5.5)
O:exc 2|1.1|17|20|23|20:21(qc=8){1s1}    # normal O string + core-hole suffix
%ENDBLOCK SPECIES_POT
```
`build_cells.py --corehole` writes these automatically (see §4). **No `SYMMETRY_GENERATE`** in a
core-hole cell (the hole breaks symmetry).

**CASTEP-24.1 default OTFG strings** (from the M0 smoke-run pseudopotential report):
| element | OTFG string |
|---|---|
| O  | `2\|1.1\|17\|20\|23\|20:21(qc=8)` |
| Ti | `3\|1.8\|9\|10\|11\|30U:40:31:32(qc=5.5)` |
| Pb | **TODO** — `grep -A1 "Element: Pb" runs/smoke/smoke.castep` (needed for PbTiO₃ M4) |
| Sr | TODO (only if an SrTiO₃ core hole is run) |

**Core-hole occupancy suffix per edge** (reduce the ionised shell by one electron):
`O K → {1s1}` · `Ti L₂,₃ → {2p5}` · `Pb M₄,₅ → {3d9}` · `Sr L₂,₃ → {2p5}`.

**`.param`** (core hole): `task:spectral`, `spectral_task:coreloss`, `xc_functional:PBE`,
`cut_off_energy:800 eV`, **`charge:+1`** (= full core hole), `nextra_bands:100`. (`templates/coreloss.param`.)

**`.odi`** (OptaDOS): `task:core`, `core_type:absorption`, adaptive broadening + LAI widths.
For the **anisotropy** (M4) add `core_geom:polarized` + `core_qdir`:
- `coreloss.odi` — polycrystalline (isotropic) → the M2b benchmark spectrum.
- `coreloss_qc.odi` — `core_qdir 0 0 1` = q∥c = **P along beam**.
- `coreloss_qperp.odi` — `core_qdir 1 0 0` = q⊥c. **Dichroism = spectrum(qc) − spectrum(qperp).**
  (For `tet_Px` the axes swap → use `1 0 0` for ∥, `0 0 1` for ⊥.)

Absolute edge onset (optional, for comparing to experiment) = the Mizoguchi chemical shift from
a separate no-hole `singlepoint` + `optados/tools/miz_chemical_shift`. **Not needed** for the
dichroism (a difference).

---

## 4. The `eels/` file map

| File | Role |
|---|---|
| `config.py` | reference crystallography, structure presets, edges, **OTFG strings + core-hole occupancies**, q-dirs, optics |
| `build_cells.py` | build/validate cells → CASTEP `.cell`; `--corehole` writes complete core-hole cells (SPECIES_POT). Needs ase (+ spglib/pymatgen optional) |
| `analyze_elnes.py` | M4–M6: dichroism, calibration, magic-angle geometry model, SNR (`--selftest`) |
| `simulate_eels.py` | M6b: abtem multislice STEM-EELS (Channel B); `--selftest` (no gpaw), core-loss needs gpaw |
| `test_eels.py` | CASTEP-independent unit suite — **run after any edit** (`python test_eels.py`, 11/11) |
| `templates/groundstate.param`,`geomopt.param` | M2a SCF / relaxation |
| `templates/coreloss.param` | core-hole spectral SCF (`charge:+1`) |
| `templates/coreloss.odi`,`coreloss_qc.odi`,`coreloss_qperp.odi` | OptaDOS core task (benchmark / q∥ / q⊥) |
| `run_pyeels.slurm` | CASTEP-independent preflight (uses the abtem env by absolute path) |
| `run_castep.slurm` | generic single CASTEP job (M2a) |
| `run_coreloss.slurm` | CASTEP core-hole SCF **then** OptaDOS (M2b/M3/M4) |
| `structures/*.cell` | generated geometries (`_1cell` for M2a, `_222` for core holes) |
| `structures/corehole/*.cell` | one excited site per file; TiO₂ have full SPECIES_POT, PbTiO₃ geometry-only until Pb OTFG captured |
| `runs/` | **gitignored** — all CASTEP/OptaDOS job dirs + outputs live here on `$SHARE` |

**The cells** (all beam = +z, matching `sim/simulate_4dstem.py`'s `rotate(-90°,y)`):
`cubic` (Pm-3m, null test) · `tet_Pz` (P4mm, c∥z = **P along beam**) · `tet_Px` (c∥x = P
perpendicular) · `scan_0.00…1.00` (displacement ladder, M2a double-well + M5 calibration) ·
`real` (labyrinth-derived tilted P) · `srtio3`,`tio2` (M2b benchmarks).

---

## 5. How to run each step (copy-paste-safe; no `<placeholders>`)

All from `cd $SHARE/phucrh/ptyco_baseline/ptyco_pipeline/eels` after `git pull origin main`.

**Preflight (CASTEP-independent):**
```bash
$SHARE/phucrh/envs/abtem/bin/python test_eels.py        # 11/11
$SHARE/phucrh/envs/abtem/bin/python build_cells.py --corehole
# or via SLURM: sbatch -p compute run_pyeels.slurm
```

**M2(a) double-well (CASTEP only) — DONE, rerun template:**
```bash
for s in 0.00 0.25 0.50 0.75 1.00; do
  d="runs/dw_$s"; mkdir -p "$d"
  cp structures/scan_${s}_1cell.cell "$d/dw_$s.cell"; cp templates/groundstate.param "$d/dw_$s.param"
  ( cd "$d" && sbatch -p compute --parsable --export=ALL,SEED=dw_$s ../../run_castep.slurm )
done
grep -H "Final energy" runs/dw_*/dw_*.castep   # energy must fall from s=0 to s=1
```

**M2(b) benchmark ELNES (CASTEP + OptaDOS):**
```bash
d=runs/tio2_OK; mkdir -p "$d"
cp structures/corehole/tio2_O.cell "$d/tio2_O.cell"
cp templates/coreloss.param        "$d/tio2_O.param"
cp templates/coreloss.odi          "$d/tio2_O.odi"
( cd "$d" && sbatch -p compute --parsable --export=ALL,SEED=tio2_O ../../run_coreloss.slurm )
# outputs: tio2_O*.dat / .agr (spectrum), tio2_O.odo (OptaDOS log)
```

**M3 null test** (cubic — dichroism MUST vanish): same as M2b but a `cubic_*` core-hole cell,
run OptaDOS twice with `coreloss_qc.odi` and `coreloss_qperp.odi`; the two spectra must coincide.

**M4 dichroism** (the result): `tet_Pz_Oap`/`tet_Pz_Ti`/… core-hole cells (need Pb OTFG for Pb),
one CASTEP SCF each, then OptaDOS with qc and qperp:
```bash
( cd runs/tetPz_Oap && sbatch -p compute --parsable \
  --export=ALL,SEED=tet_Pz_Oap,ODI2=coreloss_qperp.odi ../../run_coreloss.slurm )
```
Cross-check: `tet_Pz` q⊥ spectrum must equal `tet_Px` q∥ spectrum (rotational invariance).

**Analysis (pull the `.dat` back to the Mac, or run on Blythe):**
```bash
$SHARE/phucrh/envs/abtem/bin/python analyze_elnes.py --seed tet_Pz_Oap --edge O_K
```
(expects `<seed>.qc.*core*.dat` / `<seed>.qperp.*core*.dat` in `eels/runs/`).

---

## 6. SLURM & environment gotchas (hard-won — read before debugging)

- **Never paste literal `<placeholder>`** into a command — bash reads `<` as a stdin redirect
  (`-bash: cpu-partition: No such file or directory`). Substitute real values (`-p compute`).
- **`-o logs/…` fails if `logs/` doesn't exist** — SLURM won't create the parent dir, the job
  dies at 00:00 FAILED with no `.out`. Our scripts log to the **submit dir root** (`*_%j.out`)
  to avoid this. `logs/` and `runs/` are gitignored.
- **Capture the job id** with `--parsable`: `JOBID=$(sbatch -p compute --parsable … script)`.
  Then `squeue -j "$JOBID"`, `sacct -j "$JOBID" --format=JobID,State,ExitCode,Elapsed`.
- **`compute` is often full** (0 idle) → `PD (Resources)`. For tiny test jobs use `-p hmem` or
  `-p int` (idle nodes) for an instant start.
- **Python env by absolute path**, never `conda activate`, in batch (§2).
- **CASTEP↔OptaDOS endianness** must match (both big-endian here). If OptaDOS reads garbage or
  errors on a `.bands`/`.elnes_bin`, that's the culprit — rebuild one to match the other.
- **PBE over-tetragonalises PbTiO₃** ("supertetragonal") — expect c/a slightly high if you relax;
  we build at experimental lattice by default.

---

## 7. Key results so far (see RESULTS.md for detail)

- **M2(a) ferroelectric double-well:** polar P4mm sits **203 meV/f.u.** below centrosymmetric
  P4/mmm at fixed experimental strain, monotonic in displacement, **minimum at the experimental
  distortion** → CASTEP+PBE captures the ferroelectricity correctly. (Deeper than the ~50–100
  meV cubic-referenced value because the strained-non-polar reference is itself high-energy.)
- **Real-geometry finding (from the labyrinth):** viewed side-on (the pipeline's zone axis) the
  vortex polarisation is **mostly in-plane** — median angle-from-beam 82°, only **~14% along
  beam**. So the along-beam signal is intrinsically small here; `tet_Pz` (full along-beam) is
  the upper bound. A different zone axis would raise the along-beam share.
- **Detectability (M6 model):** magic angle ≈ **3.98 θ_E** (≈4.3 mrad for O-K at 300 keV). The
  pipeline's **100 mrad convergence washes the anisotropy out/inverts it** (surviving −0.33) →
  a real EELS measurement needs a small collection aperture (≲2 mrad).

**Two-channel framing (state in any writeup):** Channel A = spectroscopic (anisotropic DOS →
CASTEP+OptaDOS, the ELNES *shape* change); Channel B = dynamical/channelling (→ abtem multislice,
what the microscope records). CASTEP answers "does the spectrum change"; abtem answers "what's
measured". Dipole limit: ELNES sees the *magnitude* of along-beam P, **not its sign**.

---

## 8. Open items / next actions

1. **M2(b):** read TiO₂ O-K spectrum (`runs/tio2_OK/tio2_O*.dat`), confirm endianness + shape vs
   textbook rutile O-K; this locks the recipe. Repeat for Ti-L / SrTiO₃ if desired.
2. **Capture Pb (and Sr) OTFG strings** → add to `config.OTFG` → PbTiO₃/SrTiO₃ core-hole cells
   become turnkey (rerun `build_cells.py --corehole`).
3. **M3 null test** (cubic q∥ = q⊥) — the essential artifact check before trusting M4.
4. **M4 dichroism** on `tet_Pz` (O-K primary; Ti-L for anisotropy only; Pb exploratory) + the
   `tet_Pz` q⊥ = `tet_Px` q∥ rotational-invariance cross-check.
5. **M5 calibration** across the `scan_*` ladder; **M6** maps it onto the 300 keV collection
   geometry for the detectability answer.
6. **M6b (optional, needs gpaw):** dynamical forward model — install gpaw into the abtem env
   (`pip install gpaw` + setups) to run `simulate_eels.py --edge O_K --thickness-series`.

**Edge priority + honesty:** O-K is primary (reliable single-particle DFT, strong polar
sensitivity). Ti-L₂,₃ multiplet/spin-orbit isn't captured by single-particle DFT → use it for
the *anisotropy* of the e_g/t₂g DOS, not absolute lineshape. Pb-M is exploratory.
