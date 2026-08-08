# Targeted literature review — EELS/ELNES, atomic displacements, and the DFT pipeline

Scope: the background needed for the deliverable *"following a targeted literature review, CASTEP
will be used to run DFT calculations on a simple unit cell … investigate how atomic displacements
affect EELS spectra."* Focus: (i) what ELNES measures, (ii) the CASTEP+OptaDOS method we use,
(iii) reproducible small-system benchmarks, (iv) prior results on displacement→ELNES in titanates,
(v) orientation/anisotropy (the along-beam angle), and (vi) where our study is novel.

---

## 1. What EELS/ELNES measures, and why displacements matter

In core-loss EELS the beam electron ionises an atomic core level; the **energy-loss near-edge
structure (ELNES)**, the first ~30–40 eV above the edge onset, mirrors the symmetry- and
site-projected **unoccupied density of states** reached by a dipole transition from that core
level. It is therefore a local probe of coordination, bonding and — crucially here — of small
**atomic displacements** that change the local symmetry and hybridisation. In an anisotropic
crystal ELNES also depends on the **momentum-transfer direction q** (≈ the beam for a small
on-axis aperture), which is the knob that lets it sense *directional* structure (§5).

## 2. DFT-ELNES methodology (the pipeline we run)

We use the plane-wave pseudopotential **core-hole** approach in **CASTEP** with **OptaDOS**
post-processing:
- **Core-hole plane-wave method:** Gao, Pickard, Payne, Yuan & Zhu, *Core-level spectroscopy
  calculation and the plane-wave pseudopotential method* — the basis for CASTEP ELNES (excited
  atom carries a pseudopotential with a reduced core occupancy; the conduction states relax
  around the hole). [ref](https://www.researchgate.net/publication/51547996)
- **OptaDOS:** Morris, Nicholls, Pickard & Yates, *OptaDOS: A tool for obtaining density of
  states, core-level and optical spectra …*, Comput. Phys. Commun. **185**, 1477 (2014); and
  Nicholls, Morris, Pickard & Yates, *OptaDOS — a new tool for EELS calculations*, J. Phys.
  Conf. Ser. **371**, 012062 (2012). Adaptive/linear broadening + the core-loss task we use.
  [CPC](https://www.sciencedirect.com/science/article/abs/pii/S0010465514000460) ·
  [JPCS](https://iopscience.iop.org/article/10.1088/1742-6596/371/1/012062)
- **Practical CASTEP-ELNES workflow + parser:** the Mizoguchi-group tutorial and
  `castep_elnes_parser` document the exact `{1s1}`-style core-hole OTFG syntax and multi-atom
  output we reverse-engineered (see HANDOVER §3). [tutorial](https://www.edge.iis.u-tokyo.ac.jp/CASTEP-ELNES-manual-2020-English.htm)
  · [parser](https://github.com/nmdl-mizo/castep_elnes_parser) · [CASTEP core-loss](https://www.castep.org/features/capabilities/core-loss-and-optical-spectroscopies)

**Known method limits (shape our edge choice):** single-particle DFT omits the **2p spin-orbit +
multiplet** physics that dominates the **Ti L₂,₃** lineshape — TDDFT/multiplet or multichannel
multiple-scattering are needed for quantitative Ti-L (§3). Hence we lead with **O-K** (dipole
1s→2p, well described at this level) and treat Ti-L for *anisotropy trends* only.

## 3. Reproducible small-system benchmarks (pipeline validation)

- **Rutile/anatase TiO₂ — Ti L₂,₃ and O K:** the standard titanium-oxide ELNES benchmark.
  Experiment-vs-calculation for the Ti-L fine structure: *New fine structures resolved at the
  ELNES Ti-L₂,₃ edge spectra of anatase and rutile*, Ultramicroscopy (2010) — and it explicitly
  documents that **single-particle DFT struggles with Ti-L** while O/Ti-K are well reproduced.
  [Ti-L exp/calc](https://www.sciencedirect.com/science/article/abs/pii/S0304399110000793) ·
  [rutile L₂,₃ calc](https://www.sciencedirect.com/science/article/abs/pii/S003960281000484X) ·
  [multi-code Ti K benchmark](https://arxiv.org/pdf/2303.17089)
  → **We already reproduced the rutile O-K edge** (M2b): t₂g/e_g split ~2.7 eV with e_g>t₂g and
  the broad O 2p–Ti 4sp band — the textbook shape. Pipeline validated on a small system.
- **SrTiO₃ (recommended second benchmark):** cubic 5-atom perovskite titanate — the closest
  simple analogue to PbTiO₃, with extensively documented O-K/Ti-L ELNES (e.g. it is the
  paraelectric member studied in the superlattice work below). Doubles as the **non-polar
  reference** for the displacement study (STO cubic / no off-centering vs PTO tetragonal /
  displaced). Needs the Sr OTFG string (one grep from an STO run).
- **OptaDOS `Si2_CORE` example** (ships with the code we built): the most bulletproof pure-code
  reproducibility check — identical inputs → the published Si L₂,₃ result.

## 4. Atomic displacements → ELNES: prior results (the core objective)

- **PbTiO₃/SrTiO₃ superlattices — our exact material:** Torres-Pardo, Gloter, Zubko, Jecklin,
  Lichtensteiger, Colliex, Triscone & Stéphan, *Spectroscopic mapping of local structural
  distortions in ferroelectric PbTiO₃/SrTiO₃ superlattices at the unit-cell scale*, Phys. Rev. B
  **84**, 220102(R) (2011). Unit-cell-resolved Ti-L/O-K ELNES tracks the tetragonal distortion
  and the reduced Ti 3d–Pb 6sp–O 2p hybridisation across the interfaces.
  [PRB](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.220102) ·
  [arXiv:1112.4953](https://arxiv.org/abs/1112.4953)
- **O-K probes the off-centering (mechanism):** Bugnet et al., high-energy-resolution EELS of
  ferroelectric vs paraelectric BaTiO₃ — the **lowest-energy O-K fine structure broadens/becomes
  asymmetric as the Ti⁴⁺ off-centering grows**, via the O-site symmetry's effect on the core-hole
  screening. Establishes that **O-K is as sensitive as Ti-L to the ferroelectric displacement** —
  the direct precedent for using O-K to read out the polar distortion.
  [Wiley/EMC](https://onlinelibrary.wiley.com/doi/abs/10.1002/9783527808465.EMC2016.6645) ·
  [temp-dependent](https://www.researchgate.net/publication/289528630)
- **Atomic-resolution polarisation via EELS:** BaTiO₃/manganite interface studies map the polar
  state atom-column-by-column, confirming ELNES fingerprints of the displacement at STEM
  resolution. [ref](https://www.researchgate.net/publication/311522784)

**Takeaway:** it is *established* that ferroelectric off-centering imprints on Ti-L and O-K ELNES,
and that O-K is a reliable, DFT-tractable probe of it. This is exactly what our **M5 scan**
(ELNES vs displacement amplitude) quantifies from first principles.

## 5. Orientation / anisotropy dependence (the along-beam angle)

- **Orientation-dependent ELNES theory/simulation:** *The orientation-dependent simulation of
  ELNES* and *Orientation dependence of ionization edges in EELS* — ELNES depends on q relative
  to the crystal axes; simulation must project the unoccupied DOS onto q and integrate over the
  collection/convergence aperture.
  [sim](https://www.sciencedirect.com/science/article/abs/pii/S0304399199001680) ·
  [orientation](https://www.sciencedirect.com/science/article/abs/pii/S030439910000125X)
- **O-K anisotropy from band structure (direct analogue of our method):** *Anisotropy and
  collection-angle dependence of the oxygen K ELNES in V₂O₅* — computes the q-resolved O-K and its
  aperture averaging, precisely the calculation we do with OptaDOS `core_qdir`.
  [V₂O₅](https://www.sciencedirect.com/science/article/abs/pii/S0968432803000313)
- **Magic angle:** *ELNES at magic-angle conditions* — at a particular collection/convergence
  combination the anisotropic term cancels; strongly voltage-dependent. This is our M6 concern
  (the pipeline's 100 mrad convergence sits far past the magic angle, ≈4·θ_E, and averages the
  anisotropy out). [MAC](https://www.sciencedirect.com/science/article/abs/pii/S030439910600115X)

## 6. Positioning — what is known vs. our contribution

Known: (a) DFT core-hole ELNES reproduces O-K/Ti-K edges of Ti–O compounds; (b) ferroelectric
off-centering imprints on Ti-L/O-K ELNES; (c) that imprint has been mapped at unit-cell scale in
**our exact PTO/STO superlattice** (Torres-Pardo 2011); (d) ELNES is orientation-dependent, with a
known magic-angle averaging.

**Gap we address:** all the distortion-mapping above reads the **projected / in-plane** structure —
the same information projected STEM imaging (and our ptychography) already recovers. The
**component of the polar displacement *along the beam*** is invisible to projection. We ask, from
first principles, whether the **q-resolved ELNES dichroism (q∥c vs q⊥c)** carries that along-beam
component — turning ELNES into a route to the *full 3-D* polarisation of a vortex. (Honest limit:
in the dipole approximation ELNES sees the *magnitude*, not the *sign*, of along-beam P.)

## 7. How the review shapes the runs

- **Edge priority:** O-K primary (DFT-reliable *and* a validated off-centering probe, Bugnet;
  reproduced at M2b); Ti-L for anisotropy trends only (multiplet caveat, §3).
- **Benchmarks:** rutile TiO₂ O-K ✅ done; add **SrTiO₃** (analogue + paraelectric reference) and
  optionally the OptaDOS `Si2_CORE` exact-reproducibility check.
- **Displacement study (the objective):** M5 scan = ELNES vs the Ti off-centering amplitude —
  the first-principles version of Bugnet's "broadening ∝ off-centering", resolved by q.
- **Along-beam (our angle):** M4 q∥c vs q⊥ dichroism + M6 magic-angle/aperture averaging.

---

### Reference list (grouped; links above)
Method: Gao et al. (core-hole PW); Morris/Nicholls et al. OptaDOS (CPC 185, 1477, 2014; JPCS 371,
012062, 2012); Mizoguchi CASTEP-ELNES tutorial. · Benchmarks: TiO₂ Ti-L exp/calc (Ultramicroscopy
2010); rutile L₂,₃ calc; Ti-K multi-code (arXiv 2303.17089). · Displacement→ELNES: Torres-Pardo et
al. PRB 84, 220102(R) (2011); Bugnet et al. (BaTiO₃ O-K); BaTiO₃/manganite interface EELS. ·
Anisotropy: orientation-dependent ELNES sim; V₂O₅ O-K anisotropy; ELNES at magic-angle conditions.
