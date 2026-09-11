#!/usr/bin/env python
"""Plan the aberration knobs for the ROUND alpha-sweep campaign.

Physics: the microscope carries a FIXED residual C5 = 1 mm (5th-order spherical, the
uncorrectable leftover of a C3-corrector). As you open the aperture alpha, its effect
grows as alpha^6. To keep the probe a constant ~4 Å (one unit cell) — so probe size,
overlap and scan step are held fixed and alpha is the ONLY variable — the operator has two
knobs: C3 (Cs, quantised to 1 µm) and C1 (defocus, fine).

Objective, in priority order (as an operator sets up a Cs-corrected scope):
  1. d90 = TARGET (4 Å).
  2. the flattest Ronchigram that allows it. Flatness = peak-to-valley of the NON-DEFOCUS
     aberration phase across the aperture (C3 theta^4 + C5 theta^6 after removing its best-fit
     theta^2 part): a pure-defocus Ronchigram is a uniform-magnification shadow, which is what an
     operator calls flat. "Flat" = <= pi/2 (lambda/4, Rayleigh).
Two regimes follow:
  * FREE (4 Å reachable; low alpha): the traditional recipe. Keep Cs at a realistic corrector
    residual (CS_RESID, +1 µm) if that is already flat; otherwise retune Cs to the flattest
    setting. Then spread the probe to 4 Å with DEFOCUS alone.
  * FLOOR (4 Å unreachable; high alpha): no freedom -- C3 and C1 must both fight C5 for the
    smallest probe. Grid-search (C3 @1µm around the ray-edge cancel, C1 fine), rank by d90.
The pre-2026-09-11 planner used the FLOOR ranking everywhere with a least-defocus tie-break,
which at 50/70 mrad made Cs (-4 µm, C1 ~ 0) do the enlarging -- the opposite of the recipe.
Also reports the aberration-FREE defocus (df_perf) for a matched 4 Å reference probe.

Run LOCALLY (needs abtem); commit the .tsv it writes. run_campaign.sh / run_thin_atomfind.sh
read it by leading column position; flat_pv is appended last.
    ~/hyperspy-bundle/bin/python campaign/plan_probe.py
"""
import argparse, numpy as np, sys

LAM = 0.0196877              # 300 keV electron wavelength [Å]
CS_RESID_UM = 1              # realistic Cs-corrector residual (µm): kept if it is already flat
FLAT_TOL = np.pi / 2         # non-defocus phase P-V (rad) counted as a flat Ronchigram (lambda/4)
BIN1_WINDOW = 1426 * 0.0492  # BIN=1 recon probe window (Å); a probe with d99 beyond it cannot fit
# rows written commented-out, with the reason (kept here so a re-plan does not lose it)
SKIP = {30: ["a030 SKIPPED: at 30 mrad the BF disc fills only ~15% of the 200 mrad detector -> the recon is",
             "noise-dominated even with the KNOWN probe (lattice peak/bg 2.6 vs 100+ at >=50), and NL=1 carries",
             "no depth info. Fixing it needs a low-alpha detector crop, not a probe change. Uncomment to include."]}


def nondefocus_pv(alpha, C3, C5):
    """Peak-to-valley (rad) of the aberration phase across the aperture after removing its
    best-fit piston + defocus (theta^2) part, area-weighted. Independent of the C1 applied."""
    th = np.linspace(0, alpha / 1000.0, 600); sw = np.sqrt(th)
    W = C3 * th**4 / 4 + C5 * th**6 / 6
    A = np.c_[np.ones_like(th), th**2]
    coef = np.linalg.lstsq(A * sw[:, None], W * sw, rcond=None)[0]
    r = W - A @ coef
    return 2 * np.pi / LAM * (r.max() - r.min())


def flat_c3_um(alpha, C5):
    """The FREE-regime Cs: the realistic residual if already flat, else the flattest 1 µm step."""
    if nondefocus_pv(alpha, CS_RESID_UM * 1e4, C5) <= FLAT_TOL:
        return CS_RESID_UM
    return min(range(-20, 6), key=lambda c: nondefocus_pv(alpha, c * 1e4, C5))

def sizes(P, ext):
    I = (np.abs(P)**2).astype(np.float64); dx = ext/P.shape[0]   # float64: a float32 cumsum over
    In = I/I.sum(); c = np.unravel_index(I.argmax(), I.shape); yy,xx = np.indices(I.shape)  # 4M px never reaches 0.99
    r = np.hypot(yy-c[0], xx-c[1])*dx; o = np.argsort(r.ravel()); cs = np.cumsum(In.ravel()[o])
    d = lambda q: 2*r.ravel()[o][min(np.searchsorted(cs, q), cs.size - 1)]
    fwhm = 2*np.sqrt((I >= 0.5*I.max()).sum()*dx*dx/np.pi)
    return fwhm, d(0.5), d(0.9), d(0.99)

def build(abtem, alpha, C3, C1, C5, ext, N):
    ab = {"C30": C3, "C50": C5} if C5 else {"C30": C3}
    return np.asarray(abtem.Probe(energy=300e3, semiangle_cutoff=alpha, extent=ext,
                                  gpts=N, defocus=C1, aberrations=ab).build().compute().array)

def best_c1(abtem, alpha, C3, C5, target, ext, N):
    """1-D fine search over C1 (defocus) minimising |d90 - target|; tie-break small |C1|."""
    best = None
    grid = list(range(-260, 261, 10))
    for _pass in range(2):
        for C1 in grid:
            _,_,d90,_ = sizes(build(abtem, alpha, C3, C1, C5, ext, N), ext)
            key = (abs(d90-target), abs(C1))
            if best is None or key < best[0]:
                best = (key, C1, d90)
        c = best[1]; grid = list(range(c-8, c+9, 2))          # refine ±8 Å @2 Å
    return best[1], best[2]                                    # C1, d90

def plan(alphas, C5, target, thick, out):
    import abtem
    try: abtem.config.set({"local_diagnostics.progress_bar": False})
    except Exception: pass
    # depth res = LAM/alpha^2
    EXT_S, N_S = 30.0, 512            # fast search grid (maxk 8.5 A^-1 covers <=168 mrad)
    # remeasure box: 140 Å (was 60) -- a110's d99 is 51 Å and a120's ~105 Å; a 60 Å box wraps the
    # tails and under-reports both (a120 d90 read 46.8 Å, converged ~60 Å)
    EXT_M, N_M = 140.0, 2048
    rows = []
    for a in alphas:
        edge = round(-C5*(a/1000.0)**2 / 1e4)                 # ray-edge C3 cancel, in µm
        cands = sorted(set(range(edge-2, edge+3)) | {0})      # C3 grid @1µm, incl. 0
        pick = None
        for c3um in cands:
            C3 = c3um*1e4
            C1, d90 = best_c1(abtem, a, C3, C5, target, EXT_S, N_S)
            # bucket d90 error to 0.2 Å so equally-good 4 Å balances are ranked by SIMPLEST
            # knobs: least defocus (focal plane centred) then least Cs. Floor cases (d90>>4)
            # keep distinct buckets so the true minimum still wins.
            key = (int(abs(d90-target)/0.2), abs(C1), abs(c3um))
            if pick is None or key < pick[0]:
                pick = (key, c3um, C1, d90)
        c3um, C1 = pick[1], pick[2]
        regime = "floor"
        if abs(pick[3] - target) < 0.6:
            # FREE regime: 4 Å is reachable, so use the traditional recipe -- Cs flat (realistic
            # residual if possible), DEFOCUS spreads the probe. Falls back to the pick if the flat
            # Cs cannot reach the target (never seen, but the target is the primary objective).
            fc3 = flat_c3_um(a, C5)
            fC1, fd90 = best_c1(abtem, a, fc3 * 1e4, C5, target, EXT_S, N_S)
            if abs(fd90 - target) < 0.6:
                c3um, C1, regime = fc3, fC1, "free"
        pv = nondefocus_pv(a, c3um * 1e4, C5)
        # high-quality remeasure of the chosen aberrated probe + the matched perfect probe
        fw,d50,d90,d99 = sizes(build(abtem, a, c3um*1e4, C1, C5, EXT_M, N_M), EXT_M)
        # aberration-free 4 Å reference: defocus only, C3=C5=0
        dfp, dfp_d90 = best_c1(abtem, a, 0.0, 0.0, target, EXT_S, N_S)
        # window sizing: BIN4=17.5 Å, BIN2=35 Å, BIN1=70 Å real-space window
        binf = 4 if d99 < 15 else (2 if d99 < 31 else 1)
        if abs(d90-target) < 0.6:
            note = "ok:" + ("flat" if pv <= FLAT_TOL else "flattest") if regime == "free" else "ok"
        else:
            note = "FLOOR>%.1f" % d90
        if d99 > BIN1_WINDOW:
            note += "|WINDOW:d99>%.0fA" % BIN1_WINDOW
        # Depth sampling: recon at NL layers = Nyquist of the depth resolution LAM/alpha^2
        # (slice = depthres/2). The sim slab is a fixed 0.9 A (< the ~1.94 A plane spacing, run_campaign
        # SLICE default) so the GROUND TRUTH always resolves the planes and only alpha + NL decide
        # what the recon recovers -- low alpha CAN'T see the interatomic spacing, high alpha can.
        depthres = LAM/(a/1000.0)**2
        nl = max(1, round(thick / (depthres/2.0)))
        rows.append(dict(label="a%03d"%a, alpha=a, c5=C5, c3=c3um*1e4, c1=float(C1),
                         df_perf=float(dfp), bin=binf, nl=nl, aber_json="-",
                         d50=d50, d90=d90, d99=d99, note=note, flat_pv=round(pv, 2)))
        print(f"  alpha={a:3d}  C3={c3um:+3d}um  C1={C1:+4.0f}A  d90={d90:.1f} d99={d99:.1f} BIN={binf}  "
              f"non-defocus P-V {pv:.2f} rad  [{regime}; {note}]", flush=True)
    # write tsv (flat_pv LAST: the shell readers take leading columns by position)
    cols = ["label","alpha","c5","c3","c1","df_perf","bin","nl","aber_json","d50","d90","d99","note","flat_pv"]
    with open(out, "w") as f:
        f.write("# round alpha-sweep probe plan | C5=%.3g A (1mm) | target d90=%.1f A\n" % (C5, target))
        f.write("# alpha=semiangle(mrad) c3=C30/Cs(A) c1=defocus(A) df_perf=aberr-free 4A defocus(A)\n")
        f.write("# low alpha (note ok:flat/ok:flattest): Cs flat (residual +%d um if flat), DEFOCUS spreads to 4 A;\n"
                "# FLOOR: 4 A unreachable, C3+C1 fight C5. flat_pv = non-defocus aberration P-V (rad), flat <= pi/2\n"
                % CS_RESID_UM)
        f.write("\t".join(cols)+"\n")
        for r in rows:
            line = "\t".join(str(r[c]) for c in cols)
            if r["alpha"] in SKIP:
                f.write("".join(f"# {s}\n" for s in SKIP[r["alpha"]]) + "#" + line + "\n")
            else:
                f.write(line + "\n")
    print("wrote", out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--alphas", type=int, nargs="+", default=[30,50,70,90,100,110,120])
    ap.add_argument("--c5", type=float, default=1e7)          # 1 mm
    ap.add_argument("--target", type=float, default=4.0)      # d90 Å
    ap.add_argument("--thick", type=float, default=11.715)    # slab beam thickness [Å] (3 cells) for NL
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    import os
    out = a.out or os.path.join(os.path.dirname(__file__), "round_sweep.tsv")
    print(f"planning C5={a.c5:.3g} A target d90={a.target} A thick={a.thick} A alphas={a.alphas}")
    plan(a.alphas, a.c5, a.target, a.thick, out)
