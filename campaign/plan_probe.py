#!/usr/bin/env python
"""Plan the aberration knobs for the ROUND alpha-sweep campaign.

Physics: the microscope carries a FIXED residual C5 = 1 mm (5th-order spherical, the
uncorrectable leftover of a C3-corrector). As you open the aperture alpha, its effect
grows as alpha^6. To keep the probe a constant ~4 Å (one unit cell) — so probe size,
overlap and scan step are held fixed and alpha is the ONLY variable — we balance it with
two knobs, exactly as an operator would:
    C3  (Cs)  = COARSE knob, quantised to 1 µm   (a negative Cs bends the mid-aperture)
    C1  (defocus) = FINE knob                     (brings the balanced disc to 4 Å)
One alpha^2 knob (C1) cannot flatten alpha^6 across the whole aperture, hence the C3 coarse
knob; see run_thin_aberration.sh. This script grid-searches (C3 @1µm, C1 fine) to land the
90%-enclosed probe diameter d90 on TARGET Å, and also reports the aberration-FREE defocus
(df_perf) that makes a matched 4 Å reference probe at each alpha.

Run LOCALLY (needs abtem); commit the .tsv it writes. run_campaign.sh reads that .tsv.
    ~/hyperspy-bundle/bin/python campaign/plan_probe.py
"""
import argparse, numpy as np, sys

def sizes(P, ext):
    I = np.abs(P)**2; dx = ext/P.shape[0]
    In = I/I.sum(); c = np.unravel_index(I.argmax(), I.shape); yy,xx = np.indices(I.shape)
    r = np.hypot(yy-c[0], xx-c[1])*dx; o = np.argsort(r.ravel()); cs = np.cumsum(In.ravel()[o])
    d = lambda q: 2*r.ravel()[o][np.searchsorted(cs, q)]
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

def plan(alphas, C5, target, out):
    import abtem
    try: abtem.config.set({"local_diagnostics.progress_bar": False})
    except Exception: pass
    EXT_S, N_S = 30.0, 512            # fast search grid (maxk 8.5 A^-1 covers <=168 mrad)
    EXT_M, N_M = 60.0, 1024           # high-quality remeasure of the winner (captures tails)
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
        # high-quality remeasure of the chosen aberrated probe + the matched perfect probe
        fw,d50,d90,d99 = sizes(build(abtem, a, c3um*1e4, C1, C5, EXT_M, N_M), EXT_M)
        # aberration-free 4 Å reference: defocus only, C3=C5=0
        dfp, dfp_d90 = best_c1(abtem, a, 0.0, 0.0, target, EXT_S, N_S)
        # window sizing: BIN4=17.5 Å, BIN2=35 Å, BIN1=70 Å real-space window
        binf = 4 if d99 < 15 else (2 if d99 < 31 else 1)
        note = "ok" if abs(d90-target) < 0.6 else "FLOOR>%.1f" % d90
        rows.append(dict(label="a%03d"%a, alpha=a, c5=C5, c3=c3um*1e4, c1=float(C1),
                         df_perf=float(dfp), bin=binf, aber_json="-",
                         d50=d50, d90=d90, d99=d99, note=note))
        print(f"  alpha={a:3d}  C3={c3um:+3d}um  C1={C1:+4.0f}A  d50={d50:.1f} d90={d90:.1f} "
              f"d99={d99:.1f}  BIN={binf}  df_perf={dfp:+.0f}  [{note}]", flush=True)
    # write tsv
    cols = ["label","alpha","c5","c3","c1","df_perf","bin","aber_json","d50","d90","d99","note"]
    with open(out, "w") as f:
        f.write("# round alpha-sweep probe plan | C5=%.3g A (1mm) | target d90=%.1f A\n" % (C5, target))
        f.write("# alpha=semiangle(mrad) c3=C30/Cs(A) c1=defocus(A) df_perf=aberr-free 4A defocus(A)\n")
        f.write("\t".join(cols)+"\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols)+"\n")
    print("wrote", out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--alphas", type=int, nargs="+", default=[30,50,70,90,100,110,120])
    ap.add_argument("--c5", type=float, default=1e7)          # 1 mm
    ap.add_argument("--target", type=float, default=4.0)      # d90 Å
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    import os
    out = a.out or os.path.join(os.path.dirname(__file__), "round_sweep.tsv")
    print(f"planning C5={a.c5:.3g} A target d90={a.target} A alphas={a.alphas}")
    plan(a.alphas, a.c5, a.target, out)
