#!/usr/bin/env python
"""@file make_ronchigram_fig.py
@brief Probe evolution as a Cs-corrected-at-30-mrad scope is opened up: Ronchigram, aperture phase,
       probe, phase-ramp profiles and probe size.

Physics story for the ROUND aberration campaign: the corrector nulls 3rd order (Cs) but NOT 5th
order, so C5 is a FIXED residual (1 mm here). As alpha opens, C5*theta^6 explodes (~8 rad at the
aperture edge at 50 mrad, ~1600 rad at 120); the operator retunes Cs (C3) + defocus (C1) -- the
(C3, C1) of campaign/round_sweep.tsv -- to keep the PROBE compact (the criterion ptycho needs; it
recovers the phase itself, so a flat Scherzer chi is not required).

Probe plan (campaign/plan_probe.py, read from round_sweep.tsv): d90 = 4 A first, the flattest
Ronchigram second. Where 4 A is reachable (30-70 mrad) that is the traditional recipe: Cs flat --
the realistic +1 um corrector residual at 30, retuned to the flattest step against C5 at 50/70 --
and DEFOCUS spreads the probe to 4 A (the smallest reachable probe is far smaller; bottom right).
From 90 mrad 4 A is unreachable: the probe is the smallest possible and needs C1 = -160..-268 A to
balance C3 against C5 -- a theta^2 term that winds phase from the very centre (the bullseye rings).

Rows, one column per alpha:
  1. simulated Ronchigram -- probe x thin amorphous film, far-field intensity (what the operator
     sees). The sweet spot is where the ray displacement is stationary (infinite magnification).
  2. aperture phase chi, wrapped -- ring spacing is the local phase ramp.
  3. probe |psi| at the specimen, common +-30 A scale, circle = converged d90.
Bottom:
  left  -- phase-ramp line profile: Delta x(theta) = dW/dtheta, the lateral distance (A) at which
           the ray through aperture angle theta lands. The probe is compact only where rays land
           within ~+-2 A (the 4 A target band).
  right -- probe d90 and d99 vs alpha (converged 140 A box), the 4 A target, the smallest probe
           reachable at 30-70 mrad, the 70 A BIN=1 probe window, and Ronchigram flatness: the
           NON-defocus aberration P-V across the aperture (plan_probe.nondefocus_pv), flat <= lambda/4.
  Each chi panel carries its non-defocus P-V (green = flat).

History: the pre-2026-09-11 version built probes in a 30 A box (d90 20.7/26.1 A at 110/120 mrad,
really 24.5 / >47 A) and plotted a "flat-to-pi/4 aperture" measured against a running minimum,
which counts any monotonic slope of chi as flat (it reported the full aperture flat at 50/70).

    ~/hyperspy-bundle/bin/python analysis/make_ronchigram_fig.py   # -> ronchigram_evolution.png
"""
import os, datetime, numpy as np, abtem, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter

LAM = 0.0196877; C5 = 1e7                       # 300 keV wavelength [A]; fixed 5th-order residual (1 mm)
TSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "campaign", "round_sweep.tsv")


def load_plan(path=TSV):
    """(alpha_mrad, C3_A, C1_A) per planned probe, read from the planner's tsv -- including its
    commented-out rows (a30), which are skipped by the campaigns but belong on this figure.
    Read rather than hard-coded so the figure can never drift from the probes actually simulated."""
    import re
    pts = []
    for line in open(path):
        if re.match(r"#?a\d{3}\t", line):
            f = line.lstrip("#").rstrip("\n").split("\t")
            pts.append((int(f[1]), float(f[3]), float(f[4])))
    return sorted(pts)


PTS = load_plan()


def _planner():
    """campaign/plan_probe.py, for its flatness metric -- shared so figure and plan cannot disagree."""
    import importlib.util
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "campaign", "plan_probe.py")
    s = importlib.util.spec_from_file_location("plan_probe", p); m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m); return m


PLAN = _planner()
L_BOX, N_BOX = 140.0, 2048                      # converged probe box: holds a110's d99 (51 A) with margin
WINDOW_BIN1 = 1426 * 0.0492                     # BIN=1 recon probe window (Ndp x dx) ~ 70 A
BAND = 2.0                                      # +-A: rays landing inside make a ~4 A probe


def week_dir(sub="figs"):
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    d = os.path.join(repo, "aberration_experiment", sub, datetime.date.today().strftime("%G-W%V"))
    os.makedirs(d, exist_ok=True); return d


def W(C3, C1, th):
    """round wave aberration [A]; abtem defocus=C1 -> C10=-C1."""
    return -0.5 * C1 * th**2 + 0.25 * C3 * th**4 + (1 / 6.) * C5 * th**6


def ray_dx(C3, C1, th):
    """phase ramp as a ray displacement: Delta x = dW/dtheta [A] (= lambda/2pi dchi/dtheta)."""
    return -C1 * th + C3 * th**3 + C5 * th**5


def probe(a, C3, C1, L=L_BOX, N=N_BOX):
    return np.asarray(abtem.Probe(energy=300e3, semiangle_cutoff=a, extent=L, gpts=N, defocus=C1,
                                  aberrations={"C30": C3, "C50": C5}).build().compute().array)


def enclosed(P, L, qs=(0.9, 0.99)):
    I = (np.abs(P)**2).astype(np.float64); I /= I.sum(); N = P.shape[0]
    c = np.unravel_index(I.argmax(), I.shape); yy, xx = np.indices(I.shape)
    r = np.hypot(yy - c[0], xx - c[1]) * (L / N); o = np.argsort(r.ravel()); cs = np.cumsum(I.ravel()[o])
    return [2 * r.ravel()[o][min(np.searchsorted(cs, q), cs.size - 1)] for q in qs]


def smallest(a):
    """Smallest d90 reachable at alpha with the same C5 (coarse C3 x C1 grid; small probes, small box)."""
    best = (np.inf, 0, 0)
    for c3 in np.arange(-8e4, 3.1e4, 1e4):
        for c1 in range(-80, 81, 10):
            d = enclosed(probe(a, c3, c1, 30.0, 512), 30.0, (0.9,))[0]
            best = min(best, (d, c3, c1))
    return best


def ronchigram(P, L, a, rng):
    """|FT(probe x amorphous film)|^2, cropped to +-1.15 alpha. Film: 0.5 A-correlated weak phase."""
    N = P.shape[0]
    g = gaussian_filter(rng.standard_normal((N, N)), 0.5 / (L / N)); g /= g.std()
    R = np.abs(np.fft.fftshift(np.fft.fft2(P * np.exp(0.6j * g))))**2
    px = LAM / L * 1e3                                    # mrad per pixel
    h = min(int(round(1.15 * a / px)), N // 2 - 1); c = N // 2
    return np.sqrt(R[c - h:c + h, c - h:c + h]), h * px


def main():
    try: abtem.config.set({"local_diagnostics.progress_bar": False})
    except Exception: pass
    rng = np.random.default_rng(7)
    fig = plt.figure(figsize=(17, 12.5), constrained_layout=True)
    gs = fig.add_gridspec(4, len(PTS), height_ratios=[1, 1, 1, 1.45])
    d90s, d99s = [], []
    cols = plt.cm.viridis(np.linspace(0, 0.92, len(PTS)))
    for j, (a, C3, C1) in enumerate(PTS):
        P = probe(a, C3, C1)
        d90, d99 = enclosed(P, L_BOX); d90s.append(d90); d99s.append(d99)

        R, lim = ronchigram(P, L_BOX, a, rng)
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(R, cmap="gray", extent=[-lim, lim, -lim, lim], vmax=np.percentile(R, 99.5))
        ax.add_patch(Circle((0, 0), a, fill=False, ec="#ffcc00", lw=0.6, ls="--"))
        ax.set_title(f"{a} mrad\nC3 = {C3/1e4:+.0f} µm,  C1 = {C1:+.0f} Å", fontsize=9.5)
        ax.set_xticks([]); ax.set_yticks([])
        if j == 0: ax.set_ylabel("Ronchigram\n(amorphous film)", fontsize=10)

        n = 240; k = np.linspace(-a / 1000, a / 1000, n); KX, KY = np.meshgrid(k, k); TH = np.hypot(KX, KY)
        ch = np.ma.masked_where(TH > a / 1000, 2 * np.pi / LAM * W(C3, C1, TH))
        ax = fig.add_subplot(gs[1, j])
        ax.imshow(np.angle(np.exp(1j * ch)), cmap="twilight", vmin=-np.pi, vmax=np.pi, extent=[-a, a, -a, a])
        ax.set_xticks([]); ax.set_yticks([])
        pv = PLAN.nondefocus_pv(a, C3, C5)
        ax.set_xlabel(f"non-defocus P-V {pv:.1f} rad", fontsize=9,
                      color="tab:green" if pv <= PLAN.FLAT_TOL else "0.2")
        if j == 0: ax.set_ylabel("aperture phase χ\n(wrapped)", fontsize=10)

        ax = fig.add_subplot(gs[2, j]); h = int(round(30 / (L_BOX / N_BOX))); c = np.unravel_index(np.abs(P).argmax(), P.shape)
        crop = np.abs(P[c[0] - h:c[0] + h, c[1] - h:c[1] + h])
        ax.imshow(np.sqrt(crop), cmap="inferno", extent=[-30, 30, -30, 30])
        ax.add_patch(Circle((0, 0), d90 / 2, fill=False, ec="cyan", lw=0.8))
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(f"d90 = {d90:.1f} Å", fontsize=9.5)
        if j == 0: ax.set_ylabel("probe |ψ|  (±30 Å)\ncircle = d90", fontsize=10)

    # --- phase-ramp line profiles across the aperture diameter -------------------------
    axr = fig.add_subplot(gs[3, :4])
    for (a, C3, C1), col in zip(PTS, cols):
        th = np.linspace(-a / 1000, a / 1000, 801)
        axr.plot(th * 1000, ray_dx(C3, C1, th), color=col, lw=1.6, label=f"{a}")
    axr.axhspan(-BAND, BAND, color="green", alpha=0.13, lw=0)
    axr.axhline(0, color="0.5", lw=0.6)
    axr.set_ylim(-32, 32); axr.set_xlim(-125, 125)
    axr.set_xlabel("aperture angle θ along a diameter (mrad)")
    axr.set_ylabel("ray lands at Δx = ∂W/∂θ  (Å)")
    axr.set_title("phase-ramp line profile: where each part of the aperture throws its rays "
                  f"(green = ±{BAND:.0f} Å, a 4 Å footprint; defocus = a straight line)", fontsize=10)
    axr.legend(title="α (mrad)", fontsize=8, ncol=2, loc="upper center")
    axr.text(118, 29, "120 mrad reaches\n±56 Å at the edge", ha="right", va="top", fontsize=8, color=cols[-1])

    # --- probe size vs alpha ----------------------------------------------------------
    axs = fig.add_subplot(gs[3, 4:])
    al = np.array([p[0] for p in PTS])
    small = {a: smallest(a) for a in (30, 50, 70)}
    axs.semilogy(al, d90s, "s-", color="crimson", label="probe d90 (balanced, as simulated)")
    axs.semilogy(al, d99s, "^--", color="crimson", alpha=0.55, label="probe d99")
    axs.semilogy(list(small), [small[a][0] for a in small], "o", mfc="white", mec="crimson",
                 label="smallest d90 reachable (30–70 mrad)")
    axs.axhline(4.0, color="0.3", ls=":", lw=1); axs.text(122, 4.3, "4 Å target", ha="right", fontsize=8)
    axs.axhline(WINDOW_BIN1, color="k", ls="-.", lw=0.9)
    axs.text(84, WINDOW_BIN1 * 1.08, f"BIN=1 probe window {WINDOW_BIN1:.0f} Å", ha="left", fontsize=8)
    axs.set_ylim(0.5, 150); axs.set_xlim(25, 125)
    axs.set_xlabel("α opened to (mrad)"); axs.set_ylabel("probe diameter (Å, log)")
    # Ronchigram flatness = the NON-defocus aberration (defocus only changes the shadow magnification
    # uniformly), the quantity the planner optimises. (An "aperture within +-2 A" measure would
    # penalise the intended defocus spreading itself.)
    ax2 = axs.twinx()
    ax2.semilogy(al, [PLAN.nondefocus_pv(*p[:2], C5) for p in PTS], "o-", color="tab:green",
                 label="non-defocus aberration P-V (Ronchigram flatness)")
    ax2.axhline(PLAN.FLAT_TOL, color="tab:green", ls=":", lw=1)
    ax2.text(27, PLAN.FLAT_TOL * 0.62, "flat: λ/4", color="tab:green", fontsize=8)
    ax2.set_ylim(0.1, 400); ax2.set_ylabel("non-defocus P-V across the aperture (rad, log)", color="tab:green")
    h1, l1 = axs.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    axs.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left", bbox_to_anchor=(0.01, 0.85), framealpha=0.85)
    axs.set_title("probe grows once C3 can no longer hold C5 (C5 = 1 mm)", fontsize=10)
    for a in small:
        print(f"a{a}: smallest reachable d90 {small[a][0]:.2f} A at C3 {small[a][1]/1e4:+.0f} um, C1 {small[a][2]:+d} A")
    for (a, C3, C1), d90, d99 in zip(PTS, d90s, d99s):
        print(f"a{a}: d90 {d90:.1f} A, d99 {d99:.1f} A, non-defocus P-V {PLAN.nondefocus_pv(a, C3, C5):.2f} rad")

    fig.suptitle("Cs-corrected@30 mrad scope opened up — Ronchigram, aperture phase, probe, "
                 "and where it breaks", fontsize=14)
    fig.text(0.5, -0.005,
             "Probe plan (plan_probe.py): d90 = 4 Å first, flattest Ronchigram second. 30–70 mrad: Cs flat "
             "(realistic +1 µm residual at 30; retuned against C5 at 50/70) and DEFOCUS enlarges the probe to 4 Å — "
             "the traditional recipe. ≥ 90 mrad: 4 Å is unreachable, the probe is the smallest possible, and C3 "
             "and C1 must both fight C5 (C1 = −160 … −268 Å winds rings from the centre).",
             ha="center", va="top", fontsize=9.5, wrap=True)
    p = os.path.join(week_dir("figs"), "ronchigram_evolution.png")
    fig.savefig(p, dpi=130, bbox_inches="tight"); print("wrote", p)


if __name__ == "__main__":
    main()
