#!/usr/bin/env python
"""@file make_ronchigram_fig.py
@brief Ronchigram / wave-aberration evolution as a Cs-corrected-at-30-mrad scope is opened up.

Physics story for the ROUND aberration campaign: the corrector nulls 3rd order (Cs) but NOT 5th
order, so C5 is a FIXED residual (=1 mm here). As alpha opens, C5*alpha^6 grows; the operator
retunes the corrector's Cs (C3) + defocus (C1) — the (C3,C1) balance from campaign/round_sweep.tsv
— to keep the PROBE compact (the right criterion for ptycho, which recovers the phase; a flat
Scherzer chi is NOT required). This plots, per alpha: the Ronchigram (aperture wave-aberration
chi), the balanced real-space probe, chi(theta) radial, and where it breaks (probe d90 grows once
C3*alpha^4 can no longer balance C5*alpha^6, ~>90 mrad).

    ~/hyperspy-bundle/bin/python analysis/make_ronchigram_fig.py   # -> ronchigram_evolution.png
"""
import os, datetime, numpy as np, abtem, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

def week_dir(sub="figs"):
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    d = os.path.join(repo, "aberration_experiment", sub, datetime.date.today().strftime("%G-W%V"))
    os.makedirs(d, exist_ok=True); return d

LAM = 0.0196877; C5 = 1e7                       # 300 keV wavelength [Å]; fixed 5th-order residual (1 mm)
# (alpha_mrad, C3_A, C1_A) balanced probes — from campaign/round_sweep.tsv (a30 = its commented row)
PTS = [(30,1e4,-50),(50,-4e4,0),(70,-4e4,2),(90,-9e4,-160),(100,-11e4,-238),(110,-12e4,-240),(120,-13e4,-268)]

def chi_radial(C3, C1, th):
    """round wave aberration chi(theta) [rad]; abtem defocus=C1 -> C10=-C1."""
    return (2*np.pi/LAM)*(-0.5*C1*th**2 + 0.25*C3*th**4 + (1/6.)*C5*th**6)

def flat_aperture_mrad(a, C3, C1):
    """largest angle (from centre, contiguous) where |chi - running min| <= pi/4 (Scherzer-flat)."""
    th = np.linspace(0, a/1000, 400); ch = chi_radial(C3, C1, th); ch -= ch[0]; r = 0.0
    for i in range(len(th)):
        if abs(ch[i]-ch[:i+1].min()) <= np.pi/4: r = th[i]
        else: break
    return r*1000

def probe(a, C3, C1, ext=30.0, N=512):
    return np.asarray(abtem.Probe(energy=300e3, semiangle_cutoff=a, extent=ext, gpts=N,
                                  defocus=C1, aberrations={"C30":C3,"C50":C5}).build().compute().array), ext

def d90(P, ext):
    I=np.abs(P)**2; N=P.shape[0]; In=I/I.sum(); c=np.unravel_index(I.argmax(),I.shape); yy,xx=np.indices(I.shape)
    r=np.hypot(yy-c[0],xx-c[1])*(ext/N); o=np.argsort(r.ravel())
    return 2*r.ravel()[o][np.searchsorted(np.cumsum(In.ravel()[o]),0.9)]

def main():
    try: abtem.config.set({"local_diagnostics.progress_bar": False})
    except Exception: pass
    fig = plt.figure(figsize=(16,7.5), constrained_layout=True)
    gs = fig.add_gridspec(3, len(PTS), height_ratios=[1.2,1.2,1.4])
    flats=[]; d90s=[]
    for j,(a,C3,C1) in enumerate(PTS):
        N=240; k=np.linspace(-a/1000,a/1000,N); KX,KY=np.meshgrid(k,k); TH=np.hypot(KX,KY)
        ch=np.ma.masked_where(TH>a/1000, chi_radial(C3,C1,TH))
        ax=fig.add_subplot(gs[0,j])
        ax.imshow(np.angle(np.exp(1j*ch)),cmap="twilight",vmin=-np.pi,vmax=np.pi,extent=[-a,a,-a,a])
        ax.set_title(f"{a} mrad\nC3={C3/1e4:+.0f}µm",fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        P,ext=probe(a,C3,C1); dd=d90(P,ext); d90s.append(dd); flats.append(flat_aperture_mrad(a,C3,C1))
        axp=fig.add_subplot(gs[1,j]); h=60; c=np.unravel_index(np.abs(P).argmax(),P.shape)
        axp.imshow(np.abs(P[c[0]-h:c[0]+h,c[1]-h:c[1]+h]),cmap="inferno"); axp.set_xticks([]); axp.set_yticks([])
        axp.set_xlabel(f"d90={dd:.1f}Å",fontsize=9)
    axc=fig.add_subplot(gs[2,:3])
    for a,C3,C1 in PTS:
        th=np.linspace(0,a/1000,300); axc.plot(th*1000, chi_radial(C3,C1,th)-chi_radial(C3,C1,np.array([0.]))[0], label=f"{a}")
    axc.axhspan(-np.pi/4,np.pi/4,color="green",alpha=0.12); axc.set_ylim(-8,8)
    axc.set_xlabel("angle (mrad)"); axc.set_ylabel("χ (rad)"); axc.legend(fontsize=7,ncol=2,title="α mrad")
    axc.set_title("wave aberration χ(θ) — green = flat-to-π/4",fontsize=10)
    axm=fig.add_subplot(gs[2,3:]); al=[p[0] for p in PTS]
    axm.plot(al,flats,"o-",label="flat-to-π/4 aperture (mrad)"); axm.plot(al,al,"k:",lw=0.8,label="α (fully flat)")
    axm.plot(al,d90s,"s-",color="crimson",label="probe d90 (Å)"); axm.axvline(100,color="gray",ls=":")
    axm.set_xlabel("α opened to (mrad)"); axm.legend(fontsize=8)
    axm.set_title("flat aperture collapses & probe grows as α opens (C5=1mm)",fontsize=10)
    fig.suptitle("Cs-corrected@30mrad scope opened up: evolving Ronchigram (χ), balanced probe, and where it breaks",fontsize=13)
    p = os.path.join(week_dir("figs"), "ronchigram_evolution.png")
    fig.savefig(p, dpi=140); print("wrote", p)

if __name__ == "__main__":
    main()
