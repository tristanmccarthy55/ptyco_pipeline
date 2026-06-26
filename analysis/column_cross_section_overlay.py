#!/usr/bin/env python
"""Depth cross-sections (Pb & Ti columns, NL70) with the GROUND-TRUTH model overlaid.

In-plane alignment (locked by NCC): recon(row r, col c) -> GT physical
    X = X0 + c*dx ,  Y = Y0 + r*dx    (transpose, corr 0.79, zero shift)
Depth alignment: one consistent (sign, offset) for the whole volume, fitted by
correlating the recon column profile against the GT atom comb over PHYSICAL
windows only (kills the 3.9 A periodic aliases).

Overlays Pb / Ti / O atom positions on the recon cross-sections so the supervisor
can see the resolved planes land on the model and where the oxygens sit.

Usage: python column_cross_section_overlay.py  ->  ~/Desktop/column_cross_section_overlay.png
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import ase.io, abtem

D = os.path.expanduser("~/Desktop")
VASP = "/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/PTO6_STO6_18_18_labyrinthPoscar.vasp"
DZ, DX = 0.999, 20./404
X0, Y0 = 30., 10.
PB = (193, 125); TI = (153, 158); W = 22

# ---- recon
vol = np.load(D + "/NL70_new_vol.npy"); V = np.angle(vol).astype(float); V -= np.median(V, (1, 2), keepdims=True)
nL = V.shape[0]; zrec = (np.arange(nL) + 0.5) * DZ

# ---- GT atoms
a = ase.io.read(VASP); a.rotate(-90, "y", rotate_cell=True); a = abtem.orthogonalize_cell(a)
Lx, Ly, Lz = a.cell.lengths(); s = max(Lx, Ly); a.cell[0, 0] = s; a.cell[1, 1] = s
a.center(axis=0); a.center(axis=1); a.center(axis=2, vacuum=2.0)
pos = a.get_positions(); Z = a.get_atomic_numbers()

def column_atoms(yc, xc):
    Xc = X0 + xc * DX; Yc = Y0 + yc * DX
    sel = (np.abs(pos[:, 0] - Xc) < 1.3) & (np.abs(pos[:, 1] - Yc) < 0.7)
    return Xc, Yc, pos[sel], Z[sel]

def profile(yc, xc):
    p = V[:, yc-1:yc+2, xc-1:xc+2].mean((1, 2)); return p - p.mean()

# ---- joint depth registration (one sign+offset for the volume)
heavy = []
for yc, xc in (PB, TI):
    _, _, P, Zs = column_atoms(yc, xc); heavy.append(P[(Zs == 82) | (Zs == 22), 2])
profs = [profile(*c) for c in (PB, TI)]
def comb(zatoms, off, sgn):
    g = np.zeros(nL)
    for za in zatoms: g += np.exp(-0.5 * ((zrec - (sgn*za + off)) / 0.7)**2)
    return g - g.mean()
best = None
for sgn, lo, hi in [(+1, -8, 4), (-1, 66, 78)]:
    for off in np.linspace(lo, hi, 400):
        c = sum(np.dot(p, comb(z, off, sgn)) / (np.linalg.norm(p)*np.linalg.norm(comb(z, off, sgn))+1e-9)
                for p, z in zip(profs, heavy))
        if best is None or c > best[0]: best = (c, off, sgn)
corr, OFF, SGN = best
print(f"depth registration: z_recon = {SGN:+d}*z_GT + {OFF:.2f}   (joint corr {corr:.2f})")

# ---- data-driven sub-pixel in-plane calibration (NCC only locked to integer px)
ZMAX_SHOW = 66.0                       # trim exit-surface artifact (last ~4 A)
def blob_offsets(yc, xc):
    cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1); _, _, P, Zs = column_atoms(yc, xc)
    xax = (np.arange(2*W) - W) * DX; offs = []
    sel = (Zs == 82) | (Zs == 22)
    for X_at, z_at in zip(P[sel, 0], P[sel, 2]):
        lyr = int(round((SGN*z_at + OFF)/DZ - 0.5))
        if 0 <= lyr < nL and (SGN*z_at + OFF) < ZMAX_SHOW:
            w = np.clip(cs[lyr] - cs[lyr].min(), 0, None)
            if w.sum() > 0: offs.append((xax*w).sum()/w.sum() - (X_at - X0 - xc*DX))
    return offs
CAL_X = float(np.median(blob_offsets(*PB) + blob_offsets(*TI)))
print(f"sub-pixel in-plane calibration CAL_X = {CAL_X:+.3f} A ({CAL_X/DX:+.1f} px)")

# ---- figure
ext = [-W*DX, W*DX, zrec[-1], zrec[0]]
# colours chosen to stand out against inferno (dark purple -> orange -> white);
# every marker also gets a black outline (path effect) so it reads on hot blobs too.
cols = {82: ("#00f5ff", "o", "Pb"),    # bright cyan
        22: ("#39ff14", "x", "Ti"),    # neon green
        8:  ("#ff21ff", "D", "O"),     # magenta (pops on yellow/orange)
        38: ("#ffffff", "s", "Sr")}    # white
stroke = [pe.withStroke(linewidth=2.6, foreground="black")]
fig, axes = plt.subplots(1, 2, figsize=(11, 9))
for ax, (yc, xc), title in zip(axes, (PB, TI), ("Pb (lead) column", "Ti column")):
    cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
    shown = cs[zrec <= ZMAX_SHOW]
    ax.imshow(cs, extent=ext, aspect="auto", cmap="inferno",
              vmin=np.percentile(shown, 5), vmax=np.percentile(shown, 99.3))
    Xc, Yc, P, Zs = column_atoms(yc, xc)
    for zz, (col, mk, lab) in cols.items():
        m = Zs == zz
        if not m.any(): continue
        xin = P[m, 0] - Xc + CAL_X; zin = SGN * P[m, 2] + OFF
        keep = zin <= ZMAX_SHOW
        sc = ax.scatter(xin[keep], zin[keep], s=64, facecolors="none" if mk == "o" else col,
                        edgecolors=col, marker=mk, linewidths=1.7, label=lab)
        sc.set_path_effects(stroke)
    ax.set_xlim(-W*DX, W*DX); ax.set_ylim(ZMAX_SHOW, zrec[0])
    ax.set_xlabel("in-plane X (A)"); ax.set_ylabel("depth z (A)  [entrance -> exit]")
    ax.set_title(title); ax.legend(loc="upper right", fontsize=9)
fig.suptitle("NL70 reconstruction (heat) vs ground-truth model (markers) — depth cross-sections",
             fontsize=13)
fig.tight_layout()
fig.savefig(D + "/column_cross_section_overlay.png", dpi=140, bbox_inches="tight")
print("wrote column_cross_section_overlay.png")
