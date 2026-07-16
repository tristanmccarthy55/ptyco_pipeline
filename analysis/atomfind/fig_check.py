#!/usr/bin/env python
"""Visual sanity-check figures for the BLIND v2 atom finder (find_atoms_v2).

Two plots, both consuming the blind pipeline output (picked x, y, z + 3-axis error bars +
species) -- NO ground truth goes into the finding; GT is only overlaid to judge it.

  1) fig_column_overlay.png -- the "is that Ti really an oxygen?" check. The brightest
     (sketchiest) B-O column's RAW-PHASE cross-section (z vs in-plane x), with the blind
     atoms drawn AT THEIR PICKED (x, z) with x/z error bars, coloured by the assigned
     species, and the ground-truth atoms as open rings coloured by their TRUE species.
     A found-Ti disc sitting on a red (O) ring is a false positive, made obvious.

  2) fig_atoms_3d.png -- a small 3-D crop, fig7-style: blind atoms as solid spheres
     (Pb blue / Ti grey / O red) with 3-D error CROSSES (+/-1 sigma in x, y AND z), and
     the ground-truth model ghosted translucent over the top. Side-on + down-the-beam.

Run:  ~/hyperspy-bundle/bin/python atomfind/fig_check.py
"""
from __future__ import annotations
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import config, align, psf as psfmod, find

# species -> (solid colour, label). Pb blue, Ti grey, O red (match the old fig7 / GT).
SPC = {82: ("#3b6fd6", "Pb"), 22: ("#8a9099", "Ti"), 8: ("#e74c3c", "O")}


# ================================================================ Plot 1
def fig_column_overlay(V, dx, al, pos, Z, found, cfg, path, W=24):
    from atomfind.run_atomfind import _pick_BO_column
    r, c, on = _pick_BO_column(pos, Z, cfg, al, V)
    nL = V.shape[0]; zrec = (np.arange(nL) + 0.5) * cfg.dz
    ext = [-W*dx, W*dx, zrec[-1], zrec[0]]
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz)); l1 = int(round(cfg.trim_z_A[1]/cfg.dz))
    xc_col = c * dx                                         # column centre (recon-frame A)

    # blind found atoms on this column (within 0.8 A in-plane of the column centre)
    fxy = np.hypot(found["col"]*dx - c*dx, found["row"]*dx - r*dx)
    fon = np.where(fxy < 0.8)[0]
    # GT atoms on this column, mapped to (x-offset, z) in the recon frame
    gr, gc, gl = al.site_to_index(pos[on, 0], pos[on, 1], pos[on, 2])

    fig, ax = plt.subplots(figsize=(7.5, 10))
    cs = V[:, r-1:r+2, c-W:c+W].mean(1)
    shown = cs[l0:l1]
    ax.imshow(cs, extent=ext, aspect="auto", cmap="inferno",
              vmin=max(np.percentile(shown, 5), 0), vmax=np.percentile(shown, 99.5))

    # GT: open rings coloured by TRUE species, at their picked (x-offset, z)
    for j, (rr, cc, ll) in enumerate(zip(gr, gc, gl)):
        zz = ll * cfg.dz
        if zz > cfg.zmax_show_A or Z[on[j]] not in SPC:
            continue
        col, _ = SPC[int(Z[on[j]])]
        ax.scatter([(cc - c)*dx], [zz], s=190, facecolors="none", edgecolors=col,
                   linewidths=2.0, marker="o", zorder=3)

    # nearest-GT species for every found atom on the column (recon frame) -> flag mis-IDs
    gr_all, gc_all, gl_all = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    n_mis = 0
    # blind found: filled discs coloured by ASSIGNED species, with x & z error bars
    for i in fon:
        zz = found["layer"][i] * cfg.dz
        if zz > cfg.zmax_show_A:
            continue
        sp = int(found["species"][i])
        col, _ = SPC.get(sp, ("w", "?"))
        xoff = found["col"][i]*dx - c*dx
        # nearest GT atom to this found atom (recon frame, A)
        d = np.hypot((gc_all-found["col"][i])*dx, (gr_all-found["row"][i])*dx) \
            + np.abs(gl_all-found["layer"][i])*cfg.dz
        jg = int(np.argmin(d))
        mism = Z[jg] in SPC and Z[jg] != sp and d[jg] < 1.5   # close but wrong species
        ax.errorbar([xoff], [zz], xerr=found["sx_A"][i], yerr=found["sz_A"][i],
                    fmt="o", ms=7, color=col, ecolor=col, elinewidth=1.1, capsize=2, zorder=4,
                    markeredgecolor="black", markeredgewidth=0.5)
        if mism:
            n_mis += 1
            ax.scatter([xoff], [zz], s=240, facecolors="none", edgecolors="yellow",
                       linewidths=2.2, marker="o", zorder=5)

    ax.set_xlim(-W*dx, W*dx); ax.set_ylim(cfg.zmax_show_A, zrec[0])
    ax.set_xlabel("in-plane x  (A)"); ax.set_ylabel("depth z  (A)")
    ax.set_title(f"Sketchiest B-O column: blind found atoms (filled disc = assigned species,\n"
                 f"x/z error bars) vs ground truth (open ring = true species).  "
                 f"Yellow ring = species mis-ID ({n_mis} on this column).",
                 fontsize=10)
    handles = ([plt.Line2D([], [], marker="o", ls="", mfc=SPC[z][0], mec="k", ms=9,
                           label=f"found {SPC[z][1]}") for z in (82, 22, 8)] +
               [plt.Line2D([], [], marker="o", ls="", mfc="none", mec=SPC[z][0], mew=2, ms=11,
                           label=f"GT {SPC[z][1]}") for z in (82, 22, 8)] +
               [plt.Line2D([], [], marker="o", ls="", mfc="none", mec="yellow", mew=2.2, ms=12,
                           label="mis-ID")])
    ax.legend(handles=handles, loc="upper right", fontsize=8, ncol=2, framealpha=0.9)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {os.path.basename(path)}  ({n_mis} species mis-IDs on the shown column)")
    return r, c


# ================================================================ Plot 2 (pyvista)
def fig_atoms_3d(V, dx, al, pos, Z, found, cfg, path, crop_r=(178, 288), z_band=(24, 46)):
    import pyvista as pv
    pv.OFF_SCREEN = True
    dz = cfg.dz
    r0, r1 = crop_r
    l0 = int(round(z_band[0]/dz)); l1 = int(round(z_band[1]/dz))
    Lx = (r1 - r0) * dx
    Lz = (l1 - l0) * dz

    def in_crop_found():
        sel = ((found["row"] >= r0) & (found["row"] < r1) &
               (found["col"] >= r0) & (found["col"] < r1) &
               (found["layer"] >= l0) & (found["layer"] < l1))
        f = found[sel]
        P = np.c_[(f["col"]-r0)*dx, (f["row"]-r0)*dx, (f["layer"]-l0)*dz]
        S = np.c_[f["sx_A"], f["sy_A"], f["sz_A"]]
        return P, f["species"], S

    gr, gc, gl = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    gsel = ((gr >= r0) & (gr < r1) & (gc >= r0) & (gc < r1) & (gl >= l0) & (gl < l1))
    GP = np.c_[(gc[gsel]-r0)*dx, (gr[gsel]-r0)*dx, (gl[gsel]-l0)*dz]
    GZ = Z[gsel]

    P, SP, SIG = in_crop_found()
    # smaller spheres so the error CROSSES read: xy sigma (~0.1 A) is smaller than an atom
    # (xy is nailed); the along-beam z sigma (~0.3 A) sticks out past the sphere.
    R_SPH = {82: 0.30, 22: 0.24, 8: 0.17}
    CENTER = (Lx/2, Lx/2, Lz/2); Rc = max(Lx, Lz)

    def add_atoms(p):
        # blind found: solid spheres + 3-axis error crosses (+/-1 sigma x,y,z) with end caps
        for zz, (col, _lab) in SPC.items():
            m = SP == zz
            if not m.any():
                continue
            sph = pv.PolyData(P[m]).glyph(geom=pv.Sphere(radius=R_SPH[zz], theta_resolution=18,
                                          phi_resolution=18), scale=False, orient=False)
            p.add_mesh(sph, color=col, smooth_shading=True, specular=0.4, specular_power=12,
                       show_scalar_bar=False)
            bars = []; cap = 0.09
            for (x, y, z), (sx, sy, sz) in zip(P[m], SIG[m]):
                bars.append(pv.Line((x-sx, y, z), (x+sx, y, z)).tube(radius=0.05, n_sides=8))
                bars.append(pv.Line((x, y-sy, z), (x, y+sy, z)).tube(radius=0.05, n_sides=8))
                bars.append(pv.Line((x, y, z-sz), (x, y, z+sz)).tube(radius=0.05, n_sides=8))
                for zc in (z-sz, z+sz):               # caps on the (dominant) z bar
                    bars.append(pv.Line((x-cap, y, zc), (x+cap, y, zc)).tube(radius=0.04, n_sides=6))
            if bars:
                p.add_mesh(pv.merge(bars), color=col, show_scalar_bar=False)

    def add_ghosts(p):
        gstyle = {82: ("#8fb8ff", 0.48), 22: ("#cfd4d8", 0.36), 8: ("#f1a9a0", 0.24), 38: ("#8fe0a0", 0.46)}
        for zz, (col, rad) in gstyle.items():
            m = GZ == zz
            if not m.any():
                continue
            sph = pv.PolyData(GP[m]).glyph(geom=pv.Sphere(radius=rad, theta_resolution=14,
                                           phi_resolution=14), scale=False, orient=False)
            p.add_mesh(sph, color=col, opacity=0.22, smooth_shading=True, show_scalar_bar=False)

    def render(cam="side", win=(1700, 620), zoom=1.7):
        p = pv.Plotter(off_screen=True, window_size=win)
        p.set_background("white")
        p.enable_depth_peeling(12)
        add_atoms(p); add_ghosts(p)
        p.add_mesh(pv.Box(bounds=(0, Lx, 0, Lx, 0, Lz)), style="wireframe",
                   color="#b0b0b0", line_width=1, opacity=0.4)
        if cam == "side":
            p.camera.up = (1, 0, 0)
            p.camera.focal_point = CENTER
            p.camera.position = tuple(cc + Rc*o for cc, o in zip(CENTER, (0.07, 1.0, 0.09)))
            p.reset_camera(); p.camera.zoom(zoom)
        else:
            p.view_xy(); p.camera.zoom(1.3)
        img = p.screenshot(return_img=True); p.close()
        return _trim(img)

    img_side = render("side")
    img_top = render("top", win=(900, 900))

    fig = plt.figure(figsize=(15, 6.4)); fig.patch.set_facecolor("white")
    axL = fig.add_axes([0.01, 0.14, 0.66, 0.74]); axR = fig.add_axes([0.68, 0.14, 0.31, 0.74])
    for ax, im, ttl in [(axL, img_side, "Side view — along the columns"),
                        (axR, img_top, "End-on — down the beam")]:
        ax.imshow(im); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(ttl, fontsize=13)
    fig.suptitle(f"Blind v2 atoms ({len(P)} found) + 3-axis error bars   vs   ground-truth ghosts"
                 f"   (crop {Lx:.0f}x{Lz:.0f} A)", fontsize=14, y=0.99)
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=SPC[z][0], mec="none", ms=11,
                          label=SPC[z][1]) for z in (82, 22, 8)]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.02))
    fig.text(0.5, 0.075, "solid = blind reconstruction (bars = +/-1 sigma x,y,z)   ·   translucent = ground truth",
             ha="center", color="0.4", fontsize=10)
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)
    print("wrote", os.path.basename(path))


def _trim(img, pad=8):
    mask = (img < 248).any(2)
    ys, xs = np.where(mask)
    if not len(ys):
        return img
    return img[max(ys.min()-pad, 0):ys.max()+pad+1, max(xs.min()-pad, 0):xs.max()+pad+1]


def main():
    cfg = config.preset("NL70_coherent")
    os.makedirs(cfg.out_dir, exist_ok=True)
    V, dx = align.load_phase(cfg)
    pos, Z = align.load_gt(cfg)
    al = align.register(V, dx, pos, Z, cfg)
    kernels = psfmod.species_kernels(cfg, dx)
    found, seeds = find.find_atoms_v3(V, cfg, dx, kernels)
    print(f"v3 finder: {len(found)} atoms ({(found['guided']==1).sum()} guided) on {len(seeds)} columns")
    fig_column_overlay(V, dx, al, pos, Z, found, cfg,
                       os.path.join(cfg.out_dir, "fig_column_overlay.png"))
    fig_atoms_3d(V, dx, al, pos, Z, found, cfg,
                 os.path.join(cfg.out_dir, "fig_atoms_3d.png"))


if __name__ == "__main__":
    main()
