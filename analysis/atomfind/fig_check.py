#!/usr/bin/env python
"""@file fig_check.py
@brief Visual sanity-check figures for the blind atom finder (find.py).

Both plots consume the blind pipeline output (picked x, y, z + 3-axis error bars + species);
NO ground truth goes into the finding -- GT is only overlaid to judge it.

  1) fig_column_overlay.png -- the "is that Ti really an oxygen?" check: the sketchiest B-O
     column's raw-phase cross-section (z vs in-plane x), blind atoms drawn at their picked
     (x, z) with error bars and species colour, GT atoms as open rings by TRUE species. A
     found-Ti disc on a red (O) ring is an obvious false positive.
  2) fig_atoms_3d.png -- a small 3-D crop (fig7-style): blind atoms as spheres with +/-1
     sigma error crosses, GT model ghosted over the top; side-on + down-the-beam.

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


# ================================================================ Plot 1b: 2-row cross-sections
def _brightest_col(V, al, pos, Z, cfg, species, near_row=None, max_drow=None):
    """Brightest column of a given species in the recon depth-mean; optionally near a row."""
    dm = np.clip(V, 0, None).mean(0)
    win = align.in_window(pos, cfg)
    idx = np.where(win & (Z == species))[0]
    best = None
    for j in idx:
        rf, cf, _ = al.site_to_index(pos[j, 0], pos[j, 1], pos[j, 2])
        r, c = int(round(rf)), int(round(cf))
        if not (3 <= r < V.shape[1]-3 and 3 <= c < V.shape[2]-3):
            continue
        if near_row is not None and abs(r - near_row) > max_drow:
            continue
        b = dm[r-1:r+2, c-1:c+2].mean()
        if best is None or b > best[0]:
            best = (b, r, c)
    return best[1], best[2]


def fig_cross_sections(V, dec, dx, al, pos, Z, found, cfg, path, hw=None, hw_level="95%"):
    """Two horizontal cuts x depth, raw AND RL-deconvolved, with the blind atoms drawn on.

    hw = {"x":arr,"z":arr} per-atom CALIBRATED conformal half-widths (uncertainty.apply). If
    given, error bars show these (default the pessimistic 95% interval). If None, falls back
    to the raw model sigma -- which is the PRE-calibration CRLB and is far too small to be an
    honest error bar (median sigma_z ~0.03 A vs true |z-err| p95 ~1 A); only used if no GT.

    Row 1 -- B-site row: shows Ti-O / pure-O / Ti-O columns; Ti AND O drawn (O belongs).
    Row 2 -- A-site row: shows two Pb columns; ONLY Pb drawn. The pure-O columns physically
      interleave the Pb along this row (they sit at the same y), but a pure-Pb column carries
      no oxygen, so oxygen is deliberately not marked in the Pb panel -- the earlier version
      drew it and read as if O sat on the Pb sublattice.
    Found atoms = filled disc (species colour) + x/z 1-sigma error bars; GT = open rings."""
    from atomfind.run_atomfind import _pick_BO_column
    r_ti, c_ti, _ = _pick_BO_column(pos, Z, cfg, al, V)
    r_pb, c_pb = _brightest_col(V, al, pos, Z, cfg, 82, near_row=r_ti,
                                max_drow=int(round(3.0/dx)))
    # centre the Pb panel BETWEEN this Pb column and its nearest Pb neighbour on the same
    # row, so both columns sit symmetrically in view (a single-column centre leaves the
    # neighbour clipped by its depth-dependent lean).
    dm = np.clip(V, 0, None).mean(0)
    win = align.in_window(pos, cfg)
    pbc = []
    for j in np.where(win & (Z == 82))[0]:
        rf, cf, _ = al.site_to_index(pos[j, 0], pos[j, 1], pos[j, 2])
        rr, cc = int(round(rf)), int(round(cf))
        if abs(rr - r_pb) <= int(round(0.6/dx)) and 3 <= cc < V.shape[2]-3:
            pbc.append(cc)
    # nearest DISTINCT Pb column: 2.5-5.5 A away (Pb-Pb ~3.9 A). A smaller gap is the same
    # leaning column reappearing at a shifted pixel index, not a second column.
    pbc = np.array(sorted(set(pbc)))
    gap = np.abs(pbc - c_pb) * dx
    cand = pbc[(gap > 2.5) & (gap < 5.5)]
    c_pb2 = int(round(0.5*(c_pb + cand[np.argmin(np.abs(cand - c_pb))]))) if cand.size else c_pb
    nL = V.shape[0]
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz)); l1 = int(round(cfg.trim_z_A[1]/cfg.dz))
    strip_half_A = 0.35                                  # atoms within this of the cut row

    # (row, col, title, half-width A, species to DRAW)
    cuts = [(r_ti, c_ti, "B-site row:  Ti–O / pure-O / Ti–O columns", 5.0, (22, 8)),
            (r_pb, c_pb2, "A-site row:  two Pb columns  (no oxygen on a pure-Pb column)", 3.2, (82,))]
    fig, axes = plt.subplots(2, 2, figsize=(15, 13), sharey=True)
    n_mis_total = 0
    gr_all, gc_all, gl_all = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    for row_i, (r, c, ttl, wA, draw_sp) in enumerate(cuts):
        W = int(round(wA / dx))
        x0, x1 = max(c - W, 0), min(c + W, V.shape[2])   # clamp to the field
        for col_i, (img, kind) in enumerate([(V, "raw phase"), (dec, "RL-deconvolved")]):
            ax = axes[row_i][col_i]
            cs = img[:, r-1:r+2, x0:x1].mean(1)
            shown = cs[l0:l1]
            # extent maps PIXEL EDGES: put pixel CENTRES on integer index coordinates so
            # markers drawn at (col-c)*dx / layer*dz land exactly on their pixels.
            ax.imshow(cs, extent=[(x0-0.5-c)*dx, (x1-0.5-c)*dx,
                                  (cs.shape[0]-0.5)*cfg.dz, -0.5*cfg.dz],
                      aspect="auto", cmap="inferno", vmin=max(np.percentile(shown, 5), 0),
                      vmax=np.percentile(shown, 99.5))
            # GT rings in the strip (only the species this panel draws)
            gsel = np.where((np.abs(gr_all - r)*dx < strip_half_A + 0.25) &
                            (gc_all > x0) & (gc_all < x1))[0]
            for j in gsel:
                zz = gl_all[j]*cfg.dz
                if zz > cfg.zmax_show_A or int(Z[j]) not in draw_sp:
                    continue
                ax.scatter([(gc_all[j]-c)*dx], [zz], s=18, facecolors="none",
                           edgecolors=SPC[int(Z[j])][0], linewidths=1.1, zorder=3,
                           alpha=0.55)
            # found atoms in the strip: small discs + x AND z (1-sigma) error bars
            fsel = np.where((np.abs(found["row"] - r)*dx < strip_half_A) &
                            (found["col"] > x0) & (found["col"] < x1))[0]
            for i in fsel:
                sp = int(found["species"][i])
                zz = found["layer"][i]*cfg.dz
                if zz > cfg.zmax_show_A or sp not in draw_sp:
                    continue
                col = SPC.get(sp, ("w",))[0]
                exi = hw["x"][i] if hw is not None else found["sx_A"][i]
                ezi = hw["z"][i] if hw is not None else found["sz_A"][i]
                ax.errorbar([(found["col"][i]-c)*dx], [zz], xerr=exi, yerr=ezi,
                            fmt="o", ms=4.0, color=col, ecolor=col, elinewidth=1.0,
                            capsize=2.2, zorder=4, markeredgecolor="black",
                            markeredgewidth=0.35, alpha=0.7)
                # flag species mis-IDs (nearest GT within 1.5 A of another species)
                d = np.hypot((gc_all-found["col"][i])*dx, (gr_all-found["row"][i])*dx) \
                    + np.abs(gl_all-found["layer"][i])*cfg.dz
                jg = int(np.argmin(d))
                if Z[jg] in SPC and Z[jg] != sp and d[jg] < 1.5:
                    if col_i == 0:
                        n_mis_total += 1
                    ax.scatter([(found["col"][i]-c)*dx], [zz], s=60, facecolors="none",
                               edgecolors="yellow", linewidths=1.5, zorder=5)
            ax.set_xlim((x0-0.5-c)*dx, (x1-0.5-c)*dx)
            ax.set_ylim(cfg.zmax_show_A, -0.5*cfg.dz)
            if row_i == 1:
                ax.set_xlabel("in-plane x along the cut  (A)")
            if col_i == 0:
                ax.set_ylabel("depth z  (A)")
            ax.set_title(f"{ttl}\n{kind}", fontsize=11)
    handles = ([plt.Line2D([], [], marker="o", ls="", mfc=SPC[z][0], mec="k", ms=8,
                           label=f"found {SPC[z][1]}") for z in (82, 22, 8)] +
               [plt.Line2D([], [], marker="o", ls="", mfc="none", mec=SPC[z][0], mew=1.8,
                           ms=10, label=f"GT {SPC[z][1]}") for z in (82, 22, 8)] +
               [plt.Line2D([], [], marker="o", ls="", mfc="none", mec="yellow", mew=2,
                           ms=11, label="species mis-ID")])
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
               bbox_to_anchor=(0.5, 0.005), fontsize=9)
    barlbl = f"error bars = {hw_level} conformal interval" if hw is not None else \
             "error bars = model sigma (UNCALIBRATED — too small)"
    fig.suptitle("Blind v3 atoms on raw and PSF-deconvolved cross-sections  "
                 "(top: B-site Ti–O + pure-O columns;  bottom: A-site Pb columns, O suppressed)\n"
                 + barlbl, fontsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {os.path.basename(path)}  ({n_mis_total} mis-IDs in the drawn atoms)")


# ================================================================ Plot 1 (single column, kept)
def fig_column_overlay(V, dx, al, pos, Z, found, cfg, path, W=24, hw=None):
    from atomfind.run_atomfind import _pick_BO_column
    r, c, on = _pick_BO_column(pos, Z, cfg, al, V)
    nL = V.shape[0]; zrec = (np.arange(nL) + 0.5) * cfg.dz
    # pixel-centred extent (edges at +-0.5 px/layer) so markers at (col-c)*dx / layer*dz
    # land exactly on their pixels
    ext = [(-W-0.5)*dx, (W-0.5)*dx, (nL-0.5)*cfg.dz, -0.5*cfg.dz]
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
        exi = hw["x"][i] if hw is not None else found["sx_A"][i]
        ezi = hw["z"][i] if hw is not None else found["sz_A"][i]
        ax.errorbar([xoff], [zz], xerr=exi, yerr=ezi,
                    fmt="o", ms=7, color=col, ecolor=col, elinewidth=1.1, capsize=2, zorder=4,
                    markeredgecolor="black", markeredgewidth=0.5)
        if mism:
            n_mis += 1
            ax.scatter([xoff], [zz], s=240, facecolors="none", edgecolors="yellow",
                       linewidths=2.2, marker="o", zorder=5)

    ax.set_xlim((-W-0.5)*dx, (W-0.5)*dx); ax.set_ylim(cfg.zmax_show_A, -0.5*cfg.dz)
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
def fig_atoms_3d(V, dx, al, pos, Z, found, cfg, path, crop_r=(152, 348), z_band=(12, 54),
                 hw=None, hw_level="95%"):
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
        if hw is not None:
            S = np.c_[hw["x"][sel], hw["y"][sel], hw["z"][sel]]
        else:
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
            # PARALLEL projection + a small rotation about the column (z) axis: each
            # depth-plane of columns lands as its own separated horizontal row instead of
            # perspective-smearing into the others -- columns read as clean rows.
            p.enable_parallel_projection()
            p.camera.up = (1, 0, 0)
            p.camera.focal_point = CENTER
            p.camera.position = tuple(cc + Rc*o for cc, o in zip(CENTER, (0.38, 1.0, 0.0)))
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
    fig.suptitle(f"Blind v3 atoms ({len(P)} found) + 3-axis error bars   vs   ground-truth ghosts"
                 f"   (crop {Lx:.0f}x{Lz:.0f} A)", fontsize=14, y=0.99)
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=SPC[z][0], mec="none", ms=11,
                          label=SPC[z][1]) for z in (82, 22, 8)]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.02))
    fig.text(0.5, 0.075, f"solid = blind reconstruction (bars = {hw_level} conformal interval x,y,z)   ·   translucent = ground truth"
             "   ·   crop is DISPLAY-ONLY (atoms are found on the full field first -- no edge effects)",
             ha="center", color="0.4", fontsize=9)
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
    al = align.refine_with_atoms(al, found, pos, Z, cfg)   # fiducial map refinement
    print(f"v3 finder: {len(found)} atoms ({(found['guided']==1).sum()} guided) on {len(seeds)} columns")
    # CALIBRATED error bars for the figures: the raw model sigma is the pre-conformal CRLB
    # (median sigma_z ~0.03 A, far too small); draw the 95% conformal interval instead.
    from atomfind import deconv, validate, uncertainty
    _, m = validate.finder_report(found, pos, Z, al, cfg)
    qtab = uncertainty.calibrate(found, m, cfg, alphas=cfg.uq_alphas, min_n=cfg.uq_min_stratum)
    alpha = 0.05                                            # 95% (pessimistic) for the figures
    hw = uncertainty.apply(found, qtab, cfg, alpha)
    print(f"error bars = {int((1-alpha)*100)}% conformal interval; median z half-width "
          f"{np.median(hw['z']):.2f} A (raw model sigma_z was {np.median(found['sz_A']):.3f} A)")
    dec_path = os.path.join(cfg.out_dir, "deconvolved_vol.npy")
    if os.path.exists(dec_path):
        dec = np.load(dec_path)
    else:
        dec, _ = deconv.richardson_lucy_3d(V, psfmod.empirical_psf(cfg, dx), cfg)
    fig_cross_sections(V, dec, dx, al, pos, Z, found, cfg,
                       os.path.join(cfg.out_dir, "fig_cross_sections.png"), hw=hw)
    fig_column_overlay(V, dx, al, pos, Z, found, cfg,
                       os.path.join(cfg.out_dir, "fig_column_overlay.png"), hw=hw)
    fig_atoms_3d(V, dx, al, pos, Z, found, cfg,
                 os.path.join(cfg.out_dir, "fig_atoms_3d.png"), hw=hw)


if __name__ == "__main__":
    main()
