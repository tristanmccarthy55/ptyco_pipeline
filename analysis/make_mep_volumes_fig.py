#!/usr/bin/env python
"""@file make_mep_volumes_fig.py
@brief MEP phase-volume figure for the thin aberration campaign: in-plane + two x-z sections per alpha.

Per convergence angle, from the full-box KNOWN-probe lab reconstruction (recon_af_a<A>_lab_NL<NL>):
  (1) in-plane phase summed over the atomic slab, scan field only;
  (2) x-z section along a row of A-site columns -- Pb on the AO planes, O columns on the BO2
      planes 1.95 A away, so resolving the two needs depth resolution of that order;
  (3) x-z section along a row of B-site columns -- Ti and apical O ALTERNATE every 1.95 A along
      the SAME column, the hardest depth target (expected to resolve only at ~100 mrad).
The volume is loaded and cropped by atomfind's own code (align.load_phase + crop_to_fov with the
thin preset), so it is in exactly the frame of found_atoms.npy: atomfind's blind atoms are
overlaid by species, and the ground-truth AO/BO2 plane depths are drawn through atomfind's fitted
depth map (z_recon = SGN*z_GT + OFF, from report.json). The 4 A z-vacuum bands are shaded -- the
entrance/exit artefacts should sit there, not on the atomic planes.

  python analysis/make_mep_volumes_fig.py --recons ~/Desktop/thin_ab_af2 \
      --atomfind ~/Desktop/thin_ab_af2/out --gt ~/Desktop/thin_ab_af/gtdata

Alphas are discovered from recon_af_a<A>_lab_NL<NL>/**/*_recons.h5; an alpha with no atomfind run
is drawn without overlays. Figure defaults to aberration_experiment/figs/<ISO-week>/mep_volumes.png.
"""
import argparse, datetime, glob, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from atomfind import config, align                       # noqa: E402  (atomfind's own frame)

BOXZ, ZVAC = 27.525, 4.0                                  # full box and z-vacuum band [A]
SP = {82: ("Pb", "#00e5ff", "o"), 22: ("Ti", "#7CFC00", "s"), 8: ("O", "#ff4fd8", "D")}


def week_dir(sub):
    repo = os.path.dirname(HERE)
    d = os.path.join(repo, "aberration_experiment", sub, datetime.date.today().strftime("%G-W%V"))
    os.makedirs(d, exist_ok=True)
    return d


def discover(recons):
    out = []
    for d in glob.glob(os.path.join(recons, "recon_af_a*_lab_NL*")):
        m = re.search(r"recon_af_a(\d+)_lab_NL(\d+)$", d.rstrip("/"))
        h = glob.glob(os.path.join(d, "**", "*_recons.h5"), recursive=True)
        if m and h:
            out.append((int(m.group(1)), int(m.group(2)), sorted(h, key=os.path.getmtime)[-1]))
    return sorted(out)


def load(h5, nl, gt):
    """The volume exactly as atomfind sees it: thin preset, per-alpha dz, scan-field crop."""
    config.set_data_dir(gt)
    cfg = config.preset("thin")
    cfg.recon_vol, cfg.dz = h5, BOXZ / nl
    cfg = cfg.resolve()
    V, dx = align.load_phase(cfg)
    V = align.crop_to_fov(V, dx, cfg)
    pos, Z = align.load_gt(cfg)
    return V, dx, cfg, pos, Z


def deplane(V):
    """DISPLAY ONLY: remove the in-plane phase-ramp gauge of fixed-probe recons, with the
    analyze_thin_campaign model (one plane on the depth-summed phase, spread over the layers).
    A gauge, not structure; atom coordinates are untouched."""
    s = V.sum(0); yy, xx = np.mgrid[0:s.shape[0], 0:s.shape[1]].astype(float)
    A = np.c_[xx.ravel(), yy.ravel(), np.ones(xx.size)]
    c, *_ = np.linalg.lstsq(A, s.ravel(), rcond=None)
    return V - ((c[0] * xx + c[1] * yy) / V.shape[0])[None]


def planes(pos, Z, sp):
    """Distinct depths of one species' atomic planes in the GT (1.95 A apart for AO vs BO2)."""
    z = np.sort(pos[Z == sp, 2]); out = []
    for v in z:
        if not out or v - out[-1][-1] > 0.5:
            out.append([v])
        else:
            out[-1].append(v)
    return np.array([np.mean(g) for g in out])


def pick_row(found, sp, ny, dx):
    """Row through the species' column with the most found atoms, nearest the field centre."""
    f = found[found["species"] == sp]
    if len(f) == 0:
        return None
    ids, n = np.unique(f["col_id"], return_counts=True)
    best = None
    for cid in ids[n >= n.max() - 1]:
        r = np.median(f["row"][f["col_id"] == cid])
        c = np.median(f["col"][f["col_id"] == cid])
        d = np.hypot(r - ny / 2, c - ny / 2)
        if best is None or d < best[0]:
            best = (d, r)
    return best[1]


def section(V, row, half_px):
    r0, r1 = int(round(row - half_px)), int(round(row + half_px)) + 1
    return V[:, max(r0, 0):min(r1, V.shape[1]), :].mean(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recons", required=True, help="dir holding recon_af_a<A>_lab_NL<NL>")
    ap.add_argument("--atomfind", default=None, help="dir holding atomfind_a<A>/ (found_atoms.npy, report.json)")
    ap.add_argument("--gt", required=True, help="thin ground-truth data dir (gt_prepared.npz)")
    ap.add_argument("--half-A", type=float, default=0.25, help="half-thickness of each x-z slab (A)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    runs = discover(os.path.expanduser(a.recons))
    if not runs:
        raise SystemExit(f"no recon_af_a*_lab_NL*/**/*_recons.h5 under {a.recons}")
    print("alphas:", [r[0] for r in runs])

    fig, ax = plt.subplots(len(runs), 3, figsize=(14, 6.0 * len(runs)), squeeze=False,
                           constrained_layout=True)
    for i, (A, NL, h5) in enumerate(runs):
        V, dx, cfg, pos, Z = load(h5, NL, os.path.expanduser(a.gt))
        V = deplane(V)
        dz = cfg.dz; nL, ny, nx = V.shape
        found, al = None, None
        if a.atomfind:
            d = os.path.join(os.path.expanduser(a.atomfind), f"atomfind_a{A}")
            if os.path.exists(os.path.join(d, "found_atoms.npy")):
                found = np.load(os.path.join(d, "found_atoms.npy"))
                al = json.load(open(os.path.join(d, "report.json")))["align"]
        half_px = a.half_A / dx
        band = slice(int(round(ZVAC / dz)), nL - int(round(ZVAC / dz)))

        # (1) in-plane, summed over the atomic slab
        img = V[band].sum(0)
        a0 = ax[i, 0]
        a0.imshow(img, cmap="magma", extent=[0, nx * dx, ny * dx, 0],
                  vmin=np.percentile(img, 1), vmax=np.percentile(img, 99.7))
        a0.set_title(f"α = {A} mrad   (NL {NL}, dz {dz:.2f} Å)\nin-plane, atomic slab", fontsize=10)
        a0.set_xlabel("x (Å)"); a0.set_ylabel("y (Å)")

        rows = {}
        if found is not None:
            rows = {"A-site row: Pb on AO | O on BO₂": pick_row(found, 82, ny, dx),
                    "B-site row: Ti / apical O alternate": pick_row(found, 22, ny, dx)}
        else:                                        # no atomfind run yet: fixed rows through the centre
            rows = {"x–z, field centre": ny / 2}
        for j, (lab, row) in enumerate(list(rows.items())[:2]):
            if row is None:
                continue
            a0.axhline(row * dx, color=["#00e5ff", "#7CFC00"][j], lw=0.8, ls="--")
            s = section(V, row, half_px)
            axx = ax[i, 1 + j]
            # equal aspect: depth is the point of the figure, so 1 A of z must look like 1 A of x
            axx.imshow(s, cmap="magma", aspect="equal", origin="upper",
                       extent=[0, nx * dx, (nL - 0.5) * dz, -0.5 * dz],
                       vmin=np.percentile(s, 1), vmax=np.percentile(s, 99.7))
            axx.axhspan(-0.5 * dz, ZVAC, color="0.6", alpha=0.25, lw=0)
            axx.axhspan(BOXZ - ZVAC, (nL - 0.5) * dz, color="0.6", alpha=0.25, lw=0)
            if al is not None:
                # ground-truth plane depths as ticks OUTSIDE the right edge: a reference that does not
                # sit on the data (full-width lines every 1.95 A read as texture, not information)
                x1 = nx * dx
                for sp, col in ((82, "#00b8cc"), (22, "#3a9d23")):
                    for zg in planes(pos, Z, sp):
                        axx.plot([x1, x1 + 0.9], [al["SGN"] * zg + al["OFF"]] * 2, color=col,
                                 lw=1.6, clip_on=False, solid_capstyle="butt")
            if found is not None:
                m = np.abs(found["row"] - row) <= 2 * half_px
                for sp, (nm, col, mk) in SP.items():
                    f = found[m & (found["species"] == sp)]
                    if len(f):
                        axx.scatter(f["col"] * dx, f["z_A"], s=16, marker=mk, facecolors="none",
                                    edgecolors=col, linewidths=1.0, label=f"found {nm}")
            axx.set_xlim(0, nx * dx); axx.set_ylim((nL - 0.5) * dz, -0.5 * dz)
            axx.set_title(f"x–z, {lab}", fontsize=10)
            axx.set_xlabel("x (Å)"); axx.set_ylabel("depth z (Å)")
            if i == 0 and found is not None:
                axx.legend(loc="lower right", fontsize=7, framealpha=0.6)
        for j in range(len(rows), 2):
            ax[i, 1 + j].axis("off")

    fig.suptitle("Multislice ptychography phase volumes vs aperture — Cs-corrected@30 mrad opened up "
                 "(known probe, full box)\n"
                 "ticks at the right edge: ground-truth AO (cyan) / BO₂ (green) planes via atomfind's depth map  ·  "
                 "grey: 4 Å z-vacuum  ·  markers: atomfind blind atoms\n"
                 "x–z slabs ±0.25 Å at equal aspect  ·  in-plane phase-ramp gauge removed for display",
                 fontsize=11)
    out = a.out or os.path.join(week_dir("figs"), "mep_volumes.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
