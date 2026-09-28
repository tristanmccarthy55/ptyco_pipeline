#!/usr/bin/env python
"""@file figdata.py
@brief A reconstruction in the atom finder's own frame, for figures that draw atoms on the phase.

run_atomfind.py crops the phase to the scan field (fov_A about scan_center_xy) before it finds anything, and its
found_atoms.npy is indexed in that crop. A figure that draws the real atoms (dots) and the found ones (rings) on a
depth section must use the SAME crop, so this builds it with atomfind's own loader and config -- nothing re-derived.

  frame(h5, gt_dir, dz, dx, centre)   -> (V float32 (NL, ny, nx), dx, cfg), deplaned as the published figures are
  save(h5, gt_dir, dz, dx, centre, out)   writes out/phase_vol.npy (float16) + phase_vol.json; ~5 MB a leg
  load(out)                            -> (V, meta) from a saved copy

analyse_sweep.py calls save() for every leg it finds atoms on, so an analysis tarball carries what the figures need
while the ~1 GB reconstructions stay on Blythe.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def deplane(V):
    """Remove the in-plane phase ramp: one plane fitted to the depth sum, spread over the layers (as extract_psf)."""
    s = V.sum(0)
    yy, xx = np.mgrid[0:s.shape[0], 0:s.shape[1]].astype(float)
    dev = np.abs(s - np.median(s))
    m = dev < np.percentile(dev, 97)
    c, *_ = np.linalg.lstsq(np.c_[xx[m], yy[m], np.ones(m.sum())], s[m], rcond=None)
    return V - ((c[0] * xx + c[1] * yy) / V.shape[0])[None]


def config_for(h5, gt_dir, dz, dx, centre):
    from atomfind import config
    config.set_data_dir(os.path.expanduser(gt_dir))
    cfg = config.preset("thin")
    cfg.recon_vol, cfg.dz, cfg.dx = h5, float(dz), float(dx)
    if centre is not None:
        cfg.scan_center_xy = (float(centre[0]), float(centre[1]))
    return cfg.resolve()


def frame(h5, gt_dir, dz, dx, centre):
    from atomfind import align
    cfg = config_for(h5, gt_dir, dz, dx, centre)
    V, dxv = align.load_phase(cfg)
    V = align.crop_to_fov(V, dxv, cfg)
    V = V - np.median(V, axis=(1, 2), keepdims=True)
    return deplane(V).astype(np.float32), float(dxv), cfg


def save(h5, gt_dir, dz, dx, centre, out):
    V, dxv, _ = frame(h5, gt_dir, dz, dx, centre)
    os.makedirs(out, exist_ok=True)
    np.save(os.path.join(out, "phase_vol.npy"), V.astype(np.float16))
    json.dump(dict(dx=dxv, dz=float(dz), centre=None if centre is None else list(centre), nl=int(V.shape[0]),
                   h5=os.path.basename(h5)), open(os.path.join(out, "phase_vol.json"), "w"))
    return os.path.join(out, "phase_vol.npy")


def load(out):
    V = np.load(os.path.join(out, "phase_vol.npy")).astype(np.float32)
    return V, json.load(open(os.path.join(out, "phase_vol.json")))
