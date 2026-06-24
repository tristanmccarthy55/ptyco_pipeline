"""Montage the 8 orientation-sweep reconstructions for the phantom test.

Run on Blythe after run_phantom_test.slurm finishes, then scp the single PNG down:
    $SHARE/phucrh/envs/abtem/bin/python sim/montage_orientations.py
    # writes sim_phantom/orientation_montage.png

The correct detector orientation (custom_data_flip) is the panel whose 'F' is
upright: vertical bar on the LEFT, two arms pointing +x (right), the TOP arm higher
than the MIDDLE arm — matching the phantom (see the 'EXPECTED' panel).
"""
from __future__ import annotations
import sys, glob, os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = sys.argv[1] if len(sys.argv) > 1 else "sim_phantom/01"
FLIPS = ["000", "100", "010", "001", "110", "101", "011", "111"]
LABELS = {"000": "identity", "100": "flipLR", "010": "flipUD", "001": "transpose",
          "110": "rot180", "101": "rot90", "011": "rot270", "111": "anti-transpose"}


def load_orient(code):
    pats = sorted(glob.glob(os.path.join(BASE, f"*orient_{code}*", "O_phase_roi",
                                         "O_phase_roi_Niter*_Layer1.tiff")))
    if not pats:
        return None
    return np.array(Image.open(pats[-1])).astype(float)   # last (highest Niter)


# expected F (phantom coords, x->right y->up) for visual reference
def expected_F():
    coords = [(38, 14), (38, 17), (38, 20), (38, 23), (38, 26),
              (41, 26), (44, 26), (41, 20)]
    img = np.zeros((40, 40))
    for x, y in coords:
        img[int((y - 10) * 2), int((x - 30) * 2)] = 1
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(img, 1.2)


fig, axes = plt.subplots(2, 5, figsize=(20, 8))
axes = axes.ravel()
axes[0].imshow(expected_F(), origin="lower", cmap="inferno")
axes[0].set_title("EXPECTED F\n(upright = correct)", fontsize=10)
axes[0].axis("off")
axes[5].axis("off")
slots = [1, 2, 3, 4, 6, 7, 8, 9]
for code, slot in zip(FLIPS, slots):
    ax = axes[slot]
    o = load_orient(code)
    if o is None:
        ax.text(0.5, 0.5, f"{code}\n(missing)", ha="center", va="center"); ax.axis("off"); continue
    ax.imshow(o, origin="lower", cmap="gray")
    ax.set_title(f"flip {code}  ({LABELS[code]})", fontsize=10)
    ax.axis("off")
fig.suptitle("Orientation sweep — which custom_data_flip renders the F upright?", fontsize=13)
fig.tight_layout()
out = os.path.join(os.path.dirname(BASE) or ".", "orientation_montage.png")
fig.savefig(out, dpi=120)
print("wrote", out)
