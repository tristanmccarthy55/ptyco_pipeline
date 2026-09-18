#!/usr/bin/env python
"""@file pubstyle.py
@brief Shared style for the publication figures in analysis/paper/.

One claim per figure, one style for all of them. The exploratory scripts one level up
(`make_ronchigram_fig.py`, `make_mep_volumes_fig.py`, `collate_atomfind_depth.py`) stay as
DIAGNOSTICS: they carry every check needed while working. These are the manuscript versions.

Conventions here:
  - journal column widths (single 85 mm, 1.5 = 127 mm, double 180 mm), heights chosen per figure;
  - 8 pt text, 7 pt ticks/legends -- readable at final size, no rescaling in the manuscript;
  - Okabe-Ito colours (colour-blind safe), with ONE species colour used across every figure;
  - PDF (vector, TrueType 42 so the text stays editable) + a 600 dpi PNG for slides;
  - no gridlines, no top/right spines, ticks inward.
"""
from __future__ import annotations
import datetime, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- widths (inches) --------------------------------------------------------------
COL1, COL15, COL2 = 3.35, 5.00, 7.09        # 85 / 127 / 180 mm

# ---- colour (Okabe-Ito) -----------------------------------------------------------
BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"
SKY, YELLOW, PURPLE, GREY = "#56B4E9", "#E69F00", "#CC79A7", "#5A5A5A"
SPECIES = {82: ("Pb", BLUE), 22: ("Ti", ORANGE), 8: ("O", GREEN)}
SPECIES_NAME = {"Pb": BLUE, "Ti": ORANGE, "O": GREEN}


def use():
    """Apply the publication rcParams. Call once at the top of each figure script."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "lines.linewidth": 1.2, "lines.markersize": 3.5,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.8, "ytick.major.size": 2.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "legend.frameon": False,
        "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01, "pdf.fonttype": 42, "ps.fonttype": 42,
        "image.interpolation": "nearest",
    })


def panel(ax, letter, dx=-0.18, dy=1.02, **kw):
    """Bold lower-case panel label at the axes' top-left, outside the frame."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", **kw)


def scalebar(ax, length, label, color="white", y=0.06, x0=0.06, lw=2.0):
    """Horizontal scale bar in AXES fraction, annotating a length in data units."""
    x1 = x0 + length / (ax.get_xlim()[1] - ax.get_xlim()[0])
    ax.plot([x0, x1], [y, y], transform=ax.transAxes, color=color, lw=lw,
            solid_capstyle="butt", clip_on=False)
    ax.text((x0 + x1) / 2, y + 0.025, label, transform=ax.transAxes, color=color,
            ha="center", va="bottom", fontsize=7)


def imshow_clean(ax, img, **kw):
    """Image panel with no ticks (use a scale bar instead)."""
    h = ax.imshow(img, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.5)
    return h


def week_dir(sub="figs"):
    """<repo>/aberration_experiment/<sub>/<ISO week>/paper/ -- created on demand."""
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    d = os.path.join(repo, "aberration_experiment", sub,
                     datetime.date.today().strftime("%G-W%V"), "paper")
    os.makedirs(d, exist_ok=True)
    return d


def save(fig, name, out=None):
    """Write <name>.pdf (vector, for the manuscript) and <name>.png (600 dpi, for slides)."""
    d = out or week_dir()
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(d, f"{name}.{ext}")
        fig.savefig(p)
        paths.append(p)
    print("wrote " + "  ".join(paths))
    return paths
