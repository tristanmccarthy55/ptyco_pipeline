#!/usr/bin/env python
"""@file build_page.py
@brief Build the published logbook page: figures in, one self-contained HTML out.

The page template carries {{FIG_<name>}} placeholders; this fills each with a base64 data URI of
the corresponding PNG, downscaled to web width. Figures must be EMBEDDED, not published as separate
artifact files referenced by relative path -- that does not render (see HANDOFF_ANALYSIS.md, Traps).

Keeping the template and this script in the repo means the next agent can rebuild the page from the
current figures and republish it to the SAME artifact URL, instead of starting a second page.

    ~/hyperspy-bundle/bin/python aberration_experiment/page/build_page.py --out <dir>/logbook_built.html
    ~/hyperspy-bundle/bin/python aberration_experiment/page/build_page.py \\
        --template aberration_experiment/page/round_2026-09-22.html --out <dir>/round_built.html
"""
from __future__ import annotations
import argparse, base64, io, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(os.path.dirname(HERE), "figs")

# Per template: placeholder -> figure, in the order they appear on that page. Two pages are built
# from this one script -- the running logbook and each round's results page -- so a figure is named
# once and both stay in step with the run outputs.
SOURCES_BY_TEMPLATE = {
    "logbook.html": {
        "FIG_RONCHI":   "2026-W39/ronchigram_evolution.png",
        "FIG_DEPTH":    "2026-W37/atomfind_depth_vs_alpha.png",
        "FIG_VOLUMES":  "2026-W37/mep_volumes.png",
        "FIG_C1":       "2026-W39/meeting/fig4_c1.png",
        "FIG_STEP1":    "2026-W39/meeting/figB_step1.png",
        "FIG_FOCUS":    "2026-W39/simple/fig_focus.png",
        "FIG_DOSE":     "2026-W39/simple/fig_dose.png",
        "FIG_PROBE":    "2026-W39/simple/fig_probe.png",
        "FIG_PSF":      "2026-W39/simple/fig_psf.png",
        "FIG_GROWTH":   "2026-W39/simple/fig_growth.png",
        "FIG_P1":       "2026-W38/paper/fig1_probe.png",
        "FIG_P2":       "2026-W38/paper/fig2_depth.png",
        "FIG_P3":       "2026-W38/paper/fig3_volumes.png",
        "FIG_P4":       "2026-W38/paper/fig4_baselines.png",
    },
    "round_2026-09-22.html": {
        "FIG_PROBE":    "2026-W39/simple/fig_probe.png",
        "FIG_FOCUS":    "2026-W39/simple/fig_focus.png",
        "FIG_DOSE":     "2026-W39/simple/fig_dose.png",
        "FIG_NONROUND": "2026-W39/simple/fig_nonround.png",
        "FIG_BLAME":    "2026-W39/simple/fig_didfinderfail_0p2.png",
        "FIG_PSF":      "2026-W39/simple/fig_psf.png",
    },
}

MAXW = 1700          # web width; the PDFs in the repo stay the print copies
PNG_CAP = 420_000    # bytes: above this a figure goes out as JPEG instead


def encode(path: str, maxw: int = MAXW) -> tuple[str, int]:
    """PNG data URI, or JPEG when the PNG is heavy (the page must stay well under 16 MB)."""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    if im.width > maxw:
        im = im.resize((maxw, round(im.height * maxw / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, "PNG", optimize=True); raw, mime = buf.getvalue(), "png"
    if len(raw) > PNG_CAP:
        buf = io.BytesIO(); im.save(buf, "JPEG", quality=88, optimize=True, progressive=True)
        raw, mime = buf.getvalue(), "jpeg"
    return f"data:image/{mime};base64," + base64.b64encode(raw).decode(), len(raw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", default=os.path.join(HERE, "logbook.html"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    html = open(a.template).read()
    name = os.path.basename(a.template)
    if name not in SOURCES_BY_TEMPLATE:
        sys.exit(f"no figure list for template {name}; add one to SOURCES_BY_TEMPLATE")
    total = 0
    for key, rel in SOURCES_BY_TEMPLATE[name].items():
        src = os.path.join(FIGS, rel)
        if not os.path.exists(src):
            sys.exit(f"missing figure: {src}")
        uri, n = encode(src)
        if "{{%s}}" % key not in html:
            sys.exit(f"template has no placeholder for {key}")
        html = html.replace("{{%s}}" % key, uri)
        total += n
        print(f"  {key:12s} {rel:38s} {n/1024:7.0f} kB")
    left = [t for t in html.split("{{")[1:]]
    if left:
        sys.exit("unfilled placeholders: " + ", ".join(t.split("}}")[0] for t in left))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w").write(html)
    print(f"wrote {a.out}  ({os.path.getsize(a.out)/1e6:.2f} MB, figures {total/1e6:.2f} MB raw)")


if __name__ == "__main__":
    main()
