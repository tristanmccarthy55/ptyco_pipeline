import os, sys, json, numpy as np
sys.argv_ = sys.argv; sys.argv = ['x']
exec(open(os.path.join(os.environ['S'], 'refind.py')).read().split('if __name__')[0])
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import make_ceos_figs as cf
SP = {82: "#2a78d6", 22: "#eb6834", 8: "#1baf7a"}
keys = sys.argv_[1:]
fig, axs = plt.subplots(len(keys), 4, figsize=(13, 4.3 * len(keys)), gridspec_kw=dict(wspace=0.12, hspace=0.32))
for i, key in enumerate(keys):
    info = R.find(key); lg = mp.load_leg(info)
    ev = json.load(open(os.path.join(os.environ['S'], f"colshift_{key}.json")))
    en = np.array([v["e_null"] for v in ev.values()], float); thr = np.nanpercentile(en, 95)
    sup = {int(k): v["e"] > thr for k, v in ev.items()}
    cuts = mp.cuts_from(mp.load_leg(R.find("known_a080")))
    V, dx, dz = lg["V"], lg["dx"], lg["dz"]; nL, ny, nx = V.shape; cfg = lg["cfg"]
    zc = (np.arange(nL) + 0.5) * dz; inner = (zc > cfg.trim_z_A[0]) & (zc < cfg.trim_z_A[1])
    for j, (cut, rA) in enumerate(cuts):
        row = rA / dx; hp = 0.25 / dx
        sec = V[:, int(round(row - hp)):int(round(row + hp)) + 1, :].mean(1)
        for k, (scale, lab) in enumerate((("all", "as published: colour scale over all layers"), ("slab", "colour scale from the slab only; rings by support"))):
            ax = axs[i, 2 * j + k]
            ref = sec if scale == "all" else sec[inner]
            ax.imshow(sec, cmap="magma", aspect="equal", origin="upper", extent=[0, nx * dx, nL * dz, 0],
                      vmin=np.percentile(ref, 1), vmax=np.percentile(ref, 99.7), interpolation="nearest")
            gr, gc, gl = lg["gt"]; m = np.abs(gr - row) <= hp
            for z_, col, sp in zip(gl[m], gc[m], lg["Z"][m]):
                if sp in SP: ax.plot(col * dx, (z_ + 0.5) * dz, ".", color=SP[sp], ms=3.4, zorder=3)
            f = lg["found"]; mm = np.abs(f["row"] - row) <= (hp if scale == "slab" else 2.2 * hp)
            for a in f[mm]:
                cc = SP.get(int(a["species"]), "white")
                ok = sup.get(int(a["col_id"]), False) or scale == "all"
                ax.plot(a["col"] * dx, (a["layer"] + 0.5) * dz, "o", mfc="none", mec=cc if ok else "white",
                        mew=0.9 if ok else 0.6, ms=5.2, zorder=4, ls="none")
            for zb in cfg.trim_z_A: ax.axhline(zb, color="#7fd4ff", lw=0.8, ls=(0, (4, 3)))
            ax.set_xlim(nx * dx / 2 - 8, nx * dx / 2 + 8); ax.set_ylim(nL * dz, 0)
            ax.set_title(f"{key} · {cut}\n{lab}", fontsize=8.5, loc="left")
            if 2 * j + k: ax.set_yticklabels([])
h = [Line2D([], [], ls="none", marker="o", mfc="none", mec="#555", ms=7, label="found, column depths supported (above its run's no-depth null)"),
     Line2D([], [], ls="none", marker="o", mfc="none", mec="white", markeredgewidth=0.6, ms=7, label="found, column depths NOT supported (white ring)")]
fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=9)
out = os.path.join(os.environ['S'], "ring_support.png"); fig.savefig(out, dpi=110, bbox_inches="tight"); print(out)
