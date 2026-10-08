"""Re-run find_atoms_v3 locally on a run's figure volume (the same crop the finder saw), with its shipped kernels.
Optionally upsample depth by an integer factor (volume AND kernels, linear interpolation)."""
import sys, json, os, glob, time, copy, numpy as np
sys.path.insert(0, 'analysis'); sys.path.insert(0, 'analysis/atomfind/..')
import make_arm_page as mp, figdata, depth_honesty as dh
from atomfind import align, psf as psfmod, find, validate
from scipy.ndimage import zoom
R = mp.Runs(json.load(open(mp.MANIFEST)))
NAN_SLOTS = [0]
_orig_guided = find._fit_single_guided
def _guarded(resid, K, Kg, l_s, r_s, c_s, *a, **k):
    """harness-side guard (package untouched): a slot predicted at NaN is skipped and counted."""
    if not np.all(np.isfinite([l_s, r_s, c_s])):
        NAN_SLOTS[0] += 1
        return None
    return _orig_guided(resid, K, Kg, l_s, r_s, c_s, *a, **k)
find._fit_single_guided = _guarded

def kernels_for(info):
    an = os.path.dirname(info["af"])
    tag = os.path.basename(info["af"])[len("atomfind_"):]
    if "_kern-" in tag:                                   # analysed with another run's kernels (same-kernel analyses)
        src = tag.split("_kern-")[1]
        k = R.find("known_a" + src.split("armf_a")[1][:3]) if src.startswith("armf_a") and len(src) == 9 else None
        if k and k["af"]:
            an = os.path.dirname(k["af"]); tag = src
    pb = glob.glob(os.path.join(an, "psf", f"psf_Pb_{tag}_vol.npy")) or glob.glob(os.path.join(an, "psf", "psf_Pb_*_vol.npy"))
    ti = glob.glob(os.path.join(an, "psf", f"psf_Ti_{tag}_vol.npy")) or glob.glob(os.path.join(an, "psf", "psf_Ti_*_vol.npy"))
    return pb[0], ti[0]

def upz(A, f):
    """Linear interpolation in depth onto slices f times thinner, centred like the originals (slice i at (i+.5)dz)."""
    nl = A.shape[0]
    znew = (np.arange(nl * f) + 0.5) / f - 0.5          # in units of old slice index
    znew = np.clip(znew, 0, nl - 1)
    i0 = np.floor(znew).astype(int); i1 = np.minimum(i0 + 1, nl - 1); w = (znew - i0)[:, None, None]
    return (1 - w) * A[i0] + w * A[i1]

def upz_kernel(K, f):
    """Node-aligned: old slice i -> new slice f*i, so an odd kernel stays odd and centred (2hz+1 -> 2f*hz+1)."""
    nl = K.shape[0]
    znew = np.arange(f * (nl - 1) + 1) / f
    i0 = np.minimum(np.floor(znew).astype(int), nl - 1); i1 = np.minimum(i0 + 1, nl - 1); w = (znew - i0)[:, None, None]
    return (1 - w) * K[i0] + w * K[i1]

def run(key, f=1):
    info = R.find(key)
    V, cfg, pos, Z, found0, al0 = dh.frame(info)
    pb, ti = kernels_for(info)
    c = copy.copy(cfg)
    c.single_atom_vol, c.ti_kernel_vol = pb, ti
    if f > 1:
        tmp = os.path.join(os.environ["S"], f"k_{key}_f{f}"); os.makedirs(tmp, exist_ok=True)
        for src, nm in ((pb, "Pb"), (ti, "Ti")):
            K = np.load(src)   # complex single-atom reconstruction: interpolate its PHASE, store exp(i phase)
            np.save(os.path.join(tmp, f"{nm}.npy"), np.exp(1j * upz(np.angle(K).astype(np.float64), f)))
        c.single_atom_vol, c.ti_kernel_vol = os.path.join(tmp, "Pb.npy"), os.path.join(tmp, "Ti.npy")
        V = upz(V, f); c.dz = cfg.dz / f
    V = V.astype(np.float64)
    t = time.time()
    kern = psfmod.species_kernels(c, c.dx if hasattr(c, "dx") else al0.dx)
    found, _ = find.find_atoms_v3(V, c, al0.dx, kern)
    al = align.refine_with_atoms(align.register(V, al0.dx, pos, Z, c), found, pos, Z, c)
    return found, c, pos, Z, al, time.time() - t, found0

if __name__ == "__main__":
    for key in sys.argv[1:]:
        for f in [int(x) for x in os.environ.get("FS", "1,2").split(",")]:
            found, c, pos, Z, al, dt, found0 = run(key, f)
            rep = dh.score(found, pos, Z, al, c); st = dh.score(found, pos, Z, al, c, tol_z=0.5, species=True)
            sh = found.copy(); sh["layer"] = sh["layer"] + 1.95 / c.dz
            shs = dh.score(sh, pos, Z, al, c, tol_z=0.5, species=True)
            fm = lambda r: f"{100*r['Pb']:3.0f}/{100*r['Ti']:3.0f}/{100*r['O']:3.0f} p{r['prec']:.2f} z{r['z']:.2f}"
            print(f"{key:11s} x{f} n={len(found):4d} (shipped {len(found0)}) {dt:5.0f}s nan-slots {NAN_SLOTS[0]} guided {int(found['guided'].sum())} | loose {fm(rep)} | strict {fm(st)} | strict, shifted 1.95 {fm(shs)}", flush=True)
