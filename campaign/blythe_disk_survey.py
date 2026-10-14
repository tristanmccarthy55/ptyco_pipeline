#!/usr/bin/env python
"""@file blythe_disk_survey.py
@brief READ-ONLY survey of our space on Blythe: what takes the room, by kind and by run family, before a clean-up.

Walks the user's share dir (never following symlinks: recon dirs link their data into the sim dirs, so a followed link
would count the 4-D data twice) and prints a short report; the full detail goes to a TSV for analysis off Blythe. It
deletes nothing and writes only the TSV.

Report sections: the share's free space; every top-level entry of ROOT; totals by KIND (raw 4-D data, reconstruction
h5s, solver checkpoints, tarballs, ...); the repo's run dirs grouped into FAMILIES (lab / Pb / Ti legs and search
trials folded together: sim_out_af_armf_a080_lab_fs5 -> sim_out_af_armf_a080_*_fs5); anything of ours in the group root.
TSV rows: D (dirs >= --min-mb, or depth <= 2), F (files >= --big-mb, and EVERY tarball, with its byte size so it can be
matched against the copies already pulled), FAM (families), K (kinds).

    /springbrook/share/physics/phucrh/envs/abtem/bin/python campaign/blythe_disk_survey.py
    python campaign/blythe_disk_survey.py --root <dir> --repo <dir> --out-dir <dir>      # anywhere else (tests)
"""
import argparse
import os
import re
import shutil
import stat
import sys
import time

GB = 1e9


def kind_of(path, name):
    """What a regular file is, for the totals (path-based classes first: they override the file name)."""
    if "/analysis.prev_" in path:
        return "old analysis (analysis.prev_*)"
    if "/envs/" in path or "/conda" in path or "/.conda/" in path:
        return "conda envs"
    if "/.vscode-server/" in path or "/.cache/" in path:
        return "editor / caches"
    if "/.git/" in path:
        return "git objects"
    if name in ("data_dp.hdf5", "data_position.hdf5"):
        return "raw 4-D data (sims)"
    if name.endswith("_recons.h5"):
        return "reconstruction h5"
    if name.endswith((".h5", ".hdf5")):
        return "other h5"
    if re.match(r"Niter\d+\.mat$", name):
        return "solver checkpoints (Niter*.mat)"
    if name.endswith(".mat"):
        return "other .mat"
    if name.endswith((".tgz", ".tar.gz", ".tar")):
        return "tarballs"
    if name.endswith((".npy", ".npz")):
        return "numpy arrays"
    if name.endswith((".out", ".err", ".log")):
        return "logs"
    return "other"


def family(name):
    """A run dir's family: legs (lab / Pb / Ti), search trials (_t<i>) and the NL tag folded together."""
    n = re.sub(r"_NL\d+$", "", name)
    n = re.sub(r"_t\d+(?=_|$)", "_t*", n)
    n = re.sub(r"_(lab|Pb|Ti)(?=_|$)", "_*", n, count=1)
    return n


class Tot:
    __slots__ = ("size", "alloc", "n", "mtime")

    def __init__(self):
        self.size = 0; self.alloc = 0; self.n = 0; self.mtime = 0.0

    def add(self, size, alloc, n, mtime):
        self.size += size; self.alloc += alloc; self.n += n
        if mtime > self.mtime:
            self.mtime = mtime


def walk(root, big_bytes, progress=50000):
    """One pass over root, no symlinks followed. Returns per-dir recursive totals, per-kind totals, big files,
    tarballs, and the number of unreadable dirs."""
    dirs, kinds, big, tars = {}, {}, [], []
    denied = 0; seen = 0; t0 = time.time()
    order = []                                   # dirs in visit order, to roll totals up children-first
    stack = [root]
    while stack:
        d = stack.pop()
        dirs[d] = Tot(); order.append(d)
        try:
            it = os.scandir(d)
        except OSError:
            denied += 1
            continue
        with it:
            for e in it:
                try:
                    st = e.stat(follow_symlinks=False)
                except OSError:
                    continue
                if stat.S_ISDIR(st.st_mode):
                    stack.append(e.path)
                    continue
                if not stat.S_ISREG(st.st_mode):  # symlinks, sockets: a link is a few bytes, its target counted once
                    continue
                alloc = getattr(st, "st_blocks", 0) * 512
                dirs[d].add(st.st_size, alloc, 1, st.st_mtime)
                k = kind_of(e.path, e.name)
                kinds.setdefault(k, Tot()).add(st.st_size, alloc, 1, st.st_mtime)
                if st.st_size >= big_bytes:
                    big.append((e.path, st.st_size, st.st_mtime))
                if k == "tarballs":
                    tars.append((e.path, st.st_size, st.st_mtime))
                seen += 1
                if seen % progress == 0:
                    print(f"  ... {seen} files, {time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    for d in reversed(order):                    # children were visited after their parent: roll up in reverse
        p = os.path.dirname(d)
        if d != root and p in dirs:
            t = dirs[d]; dirs[p].add(t.size, t.alloc, t.n, t.mtime)
    return dirs, kinds, big, tars, denied, seen


def day(t):
    return time.strftime("%Y-%m-%d", time.localtime(t)) if t else "-"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="/springbrook/share/physics/phucrh")
    ap.add_argument("--repo", default=None, help="default <root>/ptyco_baseline/ptyco_pipeline")
    ap.add_argument("--group-root", default=None, help="checked for entries we own (default: the parent of --root)")
    ap.add_argument("--out-dir", default=None, help="where the TSV goes (default <repo>/logs)")
    ap.add_argument("--min-mb", type=float, default=50, help="TSV: dirs at least this big (plus every dir of depth <= 2)")
    ap.add_argument("--big-mb", type=float, default=100, help="TSV: files at least this big (plus every tarball)")
    ap.add_argument("--top", type=int, default=45, help="families printed")
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    repo = os.path.abspath(a.repo or os.path.join(root, "ptyco_baseline", "ptyco_pipeline"))
    out_dir = os.path.abspath(a.out_dir or os.path.join(repo, "logs"))
    group = os.path.abspath(a.group_root or os.path.dirname(root))
    t0 = time.time()

    du = shutil.disk_usage(root)
    print(f"survey of {root}  ({time.strftime('%Y-%m-%d %H:%M')})")
    print(f"filesystem: {du.total / 1e12:.2f} TB, used {du.used / 1e12:.2f} TB, FREE {du.free / 1e12:.2f} TB "
          f"({100 * du.free / du.total:.0f} %)")
    dirs, kinds, big, tars, denied, seen = walk(root, a.big_mb * 1e6)
    T = dirs[root]
    print(f"ours: {T.size / GB:.1f} GB in {T.n} files (allocated {T.alloc / GB:.1f} GB); walked in {time.time() - t0:.0f} s"
          f"{f'; {denied} dirs unreadable' if denied else ''}")

    print("\n== top level of ROOT (GB, files, newest)")
    kids = sorted((p for p in dirs if os.path.dirname(p) == root), key=lambda p: -dirs[p].size)
    loose = Tot()
    for p in os.scandir(root):
        if p.is_file(follow_symlinks=False):
            st = p.stat(follow_symlinks=False); loose.add(st.st_size, 0, 1, st.st_mtime)
    for p in kids[:30]:
        t = dirs[p]
        print(f"  {t.size / GB:8.1f}  {t.n:8d}  {day(t.mtime)}  {os.path.basename(p)}/")
    if len(kids) > 30:
        r = Tot()
        for p in kids[30:]:
            t = dirs[p]; r.add(t.size, 0, t.n, t.mtime)
        print(f"  {r.size / GB:8.1f}  {r.n:8d}  {day(r.mtime)}  ({len(kids) - 30} more dirs)")
    print(f"  {loose.size / GB:8.1f}  {loose.n:8d}  {day(loose.mtime)}  (files directly in ROOT)")

    print("\n== by kind, everything under ROOT (GB, files, newest)")
    for k, t in sorted(kinds.items(), key=lambda kv: -kv[1].size):
        print(f"  {t.size / GB:8.1f}  {t.n:8d}  {day(t.mtime)}  {k}")

    fams = {}
    if repo in dirs:
        R = dirs[repo]
        print(f"\n== repo {repo}: {R.size / GB:.1f} GB. Run dirs by family "
              f"(GB total / raw 4-D / recon h5 / checkpoints, dirs, newest)")
        kidsr = [p for p in dirs if os.path.dirname(p) == repo]
        # kind split per child needs a second, shallow pass over each run dir (cheap: few files each)
        for p in kidsr:
            name = os.path.basename(p)
            if not re.match(r"(sim_out|recon|search)_", name):
                continue
            raw = h5 = ck = 0
            for dp, dn, fn in os.walk(p):
                for f in fn:
                    fp = os.path.join(dp, f)
                    try:
                        st = os.lstat(fp)
                    except OSError:
                        continue
                    if not stat.S_ISREG(st.st_mode):
                        continue
                    k = kind_of(fp, f)
                    if k == "raw 4-D data (sims)":
                        raw += st.st_size
                    elif k == "reconstruction h5":
                        h5 += st.st_size
                    elif k == "solver checkpoints (Niter*.mat)":
                        ck += st.st_size
            f = fams.setdefault(family(name), dict(t=Tot(), raw=0, h5=0, ck=0, n=0))
            t = dirs[p]; f["t"].add(t.size, t.alloc, t.n, t.mtime); f["raw"] += raw; f["h5"] += h5; f["ck"] += ck; f["n"] += 1
        order = sorted(fams.items(), key=lambda kv: -kv[1]["t"].size)
        for name, f in order[:a.top]:
            print(f"  {f['t'].size / GB:7.1f} {f['raw'] / GB:7.1f} {f['h5'] / GB:7.1f} {f['ck'] / GB:6.1f}  {f['n']:4d}  "
                  f"{day(f['t'].mtime)}  {name}")
        if len(order) > a.top:
            r = Tot(); rr = rh = rc = rn = 0
            for _, f in order[a.top:]:
                r.add(f["t"].size, 0, 0, f["t"].mtime); rr += f["raw"]; rh += f["h5"]; rc += f["ck"]; rn += f["n"]
            print(f"  {r.size / GB:7.1f} {rr / GB:7.1f} {rh / GB:7.1f} {rc / GB:6.1f}  {rn:4d}  {day(r.mtime)}  "
                  f"({len(order) - a.top} smaller families)")
        other = [p for p in kidsr if not re.match(r"(sim_out|recon|search)_", os.path.basename(p))]
        o = Tot()
        for p in other:
            t = dirs[p]; o.add(t.size, 0, t.n, t.mtime)
        print(f"  {o.size / GB:7.1f}  (the repo's other dirs: code, logs, .git, ...)")

    print(f"\n== tarballs: {len(tars)}, {sum(s for _, s, _ in tars) / GB:.1f} GB (each listed in the TSV)")
    mine = []
    try:
        uid = os.getuid()
        for e in os.scandir(group):
            if os.path.abspath(e.path) == root:
                continue
            try:
                st = e.stat(follow_symlinks=False)
            except OSError:
                continue
            if st.st_uid == uid:
                mine.append(e.name)
    except OSError as err:
        print(f"\n== group root {group}: not listable ({err})")
    else:
        print(f"\n== group root {group}: entries we own outside ROOT: {', '.join(mine) if mine else 'none'}")

    os.makedirs(out_dir, exist_ok=True)
    tsv = os.path.join(out_dir, f"disk_survey_{time.strftime('%Y%m%d_%H%M%S')}.tsv")
    with open(tsv, "w") as fh:
        fh.write("#type\tpath\tbytes\talloc_bytes\tfiles\tnewest\n")
        for p, t in dirs.items():
            depth = 0 if p == root else os.path.relpath(p, root).count(os.sep) + 1
            if depth <= 2 or t.size >= a.min_mb * 1e6:
                fh.write(f"D\t{p}\t{t.size}\t{t.alloc}\t{t.n}\t{day(t.mtime)}\n")
        listed = set()
        for path, size, mt in big + tars:
            if path not in listed:
                listed.add(path); fh.write(f"F\t{path}\t{size}\t\t1\t{day(mt)}\n")
        for name, f in fams.items():
            fh.write(f"FAM\t{name}\t{f['t'].size}\t{f['raw']}|{f['h5']}|{f['ck']}\t{f['n']}\t{day(f['t'].mtime)}\n")
        for k, t in kinds.items():
            fh.write(f"K\t{k}\t{t.size}\t{t.alloc}\t{t.n}\t{day(t.mtime)}\n")
    print(f"\nTSV (full detail): {tsv}")
    print(f"done in {time.time() - t0:.0f} s; nothing was changed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
