#!/usr/bin/env python
"""@file blythe_cleanup.py
@brief Free space on Blythe by rules: plan first (the default, deletes nothing), delete only with --delete.

Categories (all on by default; --only picks some):
  pulled       files listed in a --pulled manifest (relpath under ROOT <TAB> bytes) with EXACTLY that size: a verified
               copy exists elsewhere (campaign/cleanup/*_on_mac.tsv, *_archived.tsv)
  trials       outer-search trial dirs recon_af_<leg>_t<i>_NL<n>, only when search_af_<leg>/search_best.json exists
               (the search finished; its grid, traces and choice were packed)
  checkpoints  the engines' Niter*.mat inside recon_*/01/, only when the same recon dir holds its final *_recons.h5
               (ptycho/run_synthetic_recon_ML.m: checkpoints are never used to resume; extract_psf reads them only when
               there is no h5)
  exitwaves    exit_wave*.mat anywhere in the repo (written by the GPU engine, read by nothing)
  raw          data_dp.hdf5 + data_position.hdf5 of the sim families named by --raw (a '*' stands for the lab / Pb / Ti
               leg); a family matching --keep is never touched, whatever --raw says
  prev         analysis.prev_* dirs inside recon dirs (an earlier run's output, moved aside by the driver)

Guards, checked when planning AND again just before each deletion: paths stay under ROOT, symlinks are never
followed or deleted through (a recon dir's links into a sim dir are removed as links); nothing a queued or running
job of ours names (squeue job names af_sim_/af_rec_/af_kick_/af_srch_<run> protect every dir containing <run>);
nothing modified within --min-age-hours. Every plan goes to logs/cleanup_plan_<ts>.tsv, every deletion to
logs/cleanup_deleted_<ts>.tsv.

    python campaign/blythe_cleanup.py @campaign/cleanup/2026-10-14_args.txt            # plan (the policy: one arg per line)
    python campaign/blythe_cleanup.py @campaign/cleanup/2026-10-14_args.txt --delete   # the same plan, then deleted
"""
import argparse
import os
import re
import shutil
import stat
import subprocess
import sys
import time

GB = 1e9
CATS = ("pulled", "trials", "checkpoints", "exitwaves", "raw", "prev")
TRIAL = re.compile(r"^recon_af_.+_t\d+_NL\d+$")       # an outer-search trial dir


def legs(pattern):
    """sim family pattern -> exact dir names ('*' = lab / Pb / Ti)."""
    return [pattern.replace("*", m) for m in ("lab", "Pb", "Ti")] if "*" in pattern else [pattern]


def tree_size(d):
    """(bytes, newest mtime) of the regular files under d, links not followed."""
    s, t = 0, 0.0
    for dp, dn, fn in os.walk(d):
        for f in fn:
            try:
                st = os.lstat(os.path.join(dp, f))
            except OSError:
                continue
            if stat.S_ISREG(st.st_mode):
                s += st.st_size; t = max(t, st.st_mtime)
    return s, t


def queued_runs():
    """Run names our queued / running jobs carry (None if squeue is unavailable)."""
    try:
        out = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%j"],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    runs = set()
    for j in out.stdout.split():
        m = re.match(r"af_(?:sim|rec|kick|srch|select)_(.+)$", j)
        if m:
            runs.add(m.group(1))
    return runs


class Planner:
    def __init__(self, a):
        self.a = a
        self.root = os.path.realpath(a.root)
        self.repo = os.path.realpath(a.repo or os.path.join(self.root, "ptyco_baseline", "ptyco_pipeline"))
        self.now = time.time()
        self.runs = queued_runs()
        self.items, self.skipped = [], []

    def inside(self, p):
        rp = os.path.realpath(os.path.dirname(p))           # the parent resolved; p itself must not be a link
        return (rp + os.sep).startswith(self.root + os.sep) and not os.path.islink(p)

    def protected(self, p, mtime):
        """A reason not to delete p, or ''."""
        if not self.inside(p):
            return "outside ROOT or a symlink"
        if self.now - mtime < self.a.min_age_hours * 3600:
            return f"modified within {self.a.min_age_hours:g} h"
        if self.runs is None:
            return "squeue unavailable: cannot check the queue"
        name = os.path.relpath(p, self.repo).split(os.sep)[0]
        for r in self.runs:
            if r in name:
                return f"a queued / running job names it ({r})"
        return ""

    def add(self, cat, p, kind, size, mtime):
        why = self.protected(p, mtime)
        (self.skipped if why else self.items).append((cat, p, kind, size, why))

    def repo_kids(self, pat):
        try:
            return sorted(e.path for e in os.scandir(self.repo) if e.is_dir(follow_symlinks=False) and re.match(pat, e.name))
        except OSError:
            return []

    # ------------------------------------------------------------------ categories
    def pulled(self):
        for m in self.a.pulled:
            for line in open(m):
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                p, want = os.path.join(self.root, f[0]), int(f[1])
                try:
                    st = os.lstat(p)
                except OSError:
                    continue                                   # already gone
                if not stat.S_ISREG(st.st_mode):
                    self.skipped.append(("pulled", p, "file", 0, "not a regular file")); continue
                if st.st_size != want:
                    self.skipped.append(("pulled", p, "file", st.st_size, f"size {st.st_size} != manifest {want}")); continue
                self.add("pulled", p, "file", st.st_size, st.st_mtime)

    def trials(self):
        for d in self.repo_kids(r"recon_af_.+_t\d+_NL\d+$"):
            leg = re.match(r"recon_af_(.+)_t\d+_NL\d+$", os.path.basename(d)).group(1)
            if not os.path.isfile(os.path.join(self.repo, f"search_af_{leg}", "search_best.json")):
                self.skipped.append(("trials", d, "dir", 0, "its search has no search_best.json")); continue
            s, t = tree_size(d)
            self.add("trials", d, "dir", s, t)

    def checkpoints(self):
        for d in self.repo_kids(r"recon_"):
            if TRIAL.search(os.path.basename(d)):
                continue                                       # trial dirs go whole (trials)
            has_h5 = False
            for dp, dn, fn in os.walk(os.path.join(d, "analysis")):
                dn[:] = [x for x in dn if not x.startswith("analysis.prev_")]
                if any(re.search(r"_recons(_\d+)?\.h5$", f) and os.path.isfile(os.path.join(dp, f))
                       and not os.path.islink(os.path.join(dp, f)) for f in fn):
                    has_h5 = True; break
            if not has_h5:
                continue
            for dp, dn, fn in os.walk(os.path.join(d, "01")):
                for f in fn:
                    if re.match(r"Niter\d+\.mat$", f):
                        p = os.path.join(dp, f); st = os.lstat(p)
                        if stat.S_ISREG(st.st_mode):
                            self.add("checkpoints", p, "file", st.st_size, st.st_mtime)

    def exitwaves(self):
        for dp, dn, fn in os.walk(self.repo):
            dn[:] = [x for x in dn if x != ".git" and not TRIAL.search(x)]   # trial dirs go whole (trials)
            for f in fn:
                if f.startswith("exit_wave") and f.endswith(".mat"):
                    p = os.path.join(dp, f); st = os.lstat(p)
                    if stat.S_ISREG(st.st_mode):
                        self.add("exitwaves", p, "file", st.st_size, st.st_mtime)

    def raw(self):
        keep = {n for k in self.a.keep for n in legs(k)}
        for pat in self.a.raw:
            for name in legs(pat):
                if name in keep:
                    self.skipped.append(("raw", os.path.join(self.repo, name), "family", 0, "in --keep")); continue
                for f in ("data_dp.hdf5", "data_position.hdf5"):
                    p = os.path.join(self.repo, name, "01", f)
                    try:
                        st = os.lstat(p)
                    except OSError:
                        continue
                    if stat.S_ISREG(st.st_mode):
                        self.add("raw", p, "file", st.st_size, st.st_mtime)

    def prev(self):
        for d in self.repo_kids(r"recon_"):
            if TRIAL.search(os.path.basename(d)):
                continue                                       # trial dirs go whole (trials)
            try:
                kids = [e.path for e in os.scandir(d) if e.is_dir(follow_symlinks=False) and e.name.startswith("analysis.prev_")]
            except OSError:
                continue
            for k in kids:
                s, t = tree_size(k)
                self.add("prev", k, "dir", s, t)


def remove(item, P):
    """Delete one planned item after re-checking it (returns '' on success, else the reason it was left)."""
    cat, p, kind, size, _ = item
    try:
        st = os.lstat(p)
    except OSError:
        return "gone already"
    if kind == "dir":
        if not stat.S_ISDIR(st.st_mode):
            return "no longer a real dir"
        why = P.protected(p, tree_size(p)[1])
        if why:
            return why
        shutil.rmtree(p)                                       # links inside are removed as links, never followed
    else:
        if not stat.S_ISREG(st.st_mode) or st.st_size != size:
            return "changed since planned"
        why = P.protected(p, st.st_mtime)
        if why:
            return why
        os.remove(p)
    return ""


def main():
    # @<file>: arguments one per line (--raw=<family>), so a family's '*' never meets a shell
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 fromfile_prefix_chars="@")
    ap.add_argument("--root", default="/springbrook/share/physics/phucrh")
    ap.add_argument("--repo", default=None)
    ap.add_argument("--pulled", action="append", default=[], help="manifest of verified copies (repeatable)")
    ap.add_argument("--raw", action="append", default=[], help="sim family whose raw 4-D data goes (repeatable)")
    ap.add_argument("--keep", action="append", default=[], help="sim family never touched (repeatable)")
    ap.add_argument("--only", nargs="+", choices=CATS, default=list(CATS))
    ap.add_argument("--min-age-hours", type=float, default=12)
    ap.add_argument("--delete", action="store_true", help="delete the plan (default: plan only)")
    a = ap.parse_args()
    P = Planner(a)
    if not os.path.isdir(P.repo):
        sys.exit(f"no repo at {P.repo}")
    print(f"{'DELETING' if a.delete else 'PLAN ONLY (nothing is deleted)'}: root {P.root}  ({time.strftime('%Y-%m-%d %H:%M')})")
    print(f"queue: {'squeue unavailable -> everything is protected' if P.runs is None else (', '.join(sorted(P.runs)) or 'no af_ jobs')}")
    for c in a.only:
        getattr(P, c)()
    tot = 0
    print(f"\n{'category':12s} {'items':>6s} {'GB':>8s}   examples")
    for c in a.only:
        it = [x for x in P.items if x[0] == c]; s = sum(x[3] for x in it); tot += s
        ex = ", ".join(os.path.relpath(x[1], P.root) for x in it[:2])
        print(f"{c:12s} {len(it):6d} {s / GB:8.1f}   {ex[:150]}")
    print(f"{'TOTAL':12s} {len(P.items):6d} {tot / GB:8.1f}")
    if P.skipped:
        print(f"\nleft alone: {len(P.skipped)}")
        reasons = {}
        for x in P.skipped:
            reasons.setdefault((x[0], x[4].split(' (')[0]), []).append(x)
        for (c, why), xs in sorted(reasons.items()):
            print(f"  {c:12s} {len(xs):5d}  {why}   e.g. {os.path.relpath(xs[0][1], P.root)[:110]}")
    logs = os.path.join(P.repo, "logs"); os.makedirs(logs, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    plan = os.path.join(logs, f"cleanup_plan_{ts}.tsv")
    with open(plan, "w") as fh:
        fh.write("#category\tpath\tkind\tbytes\tleft alone because\n")
        for x in P.items + P.skipped:
            fh.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\n")
    print(f"\nplan: {plan}")
    if not a.delete:
        print("nothing was deleted (add --delete to act on this plan)")
        return 0
    done = os.path.join(logs, f"cleanup_deleted_{ts}.tsv"); freed = 0; left = 0
    with open(done, "w") as fh:
        for x in P.items:
            why = remove(x, P)
            if why:
                left += 1; fh.write(f"LEFT\t{x[1]}\t{x[3]}\t{why}\n")
            else:
                freed += x[3]; fh.write(f"DELETED\t{x[1]}\t{x[3]}\t\n")
    print(f"deleted {len(P.items) - left} items, {freed / GB:.1f} GB freed; {left} left (re-check failed); log {done}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
