# Blythe clean-ups — what was removed, and where the copies are

One set of files per clean-up, named by date. `campaign/blythe_disk_survey.py` measures first (read-only);
`campaign/blythe_cleanup.py @<date>_args.txt` plans by rules, and deletes only with `--delete`.

## 2026-10-14 (our share: 1.52 TB in 243k files; the filesystem itself 12 % full, no per-user quota)

| file | what it is |
|---|---|
| `2026-10-14_args.txt` | the policy: manifests of verified copies, sim families whose raw 4-D data goes, families kept |
| `2026-10-14_on_mac.tsv` | 116 Blythe files (37 GB of tarballs) with a byte-identical copy on the Mac, and where it is |
| `2026-10-14_archive.txt` | 44 files (19 GB) with NO copy elsewhere, packed into one tar and pulled before anything goes: 28 tarballs never pulled, and the 16 reconstruction h5s whose probe was RECOVERED (kick, refined, true start, long kicks: the probe exists nowhere else) |
| `2026-10-14_archived.tsv` | written once the pulled archive is verified on the Mac; the archived files then join the clean-up |

Kept on Blythe: the raw sims of the four phase-1 points that failed (focal spread 5 and 7.29 Å, source 0.2 and 0.4 Å),
which the probe-mode legs reuse; every reconstruction h5 with a known or fixed probe (~84 GB; their volumes, kernels and
atoms are in the pulled analysis tarballs); checkpoints in recon dirs that have no h5; the analysis dirs (kernel
fallbacks); envs, dft-relax, eels, fold_slice.

Removed: tarballs with a Mac copy; outer-search trial dirs (both searches' grids, traces and choices are packed);
solver checkpoints beside a final h5; every exit_wave*.mat (read by nothing); raw 4-D data of the six phase-1 points
that passed and of the ARM full-physics sims (regenerable with the same draws: `campaign/run_thin_atomfind.sh` with the
row's launch line, ~14 GPU-h per three-leg family); analysis.prev_* dirs.
