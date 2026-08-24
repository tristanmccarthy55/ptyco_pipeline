"""@file __init__.py
@brief Atom finding with calibrated uncertainty for 3-D electron-ptychographic volumes.

Entry points: `run_atomfind.py` (localisation + report), `polarisation.py` (Ti-O6 off-centring).
See README.md for usage, METHODS.md for the method, RESULTS.md for measured numbers, and
PEER.md for the reproduction protocol. `test_atomfind.py` validates an install in ~1 s
without needing the data.
"""
__all__ = ["config", "align", "psf", "deconv", "find", "fit", "validate", "uncertainty"]
