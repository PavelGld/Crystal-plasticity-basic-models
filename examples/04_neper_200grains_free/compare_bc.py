#!/usr/bin/env python
"""compare_bc.py -- compare lateral contraction profile for the two BC modes.

Prerequisites (from this folder):
    abaqus python ../../scripts/extract_u.py model       u_free.txt
    abaqus python ../../scripts/extract_u.py model_clamp u_clamp.txt

Usage:  python compare_bc.py
Writes: BC_comparison.png
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = [("free  (4 lateral faces free)", "u_free.txt", "C0"),
         ("clamp (left face fully fixed)", "u_clamp.txt", "C3")]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
edges = np.linspace(0.0, 1.0, 11)

for label, fname, color in CASES:
    d = np.loadtxt(fname, comments='#')
    x, y, z, ux, uy, uz = d[:,1], d[:,2], d[:,3], d[:,4], d[:,5], d[:,6]
    xs, ey, ez = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = (x >= lo) & (x < hi)
        if b.sum() < 30:
            continue
        e_yy = np.linalg.lstsq(y[b, None], uy[b], rcond=None)[0][0]
        e_zz = np.linalg.lstsq(z[b, None], uz[b], rcond=None)[0][0]
        xs.append(0.5 * (lo + hi)); ey.append(e_yy); ez.append(e_zz)
    axes[0].plot(xs, ey, 'o-', color=color, label=label)
    axes[1].plot(xs, ez, 'o-', color=color, label=label)

for ax, t in zip(axes, ["lateral strain ey vs x", "lateral strain ez vs x"]):
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel("x (distance from fixed face)")
    ax.set_ylabel(t.split()[-3])
    ax.set_title(t)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

fig.suptitle("Lateral contraction: free lateral faces vs fully clamped grip (eps_x = 100%)")
fig.tight_layout()
fig.savefig("BC_comparison.png", dpi=120)
print("wrote BC_comparison.png")
