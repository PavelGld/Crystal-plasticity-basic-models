#!/usr/bin/env python
"""polefigures.py -- pole figures {100},{110},{111} from a quaternion file.

Usage: python polefigures.py <quats_txt> <out_png> [title]

Reads:  "elem qw qx qy qz" per line (see extract_quats.py)
Writes: stereographic pole figures (upper hemisphere), scatter density.
"""
import sys
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def quat_to_R(w, x, y, z):
    n = math.sqrt(w*w + x*x + y*y + z*z)
    if n < 1e-12:
        return np.eye(3)
    w, x, y, z = w/n, x/n, y/n, z/n
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
        [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)]])


def stereographic(v, axis=2):
    # project about sample 'axis' (0=x,1=y,2=z): make that axis the polar one
    v = np.array([v[(axis+1) % 3], v[(axis+2) % 3], v[axis]])
    if v[2] < 0:
        v = -v
    d = 1.0 + v[2]
    if d < 1e-9:
        return None
    return v[0]/d, v[1]/d


def main():
    qfile = sys.argv[1]
    out_png = sys.argv[2]
    title = sys.argv[3] if len(sys.argv) > 3 else qfile
    axis = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    axis_name = {0: "X (loading)", 1: "Y", 2: "Z"}.get(axis, "?")

    quats = []
    for line in open(qfile):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        p = s.split()
        if len(p) >= 5:
            quats.append((float(p[1]), float(p[2]), float(p[3]), float(p[4])))

    poles = [((1, 0, 0), "{100}"), ((1, 1, 0), "{110}"), ((1, 1, 1), "{111}")]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for ax, (hk, name) in zip(axes, poles):
        vc = np.array(hk, dtype=float)
        vc /= np.linalg.norm(vc)
        xs, ys = [], []
        for (w, x, y, z) in quats:
            R = quat_to_R(w, x, y, z)
            vl = R.dot(vc)
            pt = stereographic(vl, axis)
            if pt:
                xs.append(pt[0]); ys.append(pt[1])
        ax.scatter(xs, ys, s=2.5, alpha=0.22, c='C0', edgecolors='none')
        ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color='k', lw=1.1))
        ax.axhline(0, color='gray', lw=0.4)
        ax.axvline(0, color='gray', lw=0.4)
        ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
        ax.set_aspect('equal')
        ax.set_title(name)
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(title + "  |  projection axis: " + axis_name)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    print("wrote %s (%d orientations)" % (out_png, len(quats)))


if __name__ == "__main__":
    main()
