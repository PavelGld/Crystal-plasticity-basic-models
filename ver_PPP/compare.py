# Compare our macro stress histories with the reference dataset.
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(ROOT, "verification")
NAMES = ["tension", "compression", "simple_shear", "planar", "complex"]

SANITY = {
    "tension": 0.066848179679,
    "compression": -0.0698296403372,
    "simple_shear": 0.0108104818491,
    "planar": 0.0485238270244,
    "complex": 0.00308599513024,
}

def load(path):
    return np.loadtxt(path, skiprows=1)

def main():
    print("E=0.1 / 0.3 / 0.5 / 1.0, S11 in GPa")
    print("%-13s %8s %8s %8s %8s | %10s %10s" %
          ("traj", "E", "ref", "ours", "diff", "max|d|", "rms"))
    for name in NAMES:
        ref = load(os.path.join(REF, name, "stress_history.txt"))
        ours = load(os.path.join(ROOT, "out_%s.txt" % name))
        # ours starts at E=0.01; ref has E=0 at row 0
        ref = ref[1:1 + len(ours)]
        d = ours[:, 1:7] - ref[:, 1:7]
        maxd = np.max(np.abs(d))
        rms = np.sqrt(np.mean(d ** 2))
        for E in (0.1, 0.3, 0.5, 1.0):
            i = int(round(E / 0.01)) - 1
            if 0 <= i < len(ours):
                print("%-13s %8.2f %8.4f %8.4f %+8.4f | %10.4f %10.4f" %
                      (name, E, ref[i, 1], ours[i, 1], ours[i, 1] - ref[i, 1], maxd, rms))
        # final sanity value
        s_ref = ref[-1, 1]
        s_our = ours[-1, 1]
        print("%-13s %8s %8.4f %8.4f %+8.4f | sanity %.6f (doc %.6f)" %
              ("", "final", s_ref, s_our, s_our - s_ref, SANITY[name], SANITY[name]))

if __name__ == "__main__":
    main()
