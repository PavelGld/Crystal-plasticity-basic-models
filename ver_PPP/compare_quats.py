# Compare end-of-trajectory grain orientations, up to quaternion sign.
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))

def load_q(path):
    a = np.loadtxt(path, skiprows=1)
    q = a[:, 1:5].copy()
    q /= np.linalg.norm(q, axis=1)[:, None]
    return a[:, 0].astype(int), q

def angles(qr, qo):
    d = np.abs(np.sum(qr * qo, axis=1))
    d = np.clip(d, 0.0, 1.0)
    return 2.0 * np.degrees(np.arccos(d))

def main(name, ours):
    _, qr = load_q(os.path.join(ROOT, "verification", name, "quaternions_end.txt"))
    _, qo = load_q(ours)
    n = min(len(qr), len(qo))
    a = angles(qr[:n], qo[:n])
    ac = angles(qr[:n], (qo[:n] * np.array([1, -1, -1, -1])) )
    print("%s: grains %d" % (name, n))
    print("  direct   mean %.3f deg  median %.3f  max %.3f" % (a.mean(), np.median(a), a.max()))
    print("  conjugate mean %.3f deg  median %.3f  max %.3f" % (ac.mean(), np.median(ac), ac.max()))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
