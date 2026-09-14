#!/usr/bin/env python
"""neper2abaqus.py -- build a complete CPFEM deck from Neper's Abaqus export.

Usage: python neper2abaqus.py <neper_dir> <output.inp> [max_strain]

Uses mesh.inp (already correct Abaqus C3D10 node ordering + polyN elsets)
and mesh.ori (Bunge Euler angles per grain).  Emits a FLAT deck:
  *NODE, *ELEMENT, per-grain *ELSET + *SOLID SECTION + *MATERIAL,
  boundary NSETs, *STEP.
"""
import sys
import os
import re
import numpy as np


def parse_mesh_inp(path):
    nodes = {}
    elements = {}          # eid -> tuple(nodes)
    elsets = {}            # name -> [eids]
    section = None
    cur_elset = None
    cur_data = []

    def flush():
        if cur_elset and cur_data:
            elsets[cur_elset] = list(cur_data)

    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('**'):
                continue
            if line.startswith('*'):
                flush()
                cur_elset = None
                cur_data = []
                low = line.lower()
                if low.startswith('*node'):
                    section = 'node'
                elif low.startswith('*element'):
                    section = 'elem'
                elif low.startswith('*elset'):
                    m = re.search(r'elset\s*=\s*([^,]+)', line, re.I)
                    name = m.group(1).strip() if m else None
                    if name and name.lower().startswith('poly'):
                        cur_elset = name
                        section = 'elset'
                    else:
                        section = None
                else:
                    section = None
                continue
            if section == 'node':
                p = line.replace(',', ' ').split()
                if len(p) >= 4 and p[0].lstrip('-').isdigit():
                    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            elif section == 'elem':
                p = line.replace(',', ' ').split()
                if len(p) >= 2 and p[0].lstrip('-').isdigit():
                    elements[int(p[0])] = tuple(int(x) for x in p[1:])
            elif section == 'elset':
                for tok in line.replace(',', ' ').split():
                    if tok.lstrip('-').isdigit():
                        cur_data.append(int(tok))
    flush()
    return nodes, elements, elsets


def parse_ori(path):
    out = {}
    inori = False
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if s == '$ElsetOrientations':
                inori = True
                continue
            if s == '$EndElsetOrientations':
                break
            if inori and s and not s.startswith('$'):
                p = s.split()
                if len(p) >= 4 and p[0].lstrip('-').isdigit():
                    out[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    return out


def euler_to_quaternion(phi1, Phi, phi2):
    p1, P, p2 = np.radians(phi1), np.radians(Phi), np.radians(phi2)
    c1, s1 = np.cos(p1), np.sin(p1)
    c2, s2 = np.cos(P), np.sin(P)
    c3, s3 = np.cos(p2), np.sin(p2)
    R = np.array([
        [c1*c3 - s1*s2*s3, -c1*s3 - s1*c2*c3, s1*s2],
        [s1*c3 + c1*s2*s3, -s1*s3 + c1*c2*c3, -c1*s2],
        [s2*s3, c2*s3, c2]])
    tr = R[0,0] + R[1,1] + R[2,2]
    if tr > 0:
        s = 0.5 / np.sqrt(tr + 1.0); w = 0.25 / s
        x = (R[2,1]-R[1,2])*s; y = (R[0,2]-R[2,0])*s; z = (R[1,0]-R[0,1])*s
    elif R[0,0] > R[1,1] and R[0,0] > R[2,2]:
        s = 2.0*np.sqrt(1.0+R[0,0]-R[1,1]-R[2,2]); w = (R[2,1]-R[1,2])/s
        x = 0.25*s; y = (R[0,1]+R[1,0])/s; z = (R[0,2]+R[2,0])/s
    elif R[1,1] > R[2,2]:
        s = 2.0*np.sqrt(1.0+R[1,1]-R[0,0]-R[2,2]); w = (R[0,2]-R[2,0])/s
        x = (R[0,1]+R[1,0])/s; y = 0.25*s; z = (R[1,2]+R[2,1])/s
    else:
        s = 2.0*np.sqrt(1.0+R[2,2]-R[0,0]-R[1,1]); w = (R[1,0]-R[0,1])/s
        x = (R[0,2]+R[2,0])/s; y = (R[1,2]+R[2,1])/s; z = 0.25*s
    return w, x, y, z


def chunk(f, ids, per=16):
    ids = list(ids)
    for i in range(0, len(ids), per):
        f.write(", ".join(str(x) for x in ids[i:i+per]) + "\n")


def main():
    neper_dir, out = sys.argv[1], sys.argv[2]
    strain = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    dt = float(sys.argv[4]) if len(sys.argv) > 4 else 0.001
    dtmax = float(sys.argv[5]) if len(sys.argv) > 5 else dt
    base = os.path.join(neper_dir, "mesh")

    nodes, elements, elsets = parse_mesh_inp(base + ".inp")
    ori = parse_ori(base + ".ori")
    print("nodes=%d elements=%d poly-elsets=%d orientations=%d"
          % (len(nodes), len(elements), len(elsets), len(ori)))

    # validity check: every element node must exist; every elset elem must exist
    missing_nodes = sum(1 for c in elements.values() for n in c if n not in nodes)
    missing_elems = sum(1 for ids in elsets.values() for e in ids if e not in elements)
    print("missing node refs=%d, missing elem refs=%d" % (missing_nodes, missing_elems))

    # element volume sanity (signed) for C3D10: use corner tet
    def tet_vol(c):
        a = np.array(nodes[c[0]]); b = np.array(nodes[c[1]])
        d = np.array(nodes[c[2]]); e = np.array(nodes[c[3]])
        return np.dot(b-a, np.cross(d-a, e-a)) / 6.0
    vols = [tet_vol(c) for c in elements.values() if len(c) >= 4]
    neg = sum(1 for v in vols if v <= 0)
    print("negative/zero-volume corner tets: %d / %d" % (neg, len(vols)))

    xs = [c[0] for c in nodes.values()]
    ys = [c[1] for c in nodes.values()]
    zs = [c[2] for c in nodes.values()]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    minz, maxz = min(zs), max(zs)
    tol = 1e-3 * max(maxx-minx, maxy-miny, maxz-minz)

    left = [n for n, c in nodes.items() if abs(c[0]-minx) < tol]
    right = [n for n, c in nodes.items() if abs(c[0]-maxx) < tol]
    bottom = [n for n, c in nodes.items() if abs(c[1]-miny) < tol]
    back = [n for n, c in nodes.items() if abs(c[2]-minz) < tol]

    grain_list = sorted(int(m.group(1)) for m in
                        (re.match(r'poly(\d+)', k, re.I) for k in elsets) if m)

    with open(out, 'w') as f:
        f.write("*HEADING\nCPFEM Neper %d grains, eps=%.2f\n" % (len(grain_list), strain))

        f.write("*NODE\n")
        for nid, (x, y, z) in sorted(nodes.items()):
            f.write("%d, %.10e, %.10e, %.10e\n" % (nid, x, y, z))

        f.write("*ELEMENT, TYPE=C3D10, ELSET=ALLEL\n")
        for eid, c in sorted(elements.items()):
            f.write("%d, %s\n" % (eid, ", ".join(str(n) for n in c)))

        for g in grain_list:
            elset = "poly%d" % g
            mat = "GRAIN_%d" % g
            f.write("*ELSET, ELSET=%s\n" % elset)
            chunk(f, elsets[elset])
            f.write("*SOLID SECTION, ELSET=%s, MATERIAL=%s\n" % (elset, mat))

        for g in grain_list:
            mat = "GRAIN_%d" % g
            if g in ori:
                w, x, y, z = euler_to_quaternion(*ori[g])
            else:
                w, x, y, z = 1.0, 0.0, 0.0, 0.0
            f.write("*MATERIAL, NAME=%s\n" % mat)
            f.write("*USER MATERIAL, CONSTANTS=19, UNSYMM\n")
            f.write("0, 106750, 60410, 28340, 2, 2, 1, 6\n")
            f.write("0.588235, 50, %.10f, %.10f, %.10f, %.10f, 1, 2\n" % (w, x, y, z))
            f.write("115, 34, 2.25\n")
            f.write("*DEPVAR\n70\n")

        f.write("*NSET, NSET=LEFT\n"); chunk(f, left)
        f.write("*NSET, NSET=RIGHT\n"); chunk(f, right)
        f.write("*NSET, NSET=BOTTOM\n"); chunk(f, bottom)
        f.write("*NSET, NSET=BACK\n"); chunk(f, back)

        # Step 1: small transient past initial yield; Step 2: steady large step
        s0 = min(0.005, strain)
        f.write("*STEP, NLGEOM=YES, INC=200000, UNSYMM=YES\n")
        f.write("*STATIC\n%.6g, %.6g, 1e-9, %.6g\n" % (dt, s0, dtmax))
        f.write("*BOUNDARY\n")
        f.write("LEFT, 1, 1, 0.0\nBOTTOM, 2, 2, 0.0\nBACK, 3, 3, 0.0\n")
        f.write("*BOUNDARY, TYPE=DISPLACEMENT\n")
        f.write("RIGHT, 1, 1, %.6f\n" % (s0 * (maxx - minx)))
        f.write("*OUTPUT, FIELD, FREQUENCY=50\n")
        f.write("*ELEMENT OUTPUT\nS, E, LE, SDV\n")
        f.write("*NODE OUTPUT\nU\n")
        f.write("*END STEP\n")

        if strain > s0:
            f.write("*STEP, NLGEOM=YES, INC=200000, UNSYMM=YES\n")
            f.write("*STATIC\n%.6g, %.6g, 1e-9, %.6g\n" % (dtmax, strain - s0, dtmax))
            f.write("*BOUNDARY, TYPE=DISPLACEMENT\n")
            f.write("RIGHT, 1, 1, %.6f\n" % (strain * (maxx - minx)))
            f.write("*OUTPUT, FIELD, FREQUENCY=50\n")
            f.write("*ELEMENT OUTPUT\nS, E, LE, SDV\n")
            f.write("*NODE OUTPUT\nU\n")
            f.write("*END STEP\n")
    print("wrote %s (strain=%.3f, dt=%.4g, ~%d increments)" % (out, strain, dt, int(round(strain / dt))))


if __name__ == "__main__":
    main()
