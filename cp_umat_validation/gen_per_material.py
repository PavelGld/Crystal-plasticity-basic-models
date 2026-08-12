# -*- coding: utf-8 -*-
# Generate cube mesh with material-per-element orientations
import random, sys

n = 2   # 2x2x2 = 8 elements
seed = 42
random.seed(seed)
nx = ny = nz = n+1
n_elems = n*n*n

# Base PROPS (19 values, 8 per line)
base = [0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0,
        0.588235, 50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
        115.0, 34.0, 2.25]

def write_props(f, q):
    props = list(base)
    props[10:14] = list(q)
    f.write('*USER MATERIAL, CONSTANTS=%d, UNSYMM\n' % len(props))
    for i in range(0, len(props), 8):
        f.write(', '.join('%.8g' % v for v in props[i:i+8]) + '\n')

out = 'tension_8el_random.inp'
with open('cp_umat_validation/' + out, 'w') as f:
    f.write('*HEADING\n%d-element cube with material-per-element random orientations\n' % n_elems)
    f.write('*NODE\n')
    node_id = 1
    coords = {}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = i / float(n); y = j / float(n); z = k / float(n)
                f.write('%d, %.4f, %.4f, %.4f\n' % (node_id, x, y, z))
                coords[(i,j,k)] = node_id
                node_id += 1

    # Elements with individual ELSETs (one per element → one material per element)
    f.write('*ELEMENT, TYPE=C3D8\n')
    el = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                n1 = coords[(i,  j,  k)];  n2 = coords[(i+1,j,  k)]
                n3 = coords[(i+1,j+1,k)];  n4 = coords[(i,  j+1,k)]
                n5 = coords[(i,  j,  k+1)]; n6 = coords[(i+1,j,  k+1)]
                n7 = coords[(i+1,j+1,k+1)]; n8 = coords[(i,  j+1,k+1)]
                f.write('%d, %d,%d,%d,%d, %d,%d,%d,%d\n' % (el, n1,n2,n3,n4, n5,n6,n7,n8))
                el += 1

    # Individual ELSET per element → individual *SOLID SECTION + *MATERIAL per element
    for e in range(1, n_elems+1):
        f.write('*ELSET, ELSET=EL%d\n%d\n' % (e, e))
        f.write('*SOLID SECTION, ELSET=EL%d, MATERIAL=CP%d\n' % (e, e))

    # One *MATERIAL per element -- same base PROPS, different quaternion at [10-13]
    for e in range(1, n_elems+1):
        qw = random.gauss(0,1); qx = random.gauss(0,1)
        qy = random.gauss(0,1); qz = random.gauss(0,1)
        norm = (qw*qw+qx*qx+qy*qy+qz*qz)**0.5
        q = (qw/norm, qx/norm, qy/norm, qz/norm)  # random
        f.write('*MATERIAL, NAME=CP%d\n' % e)
        write_props(f, q)
        f.write('*DEPVAR\n70\n')

    # Boundary
    f.write('*BOUNDARY\n')
    for (i,j,k) in [(i,j,0) for i in range(nx) for j in range(ny)]:
        f.write('%d, 3, 3, 0.0\n' % coords[(i,j,k)])
    f.write('%d, 1, 2, 0.0\n' % coords[(0,0,0)])
    f.write('%d, 2, 2, 0.0\n' % coords[(n,0,0)])

    f.write('*STEP, NLGEOM=NO, INC=200\n*STATIC\n0.005, 1.0, 1e-10, 0.005\n')
    f.write('*BOUNDARY\n')
    for (i,j,k) in [(i,j,n) for i in range(nx) for j in range(ny)]:
        f.write('%d, 3, 3, 0.1\n' % coords[(i,j,k)])
    f.write('*OUTPUT, FIELD\n*NODE OUTPUT\nU\n*ELEMENT OUTPUT\nS\n*END STEP\n')

print('Generated: cp_umat_validation/' + out)
print('Elements: %d, Materials: %d, Seed: %d' % (n_elems, n_elems, seed))
