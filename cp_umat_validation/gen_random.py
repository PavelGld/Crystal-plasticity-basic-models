# Generate cube mesh with random orientations per element
import random, sys

n = 5  # elements per side
seed = 42  # reproducibility

nx, ny, nz = n+1, n+1, n+1
n_elems = n * n * n
random.seed(seed)

out = 'tension_125el_random.inp'

with open('cp_umat_validation/' + out, 'w') as f:
    f.write('*HEADING\n%d-element cube with random grain orientations\n' % n_elems)
    f.write('*NODE\n')
    node_id = 1
    coords = {}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = i / float(n)
                y = j / float(n)
                z = k / float(n)
                f.write('%d, %.4f, %.4f, %.4f\n' % (node_id, x, y, z))
                coords[(i,j,k)] = node_id
                node_id += 1

    f.write('*ELEMENT, TYPE=C3D8, ELSET=EALL\n')
    el = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                n1 = coords[(i,   j,   k)]
                n2 = coords[(i+1, j,   k)]
                n3 = coords[(i+1, j+1, k)]
                n4 = coords[(i,   j+1, k)]
                n5 = coords[(i,   j,   k+1)]
                n6 = coords[(i+1, j,   k+1)]
                n7 = coords[(i+1, j+1, k+1)]
                n8 = coords[(i,   j+1, k+1)]
                f.write('%d, %d,%d,%d,%d, %d,%d,%d,%d\n' % (el, n1,n2,n3,n4, n5,n6,n7,n8))
                el += 1

    f.write('*SOLID SECTION, ELSET=EALL, MATERIAL=CPFCC\n')
    f.write('*MATERIAL, NAME=CPFCC\n')

    # Base PROPS (19 for KM) + 1 (n_elems) + n_elems*4 (quaternions)
    f.write('*USER MATERIAL, CONSTANTS=%d\n' % (19 + 1 + n_elems * 4))

    # Line 1: dummy, c11,c12,c44, hard_mode, rot_mode, integ_mode, tc0
    f.write('0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0\n')
    # Line 2: strain_rate_infl, n_power, KA, Kq, Kh0, Kts, Kalfa, n_elems
    f.write('0.588235, 50.0, 1.0, 2.0, 115.0, 34.0, 2.25, %d\n' % n_elems)

    # Random quaternions: 6 per line (4 doubles * 6 = 24, but we do 2 quaternions = 8 values per line)
    quats = []
    for ei in range(n_elems):
        # random unit quaternion via Marsaglia method
        qw = random.gauss(0, 1)
        qx = random.gauss(0, 1)
        qy = random.gauss(0, 1)
        qz = random.gauss(0, 1)
        norm = (qw*qw + qx*qx + qy*qy + qz*qz) ** 0.5
        if norm < 1e-10:
            qw, qx, qy, qz = 1, 0, 0, 0
        else:
            qw /= norm; qx /= norm; qy /= norm; qz /= norm
        quats.append((qw, qx, qy, qz))

    # Write quaternions: 2 per line (8 values)
    idx = 0
    while idx < len(quats):
        line_parts = []
        for j in range(2):
            if idx + j < len(quats):
                q = quats[idx + j]
                line_parts.extend(['%.6f' % q[0], '%.6f' % q[1], '%.6f' % q[2], '%.6f' % q[3]])
        f.write(', '.join(line_parts) + '\n')
        idx += 2

    f.write('*DEPVAR\n70\n')

    # Boundary conditions
    f.write('*BOUNDARY\n')
    bottom = [(i, j, 0) for i in range(nx) for j in range(ny)]
    for (i,j,k) in bottom:
        f.write('%d, 3, 3, 0.0\n' % coords[(i,j,k)])

    n_corner1 = coords[(0,0,0)]
    n_corner2 = coords[(n,0,0)]
    f.write('%d, 1, 2, 0.0\n' % n_corner1)
    f.write('%d, 2, 2, 0.0\n' % n_corner2)

    f.write('*STEP, NLGEOM=NO, INC=20\n')
    f.write('*STATIC\n0.05, 1.0, 1e-08, 0.05\n')
    f.write('*BOUNDARY\n')

    top = [(i, j, n) for i in range(nx) for j in range(ny)]
    for (i,j,k) in top:
        f.write('%d, 3, 3, 0.1\n' % coords[(i,j,k)])

    f.write('*OUTPUT, FIELD\n*NODE OUTPUT\nU\n*ELEMENT OUTPUT\nS, SDV\n*END STEP\n')

print('Generated: cp_umat_validation/' + out)
print('Elements: %d, PROPS count: %d, Seed: %d' % (n_elems, 19 + 1 + n_elems*4, seed))
