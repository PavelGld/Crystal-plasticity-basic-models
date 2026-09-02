# 8 elements, all SAME non-identity orientation, EXPLICIT scheme (integ_mode=0)
n = 2
nx = ny = nz = n+1
n_elems = n*n*n

# integ_mode=0 (explicit), same quaternion for all
base = [0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 0.0, 6.0,  # integ_mode=0
        0.588235, 50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
        115.0, 34.0, 2.25]
q_same = (0.1933, 0.2319, 0.1493, -0.9417)

def write_props(f, q):
    props = list(base)
    props[10:14] = list(q)
    f.write('*USER MATERIAL, CONSTANTS=%d, UNSYMM\n' % len(props))
    for i in range(0, len(props), 8):
        f.write(', '.join('%.8g' % v for v in props[i:i+8]) + '\n')

out = 'tension_8el_same_explicit.inp'
with open('cp_umat_validation/' + out, 'w') as f:
    f.write('*HEADING\n8-element cube, all same orientation, EXPLICIT scheme\n')
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

    for e in range(1, n_elems+1):
        f.write('*ELSET, ELSET=EL%d\n%d\n' % (e, e))
        f.write('*SOLID SECTION, ELSET=EL%d, MATERIAL=CP%d\n' % (e, e))

    for e in range(1, n_elems+1):
        f.write('*MATERIAL, NAME=CP%d\n' % e)
        write_props(f, q_same)
        f.write('*DEPVAR\n70\n')

    f.write('*BOUNDARY\n')
    for (i,j,k) in [(i,j,0) for i in range(nx) for j in range(ny)]:
        f.write('%d, 3, 3, 0.0\n' % coords[(i,j,k)])
    f.write('%d, 1, 2, 0.0\n' % coords[(0,0,0)])
    f.write('%d, 2, 2, 0.0\n' % coords[(n,0,0)])

    f.write('*STEP, NLGEOM=NO, INC=200\n*STATIC\n0.001, 1.0, 1e-08, 0.001\n')
    f.write('*BOUNDARY\n')
    for (i,j,k) in [(i,j,n) for i in range(nx) for j in range(ny)]:
        f.write('%d, 3, 3, 0.1\n' % coords[(i,j,k)])
    f.write('*OUTPUT, FIELD\n*NODE OUTPUT\nU\n*ELEMENT OUTPUT\nS\n*END STEP\n')

print('Generated: cp_umat_validation/' + out)