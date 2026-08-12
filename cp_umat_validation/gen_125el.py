# Generate 5x5x5 cube mesh for CP-UMAT tension test
n = 5  # elements per side
nx, ny, nz = n+1, n+1, n+1  # nodes per side

with open('cp_umat_validation/tension_125el.inp', 'w') as f:
    f.write('*HEADING\n1000-element cube tension (5x5x5) -- CP UMAT v1\n')
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
    f.write('*USER MATERIAL, CONSTANTS=19\n')
    f.write('0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0\n')
    f.write('0.588235, 50.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0\n')
    f.write('115.0, 34.0, 2.25\n')
    f.write('*DEPVAR\n70\n')

    # Boundary: bottom face Z=0
    f.write('*BOUNDARY\n')
    bottom = [(i, j, 0) for i in range(nx) for j in range(ny)]
    for (i,j,k) in bottom:
        f.write('%d, 3, 3, 0.0\n' % coords[(i,j,k)])

    # Anti-rotation: corner node at bottom plane
    n_corner1 = coords[(0,0,0)]   # node at (0,0,0) -- pin X,Y
    n_corner2 = coords[(n,0,0)]   # node at (1,0,0) -- pin Y (stops Z-rotation)
    f.write('%d, 1, 2, 0.0\n' % n_corner1)
    f.write('%d, 2, 2, 0.0\n' % n_corner2)

    # Step
    f.write('*STEP, NLGEOM=NO, INC=200\n')
    f.write('*STATIC\n0.005, 1.0, 1e-08, 0.005\n')
    f.write('*BOUNDARY\n')

    # Top face U3 = 0.1
    top = [(i, j, n) for i in range(nx) for j in range(ny)]
    for (i,j,k) in top:
        f.write('%d, 3, 3, 0.1\n' % coords[(i,j,k)])

    f.write('*OUTPUT, FIELD\n*NODE OUTPUT\nU\n*ELEMENT OUTPUT\nS, SDV\n*END STEP\n')

print('Generated: cp_umat_validation/tension_125el.inp')
print('Nodes: %d, Elements: %d' % (nx*ny*nz, n*n*n))
