# -*- coding: utf-8 -*-
import random
import math

n = 7  # elements per side -> 343 elements, 2744 IPs
random.seed(42)

def random_quaternion():
    u1, u2, u3 = random.random(), random.random(), random.random()
    r1 = math.sqrt(1.0 - u1)
    r2 = math.sqrt(u1)
    t1 = 2.0 * math.pi * u2
    t2 = 2.0 * math.pi * u3
    qw = r2 * math.cos(t2)
    qx = r1 * math.sin(t1)
    qy = r1 * math.cos(t1)
    qz = r2 * math.sin(t2)
    norm = math.sqrt(qw*qw + qx*qx + qy*qy + qz*qz)
    return qw/norm, qx/norm, qy/norm, qz/norm

nx = ny = nz = n + 1

with open('cp_umat_validation/tension_343el_random.inp', 'w') as f:
    f.write('*HEADING\n7x7x7 cube (343 elements) random orientations, to eps=0.3\n')
    
    # Nodes
    f.write('*NODE\n')
    node_id = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = i / float(n)
                y = j / float(n)
                z = k / float(n)
                f.write('%d, %.4f, %.4f, %.4f\n' % (node_id, x, y, z))
                node_id += 1
    
    # Elements
    f.write('*ELEMENT, TYPE=C3D8\n')
    el = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                n1 = k*ny*nx + j*nx + i + 1
                n2 = n1 + 1
                n3 = n2 + nx
                n4 = n1 + nx
                n5 = n1 + ny*nx
                n6 = n2 + ny*nx
                n7 = n3 + ny*nx
                n8 = n4 + ny*nx
                f.write('%d, %d,%d,%d,%d, %d,%d,%d,%d\n' % (el, n1,n2,n3,n4, n5,n6,n7,n8))
                el += 1
    
    # Materials with random quaternions
    for e in range(1, n*n*n + 1):
        qw, qx, qy, qz = random_quaternion()
        f.write('*ELSET, ELSET=EL%d\n%d\n' % (e, e))
        f.write('*SOLID SECTION, ELSET=EL%d, MATERIAL=CP%d\n' % (e, e))
        f.write('*MATERIAL, NAME=CP%d\n' % e)
        f.write('*USER MATERIAL, CONSTANTS=19, UNSYMM\n')
        f.write('0, 106750, 60410, 28340, 2, 2, 1, 6\n')
        f.write('0.588235, 50, %.8f, %.8f, %.8f, %.8f, 1, 2\n' % (qw, qx, qy, qz))
        f.write('115, 34, 2.25\n')
        f.write('*DEPVAR\n70\n')
    
    # Boundary conditions — only ZERO magnitudes in model definition
    f.write('*BOUNDARY\n')
    # Bottom face z=0: U3=0
    for j in range(ny):
        for i in range(nx):
            nid = j*nx + i + 1
            f.write('%d, 3, 3, 0.0\n' % nid)
    # Anti-rotation: node at (0,0,0) fix U1=U2=0
    f.write('1, 1, 2, 0.0\n')
    # Node at (1,0,0) fix U2=0
    f.write('%d, 2, 2, 0.0\n' % (nx))
    
    # Step
    f.write('*STEP, NLGEOM=NO, INC=600, UNSYMM=YES\n')
    f.write('*STATIC\n0.0005, 0.3, 1e-8, 0.005\n')
    f.write('*CONTROLS, PARAMETERS=TIME INCREMENTATION\n')
    f.write('8, 10, 9, 16, 10, 4, 12, 5\n')
    f.write('*CONTROLS, PARAMETERS=LINE SEARCH\n')
    f.write('5, 1.0, 0.0001, 0.25, 0.10\n')
    f.write('*BOUNDARY\n')
    # Top face z=1: U3=0.3 (ramped in step)
    for j in range(ny):
        for i in range(nx):
            nid = (nz-1)*ny*nx + j*nx + i + 1
            f.write('%d, 3, 3, 0.3\n' % nid)
    f.write('*OUTPUT, FIELD\n*NODE OUTPUT\nU\n*ELEMENT OUTPUT\nS, SDV\n*END STEP\n')

print('Generated: cp_umat_validation/tension_343el_random.inp')
print('Elements: %d, IPs: %d' % (n*n*n, n*n*n*8))
