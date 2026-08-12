import random
random.seed(42)
qw,qx,qy,qz = random.gauss(0,1),random.gauss(0,1),random.gauss(0,1),random.gauss(0,1)
n = (qw*qw+qx*qx+qy*qy+qz*qz)**0.5
qw/=n; qx/=n; qy/=n; qz/=n
print('Quaternion: %.4f %.4f %.4f %.4f' % (qw,qx,qy,qz))

with open('cp_umat_validation/tension_random.inp','w') as f:
    f.write("""*HEADING
Single element tension with random initial orientation
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0
5, 0.0, 0.0, 1.0
6, 1.0, 0.0, 1.0
7, 1.0, 1.0, 1.0
8, 0.0, 1.0, 1.0
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1,2,3,4, 5,6,7,8
*SOLID SECTION, ELSET=EALL, MATERIAL=CPFCC
*MATERIAL, NAME=CPFCC
*USER MATERIAL, CONSTANTS=19
0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0
0.588235, 50.0, %.4f, %.4f, %.4f, %.4f, 1.0, 2.0
115.0, 34.0, 2.25
*DEPVAR
70
*BOUNDARY
1, 3, 3, 0.0
2, 3, 3, 0.0
3, 3, 3, 0.0
4, 3, 3, 0.0
1, 1, 2, 0.0
2, 2, 2, 0.0
*STEP, NLGEOM=NO, INC=200
*STATIC
0.005, 1.0, 1e-08, 0.005
*BOUNDARY
5, 3, 3, 0.1
6, 3, 3, 0.1
7, 3, 3, 0.1
8, 3, 3, 0.1
*OUTPUT, FIELD
*NODE OUTPUT
U
*ELEMENT OUTPUT
S, SDV
*END STEP
""" % (qw,qx,qy,qz))
print('Generated: cp_umat_validation/tension_random.inp')
