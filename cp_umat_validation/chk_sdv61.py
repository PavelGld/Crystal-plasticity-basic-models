from odbAccess import *
odb = openOdb('tension_8el_random_n20_debug.odb')
f1 = odb.steps['Step-1'].frames[1]
print("Frame 1 (after increment 1):")
for i in range(8):
    vals = [f1.fieldOutputs['SDV%d' % (60+j)].values[i].data for j in range(4)]
    print("  IP %d: qw=%.6f qx=%.6f qy=%.6f qz=%.6f" % (i, vals[0], vals[1], vals[2], vals[3]))
odb.close()
