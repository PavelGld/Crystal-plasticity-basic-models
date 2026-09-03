from odbAccess import *
odb = openOdb('tension_8el_random_n20_debug.odb')
f1 = odb.steps['Step-1'].frames[1]
print("SDV61-65 for IP 8-15 (element 2):")
for i in range(8, 16):
    qw = f1.fieldOutputs['SDV61'].values[i].data
    qx = f1.fieldOutputs['SDV62'].values[i].data
    qy = f1.fieldOutputs['SDV63'].values[i].data
    qz = f1.fieldOutputs['SDV64'].values[i].data
    noel = f1.fieldOutputs['SDV65'].values[i].data
    print("  IP %d (NOEL=%d): qw=%.6f qx=%.6f qy=%.6f qz=%.6f" % (i, int(noel), qw, qx, qy, qz))
odb.close()
