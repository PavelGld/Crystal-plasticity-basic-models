from odbAccess import *
odb = openOdb('tension_8el_random_n20.odb')
f = odb.steps['Step-1'].frames[-1]
# SDV61-64 = orientation, SDV65 = NOEL
print("SDV61 (qw) for first 8 IPs:")
for i in range(8):
    v = f.fieldOutputs['SDV61'].values[i].data
    print("  IP %d: %.6f" % (i, v))
odb.close()
