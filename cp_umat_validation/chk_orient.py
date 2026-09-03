from odbAccess import *
odb = openOdb('tension_8el_random.odb')
f1 = odb.steps['Step-1'].frames[1]
# SDV61-64 = orientation quaternion, SDV65 = element number
print("IP  NOEL   qw        qx        qy        qz")
for i in range(16):
    sd = f1.fieldOutputs['SDV%d' % (i+1)].values[0].data
    print("%2d  %3d  %8.4f %8.4f %8.4f %8.4f" % (i, int(sd[64]), sd[60], sd[61], sd[62], sd[63]))
odb.close()
