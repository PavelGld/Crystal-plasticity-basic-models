from odbAccess import *
odb = openOdb('tension_8el_random_n20.odb')
f = odb.steps['Step-1'].frames[-1]
print("SDV57 (PROPS qw) vs SDV61 (final qw) for first 8 IPs:")
for i in range(8):
    props_qw = f.fieldOutputs['SDV57'].values[i].data
    final_qw = f.fieldOutputs['SDV61'].values[i].data
    print("  IP %d: PROPS qw=%.6f  final qw=%.6f" % (i, props_qw, final_qw))
odb.close()
