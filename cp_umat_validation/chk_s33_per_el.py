from odbAccess import *
odb = openOdb('tension_8el_random_n20_debug.odb')
f = odb.steps['Step-1'].frames[-1]
print("S33 per element (average over 8 IPs):")
for el in range(8):
    s33 = 0
    for ip in range(8):
        idx = el * 8 + ip
        s33 += f.fieldOutputs['S'].values[idx].data[2] / 8.0
    print("  Element %d: S33 = %.2f MPa" % (el + 1, s33))
odb.close()
