from odbAccess import *
odb = openOdb('tension_8el_same.odb')
fL = odb.steps['Step-1'].frames[-1]
s = fL.fieldOutputs['S'].values
per_ip = []
for v in s:
    per_ip.append(v.data[2])
mn = 1e30; mx = -1e30; sm = 0.0
for x in per_ip:
    if x < mn: mn = x
    if x > mx: mx = x
    sm = sm + x
print("Integration points:", len(per_ip))
print("S33 min: %.4f  max: %.4f  mean: %.4f MPa" % (mn, mx, sm/len(per_ip)))
print("Element 1 (8 IPs):", ["%.4f" % x for x in per_ip[:8]])
odb.close()
