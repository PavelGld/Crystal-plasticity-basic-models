from odbAccess import *
odb = openOdb('tension_random.odb')
fL = odb.steps['Step-1'].frames[-1]
s = fL.fieldOutputs['S'].values
mn = 1e30; mx = -1e30; sm = 0.0
for v in s:
    d = v.data[2]
    if d < mn: mn = d
    if d > mx: mx = d
    sm = sm + d
print('S33: min=%.2f max=%.2f mean=%.2f MPa' % (mn, mx, sm/len(s)))
odb.close()
