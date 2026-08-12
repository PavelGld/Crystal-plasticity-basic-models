from odbAccess import *
odb = openOdb('tension_8el_random.odb')
fL = odb.steps['Step-1'].frames[-1]
s = fL.fieldOutputs['S'].values
vals = []
for v in s:
    vals.append(v.data[2])
mn = 1e30; mx = -1e30; sm = 0.0
for v in vals:
    if v < mn: mn = v
    if v > mx: mx = v
    sm = sm + v
print('S33: min=%.2f  max=%.2f  mean=%.2f MPa' % (mn, mx, sm/len(vals)))
odb.close()
