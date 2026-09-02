from odbAccess import *
odb = openOdb('tension_125el_random.odb')
fL = odb.steps['Step-1'].frames[-1]
s = fL.fieldOutputs['S'].values
vals = [v.data[2] for v in s]
f = open('tension_125el_s33.txt', 'w')
mn = 1e30; mx = -1e30; sm = 0.0
for x in vals:
    sm += x
    if x < mn: mn = x
    if x > mx: mx = x
    f.write(str(x) + '\n')
f.close()
n = len(vals)
print('Integration pts: %d' % n)
print('S33 min: %.2f  max: %.2f  mean: %.2f MPa' % (mn, mx, sm/n))
odb.close()
