from odbAccess import *
odb = openOdb('tension_8el.odb')
fL = odb.steps['Step-1'].frames[-1]
s = fL.fieldOutputs['S'].values
print('Elements:', len(odb.rootAssembly.instances['PART-1-1'].elements))
print('Integration pts:', len(s))
for i in range(min(3, len(s))):
    d = s[i].data
    print('  %d: S33=%.2f S12=%.2f' % (i, d[2], d[3]))
print('  ...')
d = s[-1].data
print('  last: S33=%.2f MPa' % d[2])
odb.close()
