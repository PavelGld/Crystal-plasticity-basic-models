from odbAccess import openOdb

odb = openOdb('model.odb')
step = odb.steps.values()[0]
frame = step.frames[-1]
print('last frame time = %.6f, frames = %d' % (frame.frameValue, len(step.frames)))

for name in ('S', 'E', 'LE'):
    if name in frame.fieldOutputs:
        v = frame.fieldOutputs[name].values[0]
        print('%s = %s' % (name, ['%.6e' % x for x in v.data]))

if 'SDV' in frame.fieldOutputs:
    sdv = frame.fieldOutputs['SDV'].values[0].data
    print('SDV[0:6]  = %s' % ['%.6e' % x for x in sdv[0:6]])
    print('SDV[6:10] = %s' % ['%.6e' % x for x in sdv[6:10]])
    print('SDV[69]   = %.6e' % sdv[69])

reg = step.historyRegions.values()[0]
for k in ('ALLSE', 'ALLPD', 'ALLWK', 'ALLIE', 'ALLKE', 'ETOTAL'):
    if k in reg.historyOutputs:
        d = reg.historyOutputs[k].data
        print('%s = %.9e' % (k, d[-1][1]))

odb.close()
