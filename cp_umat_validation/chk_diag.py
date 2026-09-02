from odbAccess import *
odb = openOdb('tension_random.odb')
f1 = odb.steps['Step-1'].frames[1]      # after 1 increment
fL = odb.steps['Step-1'].frames[-1]     # last converged
print("=== After 1 increment ===")
for i in [58,59,60,69]:
    v = f1.fieldOutputs['SDV%d' % (i+1)].values[0].data
    print("  STATEV[%d] = %s" % (i, v))
print("=== Last converged increment ===")
for i in [58,59,60,69]:
    v = fL.fieldOutputs['SDV%d' % (i+1)].values[0].data
    print("  STATEV[%d] = %s" % (i, v))
odb.close()
