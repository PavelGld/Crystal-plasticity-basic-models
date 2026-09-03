from odbAccess import *
odb = openOdb('tension_8el_random.odb')
f1 = odb.steps['Step-1'].frames[1]
keys = [k for k in f1.fieldOutputs.keys() if k.startswith('SDV')]
print("SDV keys:", keys[:5])
odb.close()
