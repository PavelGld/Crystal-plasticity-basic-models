# dump_u.py -- dump node coordinates + displacement from the last frame
# Usage: abaqus python dump_u.py <odb_basename> <out_txt>
from odbAccess import openOdb
import sys

odb = openOdb(sys.argv[1] + '.odb', readOnly=True)
step = odb.steps.keys()[-1]
frame = odb.steps[step].frames[-1]
inst = odb.rootAssembly.instances[odb.rootAssembly.instances.keys()[0]]

coords = {}
for n in inst.nodes:
    coords[n.label] = n.coordinates

u = frame.fieldOutputs['U']
rows = []
for v in u.values:
    c = coords[v.nodeLabel]
    rows.append((v.nodeLabel, c[0], c[1], c[2], v.data[0], v.data[1], v.data[2]))

with open(sys.argv[2], 'w') as f:
    f.write("# node x y z ux uy uz\n")
    for r in rows:
        f.write("%d %.8f %.8f %.8f %.8e %.8e %.8e\n" % r)

odb.close()
print("wrote %s (%d nodes)" % (sys.argv[2], len(rows)))
