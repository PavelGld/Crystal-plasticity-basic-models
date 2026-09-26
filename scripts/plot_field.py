# plot_field.py -- scatter plot of an element stress component over the mesh
# Usage: abaqus python plot_field.py <odb_basename> <out_png> <component> [title]
#   component index is the ODB 'S' order: 0=S11 1=S22 2=S33 3=S12 4=S13 5=S23
# Builds element centroids from node coordinates, averages the component over IPs,
# and makes a filled scatter (spatial distribution of the field).
from odbAccess import openOdb
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d   # register the '3d' projection

odb_name  = sys.argv[1]
out_png   = sys.argv[2]
comp      = int(sys.argv[3]) if len(sys.argv) > 3 else 2
title     = sys.argv[4] if len(sys.argv) > 4 else odb_name

odb = openOdb(odb_name + '.odb', readOnly=True)
step = odb.steps.keys()[-1]
frame = odb.steps[step].frames[-1]
inst = odb.rootAssembly.instances[odb.rootAssembly.instances.keys()[0]]

coords = {}
for n in inst.nodes:
    coords[n.label] = n.coordinates

# element centroid = mean of its node coordinates
centroid = {}
for e in inst.elements:
    conn = e.connectivity
    c = [0.0, 0.0, 0.0]
    for nl in conn:
        nc = coords[nl]
        for k in range(3):
            c[k] += nc[k]
    for k in range(3):
        c[k] /= len(conn)
    centroid[e.label] = c

field = frame.fieldOutputs['S']
acc = {}
cnt = {}
for v in field.values:
    el = v.elementLabel
    acc[el] = acc.get(el, 0.0) + v.data[comp]
    cnt[el] = cnt.get(el, 0) + 1

have = [el for el in acc if el in centroid]
if not have:
    print("ERROR: no element labels overlap between field and instance")
    print("  field sample:", sorted(acc.keys())[:5])
    print("  centroid sample:", sorted(centroid.keys())[:5])
else:
    xs = [centroid[el][0] for el in have]
    ys = [centroid[el][1] for el in have]
    zs = [centroid[el][2] for el in have]
    vals = [acc[el] / cnt[el] for el in have]   # ODB S is in model units (MPa)

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection='3d')
comp_name = {0: 'S11', 1: 'S22', 2: 'S33', 3: 'S12', 4: 'S13', 5: 'S23'}[comp]
sc = ax.scatter(xs, ys, zs, c=vals, cmap='jet', s=18)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title(title)
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
cb.set_label(comp_name + ' (MPa)')
fig.tight_layout()
fig.savefig(out_png, dpi=130)
print("wrote %s (%d elements, %s range %.1f..%.1f MPa)" %
      (out_png, len(vals), comp_name, min(vals), max(vals)))
odb.close()
