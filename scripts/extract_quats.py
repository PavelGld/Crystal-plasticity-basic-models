# extract_quats.py -- dump lattice quaternion (SDV7..10) per element to text
# Usage: abaqus python extract_quats.py <odb_basename> <out_txt>
from odbAccess import openOdb
import sys

odb_name = sys.argv[1]
out_name = sys.argv[2]
frame_idx = int(sys.argv[3]) if len(sys.argv) > 3 else -1

odb = openOdb(odb_name + '.odb', readOnly=True)
step = odb.steps.keys()[-1]
frames = odb.steps[step].frames
frame = frames[frame_idx]
print("step=%s frame=%d of %d (time=%s)" % (step, frame_idx, len(frames), frame.frameValue))
fo = frame.fieldOutputs

if 'SDV' in fo:
    rows = [(v.elementLabel, v.data[6], v.data[7], v.data[8], v.data[9])
            for v in fo['SDV'].values]
else:
    d7 = dict((v.elementLabel, v.data) for v in fo['SDV7'].values)
    d8 = dict((v.elementLabel, v.data) for v in fo['SDV8'].values)
    d9 = dict((v.elementLabel, v.data) for v in fo['SDV9'].values)
    d10 = dict((v.elementLabel, v.data) for v in fo['SDV10'].values)
    rows = [(el, d7[el], d8[el], d9[el], d10[el]) for el in sorted(d7.keys())]

# one line per element (first IP)
seen = set()
with open(out_name, 'w') as f:
    f.write("# elem qw qx qy qz\n")
    for (el, w, x, y, z) in rows:
        if el in seen:
            continue
        seen.add(el)
        f.write("%d %.10f %.10f %.10f %.10f\n" % (el, w, x, y, z))

odb.close()
print("wrote %s (%d elements)" % (out_name, len(seen)))
