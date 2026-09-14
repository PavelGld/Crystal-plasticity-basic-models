# extract_orientations.py -- pull lattice quaternion (SDV7..10) from last frame
# Usage: abaqus python extract_orientations.py <odb_basename> <out_txt>
from odbAccess import openOdb
import sys, math

odb_name = sys.argv[1]
out_name = sys.argv[2]

odb = openOdb(odb_name + '.odb', readOnly=True)
step = odb.steps.keys()[-1]
frame = odb.steps[step].frames[-1]
fo = frame.fieldOutputs
print("step=%s frame time=%s nframes=%d" % (step, frame.frameValue, len(odb.steps[step].frames)))

# Abaqus may expose SDV either as one field 'SDV' or as separate SDV1..SDVN
if 'SDV' in fo:
    def q_of(v):
        d = v.data
        return d[6], d[7], d[8], d[9]
    src = fo['SDV']
else:
    f7, f8, f9, f10 = fo['SDV7'], fo['SDV8'], fo['SDV9'], fo['SDV10']
    v7 = dict((v.elementLabel, v.data) for v in f7.values)
    v8 = dict((v.elementLabel, v.data) for v in f8.values)
    v9 = dict((v.elementLabel, v.data) for v in f9.values)
    v10 = dict((v.elementLabel, v.data) for v in f10.values)
    def q_of(v):
        el = v.elementLabel
        return v7[el], v8[el], v9[el], v10[el]
    src = f7

# first-IP quaternion per element
per_elem = {}
for v in src.values:
    el = v.elementLabel
    if el not in per_elem:
        per_elem[el] = q_of(v)
print("elements with SDV: %d" % len(per_elem))


def quat_to_euler(w, x, y, z):
    n = math.sqrt(w*w + x*x + y*y + z*z)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    w, x, y, z = w/n, x/n, y/n, z/n
    R = [[1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
         [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
         [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)]]
    Phi = math.acos(max(-1.0, min(1.0, R[2][2])))
    if abs(math.sin(Phi)) > 1e-8:
        phi1 = math.atan2(R[2][0], -R[2][1])
        phi2 = math.atan2(R[0][2], R[1][2])
    else:
        phi1 = math.atan2(-R[0][1], R[0][0])
        phi2 = 0.0
    return (math.degrees(phi1) % 360.0, math.degrees(Phi), math.degrees(phi2) % 360.0)


with open(out_name, 'w') as f:
    f.write("# elem phi1 Phi phi2 (Bunge, deg)  [from SDV7-10]\n")
    for el in sorted(per_elem.keys()):
        w, x, y, z = per_elem[el]
        p1, P, p2 = quat_to_euler(w, x, y, z)
        f.write("%d %.6f %.6f %.6f\n" % (el, p1, P, p2))

odb.close()
print("wrote " + out_name)
