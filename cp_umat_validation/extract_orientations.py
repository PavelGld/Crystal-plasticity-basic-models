# extract_orientations.py -- pull lattice quaternion (SDV7..10) from last frame
# Usage: abaqus python extract_orientations.py <odb_basename> <out_txt>
from odbAccess import openOdb
import sys, math

odb_name = sys.argv[1]
out_name = sys.argv[2]

odb = openOdb(odb_name + '.odb', readOnly=True)
step = odb.steps.keys()[-1]
frame = odb.steps[step].frames[-1]
print("step=%s frame time=%s  nframes=%d" % (step, frame.frameValue, len(odb.steps[step].frames)))

sdv = frame.fieldOutputs['SDV']

# collect first-IP quaternion per element
# Abaqus SDV is 1-indexed; C++ STATEV[6..9] -> SDV7..10 -> data[6..9]
per_elem = {}
for v in sdv.values:
    el = v.elementLabel
    if el in per_elem:
        continue
    d = v.data
    per_elem[el] = (d[6], d[7], d[8], d[9])   # w, x, y, z

print("elements with SDV: %d" % len(per_elem))


def quat_to_euler(w, x, y, z):
    # normalize
    n = math.sqrt(w*w + x*x + y*y + z*z)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    w, x, y, z = w/n, x/n, y/n, z/n
    # rotation matrix from quaternion
    R = [[1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
         [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
         [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)]]
    # Bunge Euler (ZYZ), degrees
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
