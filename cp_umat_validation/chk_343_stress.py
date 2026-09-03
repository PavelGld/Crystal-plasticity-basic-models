from odbAccess import *
odb = openOdb('tension_343el_random.odb')
f = odb.steps['Step-1'].frames[-1]
s = f.fieldOutputs['S']

# Collect all S33 values
all_s33 = []
for v in s.values:
    val = v.data[2]
    all_s33.append(val)

# Check type
print("Type of first element:", type(all_s33[0]))
print("First few values:", all_s33[:3])

# Per element average (8 IP per C3D8)
per_el = []
for el in range(343):
    vals = all_s33[el*8:(el+1)*8]
    s = 0.0
    for x in vals:
        s += float(x)
    avg = s / 8.0
    per_el.append(avg)

mn = 1e30; mx = -1e30; sm = 0.0
for x in per_el:
    x = float(x)
    if x < mn: mn = x
    if x > mx: mx = x
    sm += x
mean = sm / len(per_el)
sm2 = 0.0
for x in per_el:
    dx = float(x) - mean
    sm2 += dx*dx
var = sm2 / len(per_el)
std = var**0.5
std = var**0.5
print("S33 per element (343 elements):")
print("  min=%.2f  max=%.2f  mean=%.2f  std=%.2f MPa" % (mn, mx, mean, std))
print("  range: %.2f - %.2f MPa (%.1f%% variation)" % (mn, mx, (mx-mn)/mean*100))
odb.close()
