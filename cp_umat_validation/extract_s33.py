# extract_s33.py -- extract S33 (averaged over IPs) per element
# Usage: abaqus python extract_s33.py <odb_basename> <out_txt>
from odbAccess import openOdb
import sys

odb_name = sys.argv[1]
out_name = sys.argv[2]

odb = openOdb(odb_name + '.odb')
last_step = odb.steps.keys()[-1]
last_frame = odb.steps[last_step].frames[-1]

s33_field = last_frame.fieldOutputs['S']

# Accumulate per element
elem_s33 = {}
elem_count = {}
for val in s33_field.values:
    el = val.elementLabel
    if el not in elem_s33:
        elem_s33[el] = 0.0
        elem_count[el] = 0
    elem_s33[el] += val.data[2]  # S33 = data[2]
    elem_count[el] += 1

with open(out_name, 'w') as f:
    f.write("# elem_label S33_avg\n")
    for el in sorted(elem_s33.keys()):
        avg = elem_s33[el] / elem_count[el]
        f.write("%d %.10e\n" % (el, avg))

odb.close()
print("Wrote " + out_name + " (%d elements)" % len(elem_s33))
