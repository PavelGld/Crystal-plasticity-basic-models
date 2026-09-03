# check_odb.py -- check available field positions
from odbAccess import openOdb
import sys

odb = openOdb(sys.argv[1] + '.odb')
last_step = odb.steps.keys()[-1]
last_frame = odb.steps[last_step].frames[-1]

s33_field = last_frame.fieldOutputs['S']

positions = set()
for val in s33_field.values:
    positions.add(str(val.position))

print("Positions:", positions)
print("Num values:", len(s33_field.values))
if s33_field.values:
    v = s33_field.values[0]
    print("Example:", v.position, v.elementLabel, v.data)
odb.close()
