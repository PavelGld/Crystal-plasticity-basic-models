# check_volume.py -- report tr(LE)=ln(J) statistics (volumetric behaviour)
# Usage: abaqus python check_volume.py <odb_basename> [odb_basename2 ...]
from odbAccess import openOdb
import sys

for job in sys.argv[1:]:
    o = openOdb(job + '.odb', readOnly=True)
    f = o.steps[o.steps.keys()[-1]].frames[-1]
    tr = []
    for v in f.fieldOutputs['LE'].values:
        d = v.data
        tr.append(d[0] + d[1] + d[2])   # ln(J) for logarithmic strain
    n = len(tr)
    mean = sum(tr) / n
    print("%-14s ln(J): mean=%+.5f  min=%+.5f  max=%+.5f   (J_mean=%.5f, n=%d)"
          % (job, mean, min(tr), max(tr), pow(2.718281828, mean), n))
    o.close()
