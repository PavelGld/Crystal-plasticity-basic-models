from odbAccess import *
odb = openOdb('tension_8el_random_n20.odb')
step = odb.steps['Step-1']
print("Frames:", len(step.frames))
print("Last frame time:", step.frames[-1].frameValue)
print("Field keys:", list(step.frames[-1].fieldOutputs.keys()))
odb.close()
