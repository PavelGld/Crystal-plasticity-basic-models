# Example 03 — Neper 200-grain RVE, symmetry BCs

A real polycrystal from Neper: **200 grains, 11135 `C3D4` elements**, one
`*MATERIAL` per grain (Bunge Euler angles from `mesh_c.ori` → quaternion).
Uniaxial tension to 100% strain (`NLGEOM=YES`).

**BCs — 1/8 symmetric model:** three faces lie on symmetry planes, the opposite
face is pulled:

```
LEFT   (x=0)  ux = 0
BOTTOM (y=0)  uy = 0
BACK   (z=0)  uz = 0
RIGHT  (x=L)  ux = 1.0   (100% engineering strain over the step)
```

> For a tension test where the four lateral faces are **free** (no symmetry
> constraint), see [example 04](../04_neper_200grains_free/).

## Mesh
`../../meshes/neper_200gr_11135el/` — files `mesh_c.inp`, `mesh_c.ori`
(Neper's Abaqus export; **use the `.inp`, not the `.msh`** — node ordering of
C3D10 differs between them).

## Generate the deck
```bash
cd ../..
python scripts/neper2abaqus.py meshes/neper_200gr_11135el \
       examples/03_neper_200grains_symmetry/model.inp \
       1.0 0.0001 0.001 mesh_c sym
```
Args: `<mesh_dir> <out.inp> [strain] [dt_start] [dt_max] [mesh_basename] [bc_mode]`.

## Run
```bash
run.bat
```
Wall time on 4 threads: **~8 min** (1007 increments, `dt` grows automatically
to `dt_max=1e-3`).

## Results — FCC tension texture
`{111}` pole figure (projection axis = X, the loading direction) at 100% strain:
random at 5% → sharp `<111> ‖ X` at 100% (classic FCC tension).

```
abaqus python ../../scripts/extract_quats.py model quats.txt
python ../../scripts/polefigures.py quats.txt PF.png "200 grains, 100%" 0
```
(`0` = projection axis X.)

## Notes
- **`*UNSYMM=YES`** is required (non-symmetric consistent tangent).
- **`*CONTROLS, PARAMETERS=LINE SEARCH`** — Abaqus has line search **off** by
  default; without it the first increment fails and `dt` collapses.
- The UMAT must **not** pin `PNEWDT=1.0`, otherwise Abaqus never grows `dt`.
