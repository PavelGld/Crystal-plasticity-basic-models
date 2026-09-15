# Example 04 — Neper 200-grain RVE, honest uniaxial tension (free lateral faces)

Same 200-grain mesh as [example 03](../03_neper_200grains_symmetry/), but with
**proper uniaxial-tension BCs**: one face is pulled, the opposite face is held,
and **the other four faces are completely free** — the specimen may contract
laterally (no symmetry constraint, no artificial stresses).

```
LEFT   (x=0)  ux = 0            <- held face
node 514 (0,0,0)   uy = uz = 0  <- remove y,z translation
node 195 (0,1,0)   uz = 0       <- remove rotation about x
RIGHT  (x=L)  ux = 1.0          <- pulled face (100% engineering strain)
(faces y=0, y=L, z=0, z=L: FREE)
```

Why these extra nodes: `ux=0` on the whole LEFT face removes x-translation and
rotations about y,z. The remaining three rigid-body modes (y-translation,
z-translation, rotation about x) are removed by a minimal 3-2-1 set — one node
with `uy=uz=0`, one further node with `uz=0` at a different `y`. Nothing else
is constrained, so the four lateral faces are free.

> Contrast with example 03, where BOTTOM (`uy=0`) and BACK (`uz=0`) are entire
> symmetry planes that prevent lateral contraction.

## Generate the deck
```bash
cd ../..
python scripts/neper2abaqus.py meshes/neper_200gr_11135el \
       examples/04_neper_200grains_free/model.inp \
       1.0 0.0001 0.001 mesh_c free
```
The last argument switches the BC mode (`sym` | `free`).

## Run
```bash
run.bat
```

## Notes
- Everything else (mesh, material, `*UNSYMM=YES`, line search, output) is
  identical to example 03, so the two can be compared directly: the free-BC run
  shows the true lateral contraction and the resulting (slightly different)
  texture and stress level.
- See example 03's README for the UMAT build and the `neper2abaqus.py` notes.
