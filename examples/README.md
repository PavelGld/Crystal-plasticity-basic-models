# Examples

Self-contained CPFEM runs. Each folder has the deck (or a generator), a
`run.bat`, and a `README.md`.

| # | Folder | Mesh | Grains | Loading | Purpose |
|---|---|---|---|---|---|
| 01 | [`01_single_element`](01_single_element/) | 1×`C3D8` | 1 | tension 10% | smoke test (build/link/run) |
| 02 | [`02_cube_343_random`](02_cube_343_random/) | 7×7×7 `C3D8` | 343 | tension 30%, `NLGEOM=NO` | material-per-element polycrystal |
| 03 | [`03_neper_200grains_symmetry`](03_neper_200grains_symmetry/) | Neper `C3D4` (11k) | 200 | tension 100%, symmetry BCs | real microstructure |
| 04 | [`04_neper_200grains_free`](04_neper_200grains_free/) | Neper `C3D4` (11k) | 200 | tension 100%, **free lateral faces** | honest RVE tension |

## Prerequisites

1. **Eigen 3.4+** at `D:\eigen-3.4.0` (or edit the include path).
2. **Abaqus 2022** + **Intel oneAPI** (provides the linker Abaqus needs).
3. Build the UMAT once:
   ```bash
   cd cp_umat
   build_cp_umat.bat
   ```
   Produces `cp_umat/cp_umat.obj`.

## Running an example

Every `run.bat` is self-contained (sets the environment, points at
`../../cp_umat/cp_umat.obj`):
```bash
cd examples/01_single_element
run.bat
```

## Analysis scripts

`../scripts/`:

| Script | Purpose |
|---|---|
| `neper2abaqus.py` | Neper `mesh.inp`/`mesh.ori` → Abaqus CPFEM deck (per-grain materials, BCs) |
| `extract_s33.py` | ODB → S33 per element (text) |
| `extract_quats.py` | ODB → lattice quaternion `SDV7..10` per element (text) |
| `polefigures.py` | quaternion file → `{100}/{110}/{111}` pole figures (PNG) |

### Neper → Abaqus

```bash
python scripts/neper2abaqus.py <mesh_dir> <out.inp> \
       [strain] [dt_start] [dt_max] [mesh_basename] [bc_mode]
```
- `mesh_dir` — folder with `mesh.inp` + `mesh.ori` (Neper Abaqus export)
- `bc_mode` — `sym` (1/8 symmetry) or `free` (honest tension)
- The generator validates node references and signed element volumes.

### Texture from an ODB

```bash
abaqus python scripts/extract_quats.py <job> quats.txt [frame_index]
python scripts/polefigures.py quats.txt PF.png "title" 0
```
`frame_index`: `-1` = last (default), `0` = initial. Last arg of
`polefigures.py` is the projection axis (0=X, 1=Y, 2=Z).

## Key gotchas (learned the hard way)

1. **Build the deck from Neper's `mesh.inp`, not `mesh.msh`** — the C3D10 node
   ordering differs between the two (last two nodes swapped), which makes every
   element have negative volume and Abaqus dies in the Flattener.
2. **`*UNSYMM=YES`** is required (the consistent tangent is non-symmetric).
3. **`*CONTROLS, PARAMETERS=LINE SEARCH`** — Abaqus has line search **off** by
   default; without it convergence collapses.
4. **The UMAT must not set `PNEWDT=1.0`** unconditionally — Abaqus takes the
   minimum PNEWDT as the dt multiplier, so that pins the increment size.
5. **Mesh density drives cost.** ~230 elements/grain is overkill for texture;
   ~50 gave a ~30× speed-up with identical texture.
