#Output files (odb, msg, sta) go to the working directory — project root.
Example: D:\Pet_Projects\EigenProject\tension.odb.

 CP-UMAT — Crystal Plasticity UMAT for Abaqus/Standard

Implementation of crystal plasticity (FCC, Hutchinson power-law, implicit Newton–Raphson integration with consistent tangent) as an Abaqus/Standard UMAT subroutine.

Single grain per integration point. NLGEOM=NO (small strain).

## Build (Windows)

Requires: MSVC 2022 Community, Eigen 3.4 at `D:\eigen-3.4.0`.

```bat
cd cp_umat
build_cp_umat.bat
```

Produces `cp_umat.obj` — linkable into Abaqus via `user=cp_umat\cp_umat.obj`.

## Run

```bat
abaqus job=tension user=cp_umat\cp_umat.obj
```

## PROPS layout

`*USER MATERIAL, CONSTANTS=N` where N depends on hard_mode:

| PROPS[i] | Meaning |
|----------|---------|
| 0 | unused (reserved) |
| 1-3 | c11, c12, c44 (Pa) |
| 4 | hard_mode: 0=none, 1=isotropic, 2=Kocks-Mecking |
| 5 | rot_mode: 0=none, 1=Mandel, 2=Taylor |
| 6 | integ_mode: 0=explicit, 1=implicit (recommended) |
| 7 | tc0 — initial critical resolved shear stress (Pa) |
| 8 | strain_rate_infl — γ̇₀ reference strain rate |
| 9 | n_power — rate sensitivity exponent |
| 10-13 | initial quaternion (qw, qx, qy, qz) |
| 14-18 | (KM only) A, q, h0, ts, alfa |
| 14-15 | (iso only) a, b |

Constants for NPROPS:
- hard_mode=0: 14
- hard_mode=1: 16
- hard_mode=2: 19

**Important:** Abaqus reads at most 8 values per data line. Split PROPS into multiple lines.

Example (FCC Cu, KM hardening, implicit):
```
*USER MATERIAL, CONSTANTS=19
0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0
0.588235, 50.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0
115.0, 34.0, 2.25ticity UMAT for Abaqus/Standard

Implementation of crystal plasticity (FCC, Hutchinson power-law, implicit Newton–Raphson integration with consistent tangent) as an Abaqus/Standard UMAT subroutine.

Single grain per integration point. NLGEOM=NO (small strain).

## Build (Windows)

Requires: MSVC 2022 Community, Eigen 3.4 at `D:\eigen-3.4.0`.

```bat
cd cp_umat
build_cp_umat.bat
```

Produces `cp_umat.obj` — linkable into Abaqus via `user=cp_umat\cp_umat.obj`.

## Run

```bat
abaqus job=tension user=cp_umat\cp_umat.obj
```

## PROPS layout

`*USER MATERIAL, CONSTANTS=N` where N depends on hard_mode:

| PROPS[i] | Meaning |
|----------|---------|
| 0 | unused (reserved) |
| 1-3 | c11, c12, c44 (Pa) |
| 4 | hard_mode: 0=none, 1=isotropic, 2=Kocks-Mecking |
| 5 | rot_mode: 0=none, 1=Mandel, 2=Taylor |
| 6 | integ_mode: 0=explicit, 1=implicit (recommended) |
| 7 | tc0 — initial critical resolved shear stress (Pa) |
| 8 | strain_rate_infl — γ̇₀ reference strain rate |
| 9 | n_power — rate sensitivity exponent |
| 10-13 | initial quaternion (qw, qx, qy, qz) |
| 14-18 | (KM only) A, q, h0, ts, alfa |
| 14-15 | (iso only) a, b |

Constants for NPROPS:
- hard_mode=0: 14
- hard_mode=1: 16
- hard_mode=2: 19

**Important:** Abaqus reads at most 8 values per data line. Split PROPS into multiple lines.

Example (FCC Cu, KM hardening, implicit):
```
*USER MATERIAL, CONSTANTS=19
0.0, 106750.0e6, 60410.0e6, 28340.0e6, 2.0, 2.0, 1.0, 6.0e6, 0.588235, 50.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 115.0e6, 34.0e6, 2.25
```

## STATEV layout (NSTATV = 70)

| STATEV[i] | Content |
|-----------|---------|
| 0-5 | stress_v (Voigt, crystal frame) |
| 6-9 | quaternion (w,x,y,z) |
| 10-33 | tc[0..23] — critical resolved shear stress per slip system |
| 34-57 | shifts[0..23] — accumulated slip per system |
| 58-66 | (reserved) |
| 67 | sim_time |
| 68 | acc_def — accumulated equivalent strain |
| 69 | guard (123456.0 = initialized) |

## Slip systems

FCC {111}⟨110⟩ — 12 systems (doubled to 24 via ±b). Hardcoded (not in PROPS).

## Validation files

cp_umat_validation/ contains test cases:

| File | Description |
|------|-------------|
| 	ension.inp | Single C3D8 element, uniaxial tension to eps=0.1 |
| compression.inp | Single element, uniaxial compression |
| shear.inp | Single element, simple shear |
| cyclic.inp | Single element, cyclic tension-compression-zero |
| 	ension_8el.inp | 2x2x2 cube, 8 elements, uniaxial tension |
| 	ension_125el.inp | 5x5x5 cube, 125 elements, uniaxial tension |

**Mesh generators:** gen_1000el.py, gen_125el.py — Python scripts
(run with baqus python) that generate structured cube meshes with
elements up to 1000. Edit 
 = 10 to change size.

**ODB extraction:** extract_125el.py — example script for extracting
stress fields from Abaqus ODB.

**Important:** run multi-element tests with cpus=1 until the
static CPUMAT singleton is removed (see DEBUG_CONTEXT.md, point 4).


