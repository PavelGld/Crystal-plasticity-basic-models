# CP-UMAT — Crystal Plasticity UMAT for Abaqus/Standard

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

`cp_umat_validation/` contains 4 single-element test cases:
- `tension.inp` — uniaxial tension to ε = 0.1
- `compression.inp` — uniaxial compression to ε = 0.1
- `shear.inp` — simple shear to γ = 0.1
- `cyclic.inp` — cyclic tension→compression→zero, 5% strain amplitude

## Known limitations (v1)

- NLGEOM=NO only (small strain); no DROT/DFGRD handling
- Single grain per integration point (no RVE averaging)
- FCC only (slip systems hardcoded)
- No temperature/rate dependence beyond the power law
