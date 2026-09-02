# CP-UMAT: Crystal Plasticity UMAT for Abaqus/Standard

Single-grain per integration point. Implicit Newton–Raphson integration with consistent tangent.

**Status (2026-08-12):** Bugs #1–4 from review implemented. Single-element random orientation now converges. Multi-element random orientations still fail — see § Known Issues.
FCC {111}⟨110⟩ slip systems (hardcoded), Hutchinson power-law, Kocks-Mecking or isotropic hardening.

## Quick Start

```bat
cd cp_umat
build_cp_umat.bat                          # MSVC 2022 → cp_umat.obj
cd ..
abaqus job=tension input=cp_umat_validation\tension.inp user=cp_umat\cp_umat.obj cpus=1
```

**Important:** use `cpus=1` — the current UMAT is not thread-safe (single `Grain` singleton — see § Known Issues).

ODB output lands in the project root. To redirect: `abaqus job=X ... output=cp_umat_validation\output\`. Abaqus output files (`.odb`, `.msg`, `.sta`, `.com`, etc.) are git-ignored.

## Files

```
cp_model/                     ← Header-only crystal plasticity library
├── tensor_utils.hpp          Voigt ↔ tensor conversion, Voigt types (Vector6d, Matrix6d)
├── MaterialParams.hpp        Material parameters struct (elasticity, slip systems, hardening)
├── Grain.hpp                 Single crystal integration (explicit + implicit Newton–Raphson)
└── RVE.hpp                   Polycrystal homogenization (Taylor upper bound)

cp_umat/                      ← Abaqus UMAT interface
├── cp_umat_types.hpp         Abaqus ↔ Voigt conversion matrices (A_to_V, V_to_A)
├── cp_umat_interface.hpp     CPUMAT class (PROPS → MaterialParams, STATEV ↔ Grain)
├── cp_umat.cpp               extern "C" umat_ / umat entry + CPUMAT implementation
├── build_cp_umat.bat         MSVC build script → cp_umat.obj
├── DEBUG_CONTEXT.md          Debug log — resolved CRT + unresolved issues
└── README.md                 This file

cp_umat_validation/           ← Abaqus input files + validation scripts
├── tension.inp               Single C3D8 element, uniaxial tension to ε=0.1
├── compression.inp           Single element, uniaxial compression
├── shear.inp                 Single element, simple shear
├── cyclic.inp                Single element, cyclic tension→compression→zero
├── tension_8el.inp           2×2×2 cube, 8 elements, uniaxial tension
├── tension_125el.inp         5×5×5 cube, 125 elements
├── tension_1000el.inp        10×10×10 cube, 1000 elements
├── tension_random.inp        Single element with random initial orientation
├── gen_1000el.py             Structured cube mesh generator (10×10×10)
├── gen_random.py             Mesh generator with random element orientations (multi-element — broken, see § PROPS)
├── gen_random_single.py      Random quaternion .inp generator (single element — works)
├── gen_multi_1el.py          Multi-PROPS format test (1 element, CONSTANTS=24 — broken, see § PROPS)
├── extract_125el.py          ODB stress field extraction script
└── output/                   Abaqus output target directory
```

## PROPS Layout (Single-Element, Working)

`*USER MATERIAL, CONSTANTS=19` (FCC Cu, Kocks-Mecking hardening, implicit integration):

```
*USER MATERIAL, CONSTANTS=19
0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0
0.588235, 50.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0
115.0, 34.0, 2.25
```

| PROPS[i] | Meaning |
|----------|---------|
| 0 | unused (reserved) |
| 1-3 | c11, c12, c44 (MPa — mm-model units) |
| 4 | hard_mode: 0=none, 1=isotropic, 2=Kocks-Mecking |
| 5 | rot_mode: 0=none, 1=Mandel, 2=Taylor |
| 6 | integ_mode: 0=explicit, 1=implicit (use 1) |
| 7 | tc0 — initial CRSS (MPa) |
| 8 | strain_rate_infl — γ̇₀ |
| 9 | n_power — rate sensitivity exponent |
| 10-13 | initial quaternion (qw, qx, qy, qz) |
| 14-18 | (hard_mode=2) A, q, h0, ts, alfa |
| 14-15 | (hard_mode=1) a, b |

NPROPS = 19 (KM), 16 (iso), 14 (none).

**Critical:** Abaqus reads at most **8 values per data line** in `*USER MATERIAL`. All `.inp` files must split PROPS into multiple lines of ≤8 values each. The initial bug was all 19 values on one line → PROPS[8] (strain_rate_infl) shifted to 0 → purely elastic response.

## STATEV Layout (NSTATV = 70)

| STATEV[i] | Content |
|-----------|---------|
| 0-5 | stress_v (Voigt, crystal frame) |
| 6-9 | quaternion (w, x, y, z) |
| 10-33 | tc[0..23] — critical resolved shear stress, 24 systems |
| 34-57 | shifts[0..23] — accumulated slip per system |
| 58-66 | reserved |
| 67 | sim_time |
| 68 | acc_def — accumulated equivalent strain |
| 69 | guard = 123456.0 when initialized |

The guard (STATEV[69]) distinguishes first call per integration point (needs init from PROPS) from subsequent calls (load from STATEV). This is the standard Abaqus pattern for detecting restart vs fresh start — Abaqus preserves STATEV across restarts.

## What Works

### Single-element tests

| Test | PROPS | Quaternion | S33 (MPa) | Time | Notes |
|------|-------|------------|-----------|------|-------|
| tension.inp | 19, identity | (1,0,0,0) | 49.78 | 31s | 200 inc, 0 cutbacks |
| tension.inp | 19, random | (-0.19,-0.23,-0.15,0.94) | 25.50 | 26s | Anisotropy confirmed |

### Multi-element uniform tests (cpus=1)

| Test | Elements | S33 | Time |
|------|----------|-----|------|
| tension_8el.inp | 8 (2×2×2) | 49.78 uniform | ~2 min |
| tension_125el.inp | 125 (5×5×5) | 49.78 uniform | ~22 min |

All elements share identical identity orientation → uniform response. Correct physics.

### CRT Linking

CRT symbols (`__vcrt_initialize`, `__acrt_initialize`, `_purecall`, `_Get_deleter`) prevented Abaqus from linking the UMAT object file. Three changes fixed this:

1. **`Grain.hpp` — guarded includes:** `#ifndef ABAQUS_UMAT` around `<iostream>`, `<fstream>`, `<stdexcept>`. `#ifdef ABAQUS_UMAT` changes `MaterialParamsPtr params` (shared_ptr) to raw `const MaterialParams* params`.
2. **`cp_umat.cpp/build_cp_umat.bat` — compile flags:** `/GR- /EHsc- /DABAQUS_UMAT` disables RTTI, exception handling, and activates UMAT-specific code paths.
3. **No `std::map`, `std::vector`, `std::shared_ptr` in UMAT code path** — all replaced with raw pointers and C arrays.

The EMSC UMAT reference (`UMAT/UMAT/bin/umat.obj`) avoids CRT by using only C headers (`<cstdio>`, `<cstdlib>`, `<cmath>`) and fixed-size Eigen — no iostream, no shared_ptr, no exceptions.

## Multi-Element PROPS — Broken

**Problem:** Abaqus input file processor rejects `*USER MATERIAL` with more than 19 constants.

Attempted: `CONSTANTS=24` (19 base + 1 n_elems + 4 for 1 quaternion), `CONSTANTS=52` (19+1+32 for 8 elements). Both fail with:
```
***ERROR: THERE ARE INVALID DATA ASSOCIATED WITH THIS USER DEFINED MATERIAL DEFINITION
```

This happens during `.inp` parsing (Analysis Input File Processor), before any UMAT code runs. The single-element format (19 PROPS, 3 lines) works perfectly. The error is purely in the Abaqus parser, not in our C++ code.

Attempted fixes (none worked):
- Merging `n_elems` value with previous line (4 values instead of 3+1)
- Python-generated `.inp` (excluded PowerShell formatting issues)
- 24 constants (just 1 element) — same error

**Hypotheses:**
1. Abaqus has a compile-time limit on `*USER MATERIAL` constants — some versions cap at 20-80
2. `*USER MATERIAL` needs a flag like `UNSYMM` for large arrays
3. The `*DEPVAR` line must come before `*USER MATERIAL` (currently after)
4. `*MATERIAL, NAME=CPFCC` name conflict with other sections

**Workaround (not implemented):** Keep single-element PROPS format for material params. Pass orientations via a separate mechanism — e.g., `*INITIAL CONDITIONS, TYPE=SOLUTION` reading grain orientations from an external file (`*INITIAL CONDITIONS, TYPE=SOLUTION, INPUT=orientations.txt`).

## Known Issues

### Thread Safety (cpus>1)

`static CPUMAT cpumat` with a single `Grain* grain_` shared across all integration points. With cpus>1, Abaqus calls UMAT from multiple threads → race condition on `grain_` members → segfault.

Fix options:
1. Per-call `Grain` from STATEV (no persistent object) — fully thread-safe, small allocation cost
2. Thread-local storage (`__declspec(thread)`) for grain pointers
3. Pre-allocate `Grain* grains[MAX_IPS]` and index by `NOEL*8+NPT` — but this requires solving the static-array CRT issue (§ CRT and Static Arrays)

Current workaround: always use `cpus=1`.

### CRT and Static Arrays

Static C arrays of Eigen types (e.g., `Eigen::Quaterniond g_orientations[MAX_ELEMS]`) cause CRT initialization calls at DLL load time. The Abaqus linking environment can't resolve these.

Fix: use raw `double` arrays and construct Eigen objects on demand (not as static globals).

### DDSDDE Symmetry

The consistent tangent `J⁻¹·(P − c·gᵀ)` is slightly non-symmetric due to the `c·gᵀ` term. Abaqus uses only the symmetric part unless `*USER MATERIAL, UNSYMM` is specified. Without UNSYMM, global Newton convergence degrades from quadratic to linear. For single-element tests this is harmless (only 1 global iteration needed). For large models, add `UNSYMM` to the material definition.

### PNEWDT on Fallback

If `solve_implicit` reaches max fallback depth, it returns without updating `PNEWDT` (stays 1.0). Abaqus doesn't know the step failed → may waste iterations. Fix: set `PNEWDT = 0.25` on fallback.

### NLGEOM=NO Only

v1 supports small-strain only. `DROT` is ignored. For large deformations (texture evolution at ε > 0.1), NLGEOM=YES with `DFGRD0`/`DFGRD1` handling is needed.

### A_to_V for Shears

`DSTRAN` in Abaqus for shear components is engineering (γ = 2ε), same as our Voigt. The current `A_to_V` has an extra ×2 on shears. This is compensated by `P_v(3,3) = c44` instead of `2*c44`. The errors cancel in the elastic regime, but `computeDeviatoricDu()` overestimates `du` for shear-dominated loading, affecting plastic kinetics in `shear.inp`.

## CRT Linking Solution (Reference)

When adding new code to the UMAT path, these rules avoid CRT symbol conflicts:

1. **No `<iostream>`, `<fstream>`, `<stdexcept>`** in any file compiled with `/DABAQUS_UMAT`
2. **No `std::shared_ptr`, `std::unique_ptr`** — use raw pointers
3. **No `std::map`, `std::vector`, `std::string`** — use C arrays, `const char*`
4. **No static Eigen objects** — construct on stack or heap, not in static/global storage
5. **Compile with `/GR- /EHsc- /DABAQUS_UMAT`**
6. **Run Abaqus with Intel oneAPI `setvars.bat`** for the linker

The EMSC UMAT reference demonstrates this pattern: `umat.cpp` + pure headers with `<cstdio>`, `<cstdlib>`, fixed-size Eigen only.

## Build System

Tests (standalone, no Abaqus):
```bat
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp
tests.exe
```

UMAT (Abaqus):
```bat
cd cp_umat
build_cp_umat.bat    # uses /GR- /EHsc- /DABAQUS_UMAT
```

Requires: Visual Studio 2022 (Community), Eigen 3.4 at `D:\eigen-3.4.0`, Intel oneAPI (for Abaqus linker, `C:\Program Files (x86)\Intel\oneAPI\setvars.bat`), Abaqus 2022.
