# CP-UMAT Debug Context — RESOLVED

## Root Cause: `*USER MATERIAL` line-split

Abaqus reads at most **8 values per line** in `*USER MATERIAL` data.
All 19 PROPS were on one line, so PROPS[8] (strain_rate_infl) shifted
to 0.0 → gamma_dot = 0 → pure elastic.

**Fix:** split into 3 lines, 8 values each:
```
*USER MATERIAL, CONSTANTS=19
0.0, 106750.0, 60410.0, 28340.0, 2.0, 2.0, 1.0, 6.0
0.588235, 50.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0
115.0, 34.0, 2.25
```

**Result:** S33 = 49.78 MPa (plastic), shifts > 0, tc grows with hardening.
200 increments, 0 cutbacks. 

## Also fixed

2. **load_state on first call** — added `else` clause: init → save_state,
   skip load_state on KINC=1. Otherwise load_state zeroed the freshly
   initialized state.
3. **JSTEP** added to UMAT signature — KINC was reading JSTEP(1) instead
   of the increment number.

## Standalone vs Abaqus

Abaqus rate differs from standalone (0.1 vs 0.0017 s⁻¹ at dt=0.005).
Stress should match when run at same rate. Use compare_abaqus.cpp for
verification.

## CRT Linking Solution (unchanged)

Three changes fixed CRT linking errors:

1. **Grain.hpp** — guarded problematic includes:
   - `#ifndef ABAQUS_UMAT` / `#endif` around `<iostream>`, `<fstream>`, `<stdexcept>`
   - `#ifdef ABAQUS_UMAT` — `params` becomes raw `const MaterialParams*` instead of `shared_ptr`
   - Throw in `solve_implicit` fallback guarded: `return` instead of `throw`

2. **build_cp_umat.bat** — flags: `/GR- /EHsc- /DABAQUS_UMAT`

3. **cp_umat.cpp** — raw `new MaterialParams()` instead of `shared_ptr`
