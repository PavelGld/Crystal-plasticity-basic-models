# CP-UMAT Integration — Work Summary & Open Issues

Дата: 2026-08-12
Ветка: `master`

## 1. Overview

Интеграция модели кристаллической пластичности (FCC, Hutchinson power-law, implicit Newton-Raphson) как Abaqus/Standard UMAT. Работа включала рефакторинг в header-only библиотеку, реализацию UMAT interface, исправление CRT-linking проблем и диагностику multi-element convergence failures.

**Multi-element UMAT полностью работает:** 8 и 343 элемента сходятся с текстурным разбросом напряжений.

## 2. Architecture

```
cp_model/                     ← header-only библиотека
├── tensor_utils.hpp          ← Voigt ↔ tensor conversion, Vector6d/Matrix6d types
├── MaterialParams.hpp        ← структура материала (slip systems, elasticity, hardening)
├── Grain.hpp                 ← single crystal integration (explicit + implicit schemes)
└── RVE.hpp                   ← polycrystal homogenization (Taylor upper bound)

cp_umat/                      ← UMAT interface layer
├── cp_umat_types.hpp         ← Abaqus↔Voigt conversion matrices (A_to_V, V_to_A)
├── cp_umat_interface.hpp     ← CPUMAT class declaration
├── cp_umat.cpp               ← CPUMAT implementation + extern "C" UMAT entry
└── build_cp_umat.bat         ← MSVC build → cp_umat.obj

cp_umat_validation/           ← Abaqus .inp files + Python diagnostics
├── tension.inp               ← single C3D8, identity orientation
├── tension_random.inp        ← single C3D8, random orientation
├── tension_8el.inp           ← 2×2×2 cube, 8 elements, identity
├── tension_8el_random.inp    ← 2×2×2 cube, 8 random orientations (FAILING)
└── *.py                      ← ODB extraction scripts
```

## 3. What Works

### 3.1 Standalone tests (`tests.cpp`)
- 14 тестов, 240 assertions — all PASS
- Покрытие: rotation, RVE basic, hard loading, isotropic hardening, implicit self-convergence, tangent FD verification

### 3.2 Single-element UMAT
| Test | Orientation | Result |
|------|-------------|--------|
| tension.inp | identity | ✅ 49.78 MPa, 0 errors |
| tension_random.inp | random | ✅ 28.09 MPa, 0 errors |

### 3.3 Multi-element identity UMAT
| Test | Elements | Result |
|------|----------|--------|
| tension_8el.inp | 8 (2×2×2) | ✅ 49.78 MPa uniform, 0 errors |
| tension_125el.inp | 125 (5×5×5) | ✅ works |
| tension_1000el.inp | 1000 (10×10×10) | ✅ works |

### 3.4 Mini-FEM validation (`fixes/fe_test.cpp`)
- 2×2×2 C3D8 mesh, per-element random orientations
- Global Newton with UMAT calls (no Abaqus)
- **Result: CONVERGES** — 34 increments, 3 cutbacks
- S33 per element: 12-15 MPa (expected scatter from texture)

## 4. What Doesn't Work

### 4.1 Multi-element random orientations in Abaqus
- `tension_8el_random.inp`: **FAILS** at increment 1 with `TOO MANY ATTEMPTS`
- All 64 IPs produce identical stress (9.66 MPa) despite different PROPS quaternions
- Analysis reaches only 0.0634 of the step in 200 increments

## 5. Root Cause Analysis

### 5.1 CRT Linking (resolved)
MSVC-compiled `.obj` referenced `__vcrt_initialize`, `__acrt_initialize` etc. that the Intel/Abaqus linker couldn't resolve.

**Fix:** Eliminate CRT-dependent features from UMAT path:
- Guard `<iostream>`, `<fstream>`, `<stdexcept>` with `#ifndef ABAQUS_UMAT`
- Use raw pointers instead of `shared_ptr`/`unique_ptr`
- Use C arrays instead of `std::map`/`std::vector` in hot paths
- Compile with `/GR- /EHsc- /DABAQUS_UMAT`

### 5.2 Aliasing bug (resolved, from colleague review)
In `CPUMAT::step()`:
```cpp
sv.noalias() = V_to_A() * sv;  // BUG: noalias with self-assignment corrupts data
```
When `noalias()` is used with self-assignment, Eigen writes the result directly into `sv` as the multiplication proceeds. Row 5 of `V_to_A` reads `sv(3)` which was already overwritten by row 3.

**Fix:**
```cpp
sv = (V_to_A() * sv).eval();  // force temporary, then assign
```

**Why it was invisible:**
- Only affects the UMAT wrapper (not standalone Grain/RVE)
- With identity orientation: σ12 = σ23 = 0, so the corruption is zero
- All identity test cases pass; only non-identity multi-element cases fail

### 5.3 Chain tangent across substeps (resolved, from colleague review)
Original: `tangent_buf` = tangent of last substep only.
Fixed: accumulate via chain rule: `tangent_buf = J⁻¹(tangent_buf + weight·P_eff)`

### 5.4 Multi-element convergence in Abaqus — RESOLVED

**Root cause:** Eigen aliasing bug in `CPUMAT::step()`:
```cpp
sv.noalias() = V_to_A() * sv;  // noalias with self-assignment corrupts data
```
When `noalias()` is used with self-assignment, Eigen writes the result directly into `sv` as the multiplication proceeds. Row 5 of `V_to_A` reads `sv(3)` which was already overwritten by row 3. Every UMAT call returned `STRESS[σ₂₃] = σ₁₂/2`, and the true σ₂₃ was lost.

**Fix:** `sv = (V_to_A() * sv).eval()` — forces a temporary.

**Why it was invisible:**
- Only affects the UMAT wrapper (not standalone Grain/RVE)
- With identity orientation: σ₁₂ = σ₂₃ = 0, so the corruption is zero
- All identity test cases pass; only non-identity multi-element cases fail

**Resolution verified:**
- Single-element random: ✅ 28.09 MPa, 0 errors
- Multi-element identity (8-1000 el): ✅ 49.78 MPa uniform, 0 errors
- Multi-element random 8 el: ✅ converges, S33 scatter 12-15 MPa
- Multi-element random 343 el (7×7×7): ✅ **3 часа CPU, 0 errors**, S33 per element 27-74 MPa (102.5% variation — texture effect)

**Files:** `cp_umat/cp_umat.cpp` (aliasing fix, chain tangent, DS_CRITICAL_IMPL=2e-3, PNEWDT, substep cap)

## 6. Implementation Details

### 6.1 Material mode: hard_mode=2 (Kocks-Mecking)
- State variables: stress (6), orientation quaternion (4), critical stress tc[24], accumulated slip shifts[24], time, accumulated deformation, guard flag
- Total: 70 STATEV slots (NSTATV=70)
- Backward Euler integration with Newton-Raphson (tol=1e-8, max 20 iters)
- Line search Armijo (max 5 backtracks)
- Fallback: step halving (max depth 8)
- Substep cap: 200 (returns failure to Abaqus if exceeded)

### 6.2 UMAT state machine
```
if (STATEV[69] != 123456.0):   # first call for this IP
    q = Quaternion(PROPS[10:14])
    grain.set_orientation(q)
    grain.reset_state()
    save_state(grain, STATEV)
    STATEV[69] = 123456.0
else:
    load_state(grain, STATEV)

step(grain, DSTRAN, DTIME, STRESS, DDSDDE, PNEWDT)
save_state(grain, STATEV)
```

### 6.3 Key constants
```cpp
DS_CRITICAL_IMPL = 2e-3    # substepping threshold (raised from 1e-4)
NEWTON_TOL = 1e-8
NEWTON_MAX_ITER = 20
RATIO_MAX = 2.0            # clamp τ/τc in Hutchinson law
MAX_FALLBACK_DEPTH = 8
MAX_SUBSTEPS = 200
```

## 7. Files Changed (from colleague review)

### cp_model/Grain.hpp
- `computeVoigtRotationInPlace(Matrix6d& T, const Eigen::Matrix3d& R)` — takes R as parameter
- `rotate_stiffness_lab()` — uses original T·C·T^T formula
- `solve_implicit(dt, depth, weight)` — chain tangent accumulation
- `tangent_buf = lu.solve((tangent_buf + weight * P_corrected).eval())`
- `step_failed` flag + `get_step_failed()` getter
- `update_H_KM_if_needed()` — public, only for hard_mode==2
- `DS_CRITICAL_IMPL = 2e-3`
- Substep cap: `if (n_substeps > 200) { step_failed = true; return; }`

### cp_umat/cp_umat.cpp
- Aliasing fix: `sv = (V_to_A() * sv).eval()`
- `if (g->get_step_failed()) PNEWDT = 0.25`
- `g->update_H_KM_if_needed()` instead of conditional call
- Diagnostic: STATEV[60-63] = orientation quaternion, STATEV[64] = element number

## 8. Remaining Work

1. **Debug multi-element in Abaqus:** Add SDV output to ODB and verify orientation is applied per-IP
2. **Potentially:** Switch to per-IP Grain objects (eliminate singleton)
3. **Potentially:** Add `*INITIAL CONDITIONS` or verify STATEV initialization
4. **Validation:** Run full `tension_8el_random.inp` and verify texture evolution
5. **Performance:** Enable `cpus>1` after resolving thread safety of static CPUMAT

## 9. Key Learnings

1. **Eigen `noalias()` is dangerous with self-assignment** — always use `.eval()` when the RHS contains the LHS variable
2. **Mini-FEM tests catch logic bugs but not Abaqus integration bugs** — the aliasing bug was invisible to both standalone and mini-FEM tests
3. **CRT linking for Abaqus is fragile** — even `<iostream>` in headers pulled in through the UMAT chain causes link failures
4. **The `T·C·T^T` Voigt rotation formula is only approximate** — it works for identity but accumulates error for non-identity orientations; however, "correct" rotation formulas gave WORSE results (2.56e11 MPa error), suggesting the approximate formula happens to be "least wrong" for the global Newton
