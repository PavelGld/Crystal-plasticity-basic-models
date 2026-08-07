# CP-UMAT Integration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate crystal plasticity model as Abaqus/Standard UMAT (single grain per integration point, v1).

**Architecture:** Header-only model library in `cp_model/`, UMAT interface layer in `cp_umat/`, validation .inp files in `cp_umat_validation/`.

**Tech Stack:** C++17, Eigen 3.4 (`D:\eigen-3.4.0`), MSVC `cl`, Abaqus/Standard.

## Global Constraints

- Ветка: `implicit-integration` (текущая, поверх неявной схемы).
- Сборка тестов: `cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp` (ОДИН .cpp файл после рефакторинга).
- Сборка UMAT: `cl /EHsc /std:c++17 /O2 /openmp /I ..\..\ /I ..\..\cp_model /I D:\eigen-3.4.0 cp_umat.cpp /c /Fo:cp_umat.obj`
- Спека: `docs/superpowers/specs/2026-07-26-cp-umat-integration-design.md`
- No cout/cin in UMAT code.
- Existing 16 tests must stay green after Task A.
- RVE.h/.cpp kept as-is during refactor (merge into RVE.hpp + delete old).

---

### Task A: Header-only рефакторинг → cp_model/

**Files:** CREATE cp_model/tensor_utils.hpp (rename), cp_model/MaterialParams.hpp (merge h+cpp), cp_model/Grain.hpp (merge h+cpp), cp_model/RVE.hpp (merge h+cpp). MODIFY tests.cpp, EigenProject.cpp (#include paths). DELETE tensor_utils.h, MaterialParams.h/.cpp, Grain.h/.cpp, RVE.h/.cpp.

**Steps:**
- [ ] Create cp_model/ directory
- [ ] Copy & rename tensor_utils.h → cp_model/tensor_utils.hpp (no changes needed — already header-only)
- [ ] Merge MaterialParams.h/.cpp → cp_model/MaterialParams.hpp (copy cpp body into hpp, wrapper guard, remove `#include "MaterialParams.h"` from cpp body, all methods inline)
- [ ] Merge Grain.h/.cpp → cp_model/Grain.hpp (same pattern; keep `extern` includes like `<unsupported/Eigen/MatrixFunctions>`)
- [ ] Merge RVE.h/.cpp → cp_model/RVE.hpp (same)
- [ ] Update tests.cpp: `#include "Grain.h"` → `#include "cp_model/Grain.hpp"`; same for RVE, MaterialParams, tensor_utils
- [ ] Update EigenProject.cpp: same include changes
- [ ] Delete old files: tensor_utils.h, MaterialParams.h/.cpp, Grain.h/.cpp, RVE.h/.cpp
- [ ] Build: `cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp`
- [ ] Run: `tests.exe` → 16/16 PASS (правки только include-путей, без изменения поведения)
- [ ] Commit: `git add cp_model/ tests.cpp EigenProject.cpp && git rm old_files && git commit -m "refactor: header-only model library in cp_model/"`

### Task B: Abaqus-интерфейс (cp_umat/)

**Files:** CREATE cp_umat/cp_umat_types.hpp, cp_umat/cp_umat_interface.hpp, cp_umat/cp_umat.cpp, cp_umat/build_cp_umat.bat.

**cp_umat_types.hpp** — static constexpr 6×6 матрицы:

```cpp
#pragma once
#include <Eigen/Dense>
// Abaqus UMAT ordering: [11,22,33,12,13,23] — tensor components
// Our Voigt: [00,11,22,23,13,12] — engineering (shears ×2)

inline const Eigen::Matrix<double,6,6>& A_to_V() {  // Abaqus Δε → Voigt Δε_eng
    static const Eigen::Matrix<double,6,6> M = (
        Eigen::Matrix<double,6,6>() <<
        1,0,0,0,0,0,    // Voigt[0]=ε₁₁ ← Abaqus[0]=ε₁₁
        0,1,0,0,0,0,    // Voigt[1]=ε₂₂ ← Abaqus[1]=ε₂₂
        0,0,1,0,0,0,    // Voigt[2]=ε₃₃ ← Abaqus[2]=ε₃₃
        0,0,0,0,0,2,    // Voigt[3]=2ε₂₃ ← Abaqus[5]=ε₂₃ (×2)
        0,0,0,0,2,0,    // Voigt[4]=2ε₁₃ ← Abaqus[4]=ε₁₃ (×2)
        0,0,0,2,0,0     // Voigt[5]=2ε₁₂ ← Abaqus[3]=ε₁₂ (×2)
    ).finished();
    return M;
}
inline const Eigen::Matrix<double,6,6>& V_to_A() {  // Voigt Δσ_eng → Abaqus Δσ
    static const Eigen::Matrix<double,6,6> M = (
        Eigen::Matrix<double,6,6>() <<
        1,0,0,0,0,0,    // Abaqus[0]=σ₁₁ ← Voigt[0]
        0,1,0,0,0,0,    // Abaqus[1]=σ₂₂ ← Voigt[1]
        0,0,1,0,0,0,    // Abaqus[2]=σ₃₃ ← Voigt[2]
        0,0,0,0,0,0.5,  // Abaqus[3]=σ₁₂ ← Voigt[5]/2
        0,0,0,0,0.5,0,  // Abaqus[4]=σ₁₃ ← Voigt[4]/2
        0,0,0,0.5,0,0   // Abaqus[5]=σ₂₃ ← Voigt[3]/2
    ).finished();
    return M;
}
```

**cp_umat_interface.hpp** — CPUMAT class wrapping Grain:

Key method signatures:
```cpp
class CPUMAT {
    static constexpr int NSTATV = 70;
    MaterialParamsPtr params_;
    Grain grain_;
    bool initialized_ = false;
public:
    void init_from_props(const double* PROPS, int NPROPS);  // builds params_ + grain_ from PROPS, sets FCC systems
    void load_state(const double* STATEV);    // STATEV → grain_ members
    void save_state(double* STATEV) const;    // grain_ members → STATEV
    void step(const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT);
};
```

init_from_props must:
- Read c11,c12,c44 → construct P_v
- Read hard_mode,rot_mode,integ_mode
- Build FCC {111}⟨110⟩ 12 slip systems (doubled → 24) with M, Wpm, Ms_v, MsE_v
- Read tc0, strain_rate_infl, n_power
- Read KM/iso hardening params based on hard_mode
- Read initial quaternion Q(w,x,y,z)
- Construct MaterialParams with all the above (not from file — from hardcoded FCC + PROPS values)
- Initialize grain_ with params_ and quaternion

The hardcoded FCC slip systems (12 systems, ±b → 24):
```
// {111}<110> planes and directions
const double planes[12][3] = {{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1}};
const double dirs[12][3] = {{0,1,-1},{0,1,1},{0,1,1},{0,1,1},{1,0,-1},{1,0,1},{-1,0,1},{1,0,1},{1,-1,0},{1,1,0},{1,1,0},{-1,1,0}};
```

Hmm, actually the way MaterialParams is constructed from a file in the current code, it reads b1,b2,b3,n1,n2,n3 for each system. For the UMAT we need to hardcode these or compute them from the standard FCC systems. The {111}<110> systems are:
- slip planes: (111), (111), (111), (111) — all 4 {111} variants
- slip directions: <110> — 3 per plane: <110>, <101>, <011> (and their negatives are handled by the ±b doubling)

The standard FCC slip systems (12, numbered for clarity):
| Plane | Direction |
|-------|-----------|
| (111) | [0 1 -1] |
| (111) | [1 0 -1] |
| (111) | [1 -1 0] |
| (11-1) | [0 1 1] |
| (11-1) | [1 0 1] |
| (11-1) | [1 -1 0] |
| (1-11) | [0 1 1] |
| (1-11) | [1 0 1] |
| (1-11) | [-1 1 0] |
| (-111) | [0 1 1] |
| (-111) | [1 0 -1] |
| (-111) | [1 1 0] |

Wait, I need to list them carefully. Let me count:
4 slip planes {111}: (111), (11-1), (1-11), (-111)
3 slip directions <110> per plane: [01-1], [10-1], [1-10] on (111) etc.

Actually, let me just check what the current grain.txt has:
```
0 1 1  1 -1 1    → b=011, n=1-11 → system (1-11)[011]
0 1 1  1 1 -1    → b=011, n=11-1 → system (11-1)[011]
0 -1 1 1 1 1     → b=0-11, n=111 → system (111)[0-11]
0 -1 1 -1 1 1    → b=0-11, n=-111 → system (-111)[0-11]
-1 1 0 1 1 1     → b=-110, n=111 → system (111)[-110]
-1 1 0 1 1 -1    → b=-110, n=11-1 → system (11-1)[-110]
1 1 0 -1 1 1     → b=110, n=-111 → system (-111)[110]
1 1 0 1 -1 1     → b=110, n=1-11 → system (1-11)[110]
1 0 1 1 1 -1     → b=101, n=11-1 → system (11-1)[101]
1 0 1 -1 1 1     → b=101, n=-111 → system (-111)[101]
-1 0 1 1 1 1     → b=-101, n=111 → system (111)[-101]
-1 0 1 1 -1 1    → b=-101, n=1-11 → system (1-11)[-101]
```

These are the standard 12 FCC {111}<110>. I need to hardcode them in init_from_props. I'll provide the exact b,n arrays from grain.txt.

For the cp_umat.cpp entry point:
```cpp
extern "C" void umat_(double* STRESS, double* STATEV, double* DDSDDE,
    double& SSE, double& SPD, double& SCD, double& RPL,
    double* DDSDDT, double* DRPLDE, double* DRPLDT,
    double* STRAN, double* DSTRAN, double* TIME, double* DTIME,
    double& TEMP, double& DTEMP, double* PREDEF, double* DPRED,
    char* CMNAME, int& NDI, int& NSHR, int& NTENS, int& NSTATV,
    double* PROPS, int& NPROPS, double* COORDS, double* DROT,
    double& PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
    int& NOEL, int& NPT, int& LAYER, int& KSPT, int& KINC);

extern "C" void umat(...) { umat_(...); }
```

The CPUMAT instance should be a file-scope static (one per simulation, since UMAT calls are serial in Standard). On KINC=1 with STATEV[0]!=123456.0 → init_from_props. On every call: load_state, step, save_state.

For the build script:
```bat
@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cl /EHsc /std:c++17 /O2 /openmp /I "..\..\" /I "..\..\cp_model" /I "D:\eigen-3.4.0" /D_ENABLE_EXTENDED_ALIGNED_STORAGE cp_umat.cpp /c /Fo:cp_umat.obj
```

- [ ] Create all 4 files with content as specified
- [ ] Compile: `cd cp_umat && build_cp_umat.bat` → cp_umat.obj created without errors
- [ ] Commit: `git add cp_umat/ && git commit -m "feat: CP-UMAT Abaqus interface (cp_umat/)"`

### Task C: Валидация (.inp файлы + README)

**Files:** CREATE cp_umat_validation/tension.inp, compression.inp, shear.inp, cyclic.inp, cp_umat/README.md.

Each .inp: single C3D8 element, 1×1×1 mm, Abaqus/Standard, NLGEOM=NO, Static step, 100-200 increments.

PROPS example (FCC Cu, KM hardening, integ_mode=1, initial orientation = identity quaternion):
```
*USER MATERIAL, CONSTANTS=19
19, 106750.0e6, 60410.0e6, 28340.0e6, 2, 2, 1, 6.0e6, 0.588235, 50.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 115.0e6, 34.0e6, 2.25
```
(PROPS: NPROPS, c11,c12,c44, hard_mode=2, rot_mode=2, integ_mode=1, tc0=6e6, srf=0.588, np=50, Q(1,0,0,0), A=1, q=2, h0=115e6, ts=34e6, alfa=2.25)

Tension .inp: fix bottom face U3=0, prescribe U3 on top face (+0.1 mm displacement).
Compression: same but -0.1 mm.
Shear: top face U1 displacement.
Cyclic: U3 +0.05 → -0.05 → 0.0 displacement.

READMEMD.md: instructions for building and running with Abaqus.

- [ ] Write all 4 .inp files
- [ ] Write cp_umat/README.md (build step, run step, PROPS reference, STATEV layout)
- [ ] Commit: `git add cp_umat_validation/ cp_umat/README.md && git commit -m "docs: CP-UMAT validation files and README"`
