# CP-UMAT Debug Context — Elastic Response

## Current State

UMAT links with Abaqus. Analysis (tension.inp) completes: 200 increments, 0 cutbacks, 31s.
**Material behaves elastically** — stress 6300 MPa vs expected 33 MPa at eps=0.1.

## State Variables (from ODB)

- STATEV(70) guard = 123456.0 — init confirmed
- STATEV(11) tc[0] = 6.0 — correctly initialized from PROPS
- STATEV(35)..STATEV(58) shifts[0..23] = all zero after 200 increments — NO slip
- STATEV(1)..STATEV(6) stress_v = [0, 0, 6308.7, 0, 0, 0] MPa — fully elastic

## Comparison

Standalone (same strain rate, dt=0.005, target eps=0.1):
- S33 = 33.2 MPa (plastic, well below yield)
- Shifts accumulate normally

Abaqus UMAT:
- S33 = 6308.7 MPa (elastic-dominated, near c11 * eps)
- Shifts stay zero — no plastic activation

## CRT Linking Solution (working)

Three changes fixed CRT linking errors:

1. **Grain.hpp** — guarded problematic includes:
   - `#ifndef ABAQUS_UMAT` / `#endif` around `<iostream>`, `<fstream>`, `<stdexcept>`
   - `#ifdef ABAQUS_UMAT` — `params` becomes raw `const MaterialParams*` instead of `shared_ptr`
   - File-based constructor guarded out (`#ifndef ABAQUS_UMAT`)
   - Throw in `solve_implicit` fallback guarded: `return` instead of `throw` under ABAQUS_UMAT

2. **build_cp_umat.bat** — flags: `/GR- /EHsc- /DABAQUS_UMAT` (disable RTTI, disable exception handling)

3. **cp_umat.cpp** — raw `new MaterialParams()` to populate params, no `shared_ptr`/`make_shared`

## UMAT State Machine (verify this!)

In `extern "C" umat_()`:

```cpp
static CPUMAT cpumat;

if (KINC == 1 && std::abs(STATEV[69] - 123456.0) > 0.5) {
    cpumat.init_from_props(PROPS, NPROPS);
    cpumat.save_state(STATEV);  // write initial state (tc=6, stress=0)
    STATEV[69] = 123456.0;
} else {
    cpumat.load_state(STATEV);
}

cpumat.step(DSTRAN, *DTIME, STRESS, DDSDDE, PNEWDT);
cpumat.save_state(STATEV);
```

Note: `init_from_props` creates `new Grain(raw_params, q)`. Then on first call,
`save_state` writes initial state, then `step()` runs, then `save_state` overwrites
with step result. On subsequent calls, `load_state` restores state from previous
step, `step()` updates, `save_state` writes back.

## CPUMAT::step — Strain Path

1. `A_to_V() * DSTRAN` — converts Abaqus tensor DELTA_epsilon to Voigt engineering
2. `voigtToTensorInPlace` — Voigt → 3x3 tensor
3. `L = D_tensor / DTIME` — strain rate tensor (pure deformation, no spin — NLGEOM=NO)
4. `grain_->Simulate_Step(DTIME, L)` — crystal plasticity integration
5. `grain_->get_stress_lab()` → `V_to_A()` → `STRESS[]` — lab stress to Abaqus tensor
6. `grain_->get_tangent_lab_into(T_buf_)` → `V_to_A * T * A_to_V` → `DDSDDE[]`

## Hypotheses for Elastic Behavior

### H1: Newton never converges, state lost via fallback return

`solve_implicit` fallback path under ABAQUS_UMAT:
```cpp
if (!converged) {
    if (depth >= MAX_FALLBACK_DEPTH)
        return;  // <-- leaves Stress_v, tc, shifts UNCHANGED
    solve_implicit(0.5 * dt, depth + 1);
    solve_implicit(0.5 * dt, depth + 1);
    return;
}
```

On first plastic step: tau/tc clamped at 2.0, gamma_dot = rf * 2^50 (huge).
Residual R explodes. Newton step must pull sigma down by ~50 MPa.
Line search with Armijo may fail (rnorm_new > (1 - 1e-4*lambda) * rnorm for all lambda).
Accepted=false → break → converged=false → fallback return → state unchanged.

**Check**: add a convergence flag to STATEV. Count iterations. If depth reaches 8,
the throw (now return) happens and step is lost.

### H2: Strain rate wrong — du computed from stale D_v

`computeDeviatoricDu()` uses D_v member set in Simulate_Step from L.
L is set from DSTRAN/DTIME. If D_v ends up zero or tiny, rate_factor=0 → no plasticity.

**Check**: verify du value at first step equals standalone's du (~0.066 at Abaqus rate).

### H3: Voigt/Abaqus conversion mismatch

A_to_V and V_to_A for normal components are identity. For shear: permute + x2 / /2.
If the permute is wrong for DDSDDE, Abaqus global Newton may give wrong DSTRAN.
But STRESS output matches STATEV — internal consistency OK.

**Check**: compare DSTRAN values received by UMAT vs expected from boundary conditions.

### H4: load_state/save_state corrupts state

On first call: init writes tc=6, stress=0 via save_state.
Then step runs... but if step() internally reads STATEV (via load_state), it might
overwrite the initialized grain state. The current code skips load_state on first
call (due to if/else block), so this should NOT happen.

**Check**: verify grain state before/after step on first call matches expectation.

## Key Files

- `cp_model/Grain.hpp` — core model (guarded for UMAT)
- `cp_umat/cp_umat_types.hpp` — A_to_V / V_to_A matrices
- `cp_umat/cp_umat_interface.hpp` — CPUMAT class
- `cp_umat/cp_umat.cpp` — implementation + extern "C" entry
- `cp_umat/build_cp_umat.bat` — MSVC build to cp_umat.obj
- `cp_umat_validation/tension.inp` — single C3D8 uniaxial tension

## Quick Diagnostic

Add to `extern "C" umat_()` after `cpumat.step()` and before `cpumat.save_state()`:

```cpp
Grain* g = cpumat.grain();
STATEV[63] = (g && g->get_shifts(0) > 0) ? 1.0 : 0.0;  // any slip?
STATEV[64] = (g) ? g->get_new_step(0) : -1.0;            // system 0 active?
```

Then check SDV64 in ODB: if always 0 → Newton never converges.
