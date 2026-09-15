# CP-UMAT Debug Context — RESOLVED + Multi-Element

## Root Cause (resolved): `*USER MATERIAL` line-split

Abaqus reads at most **8 values per line** in `*USER MATERIAL` data.
All 19 PROPS were on one line, so PROPS[8] (strain_rate_infl) shifted
to 0.0 → gamma_dot = 0 → pure elastic.

**Fix:** split into 3 lines, ≤8 values each. See tension.inp for correct format.

## Also fixed

2. **load_state on first call** — added `else` clause: init → save_state,
   skip load_state on KINC=1 (was zeroing freshly initialized state).
3. **JSTEP** added to UMAT signature — was missing, KINC read garbage.
4. **Anti-rotation BC** — added second pinned node (U2=0 at corner) to
   prevent rigid-body rotation about Z.

## Result: single-element

S33 = 49.78 MPa at eps=0.1. Shifts > 0, tc grows (hardening works).
200 increments, 0 cutbacks. Abaqus vs standalone: stress matches within
rate sensitivity (strain rate differs: 0.1 vs 0.0017 s⁻¹).

## Multi-element tests

| File | Elements | Time | Result |
|------|----------|------|--------|
| tension_8el.inp  | 8 (2×2×2)  | ~2 min  | OK, cpus=1 only |
| tension_125el.inp | 125 (5×5×5) | ~22 min | OK, cpus=1 only |

ODB in project root: `tension_125el.odb`, `tension_8el.odb`.

**cpus>1 crashes** — `static CPUMAT` is a singleton shared across threads.
Fix (p.4 below) needed before parallel runs.

## Remaining issues (from review)

### 1. `static CPUMAT` — not thread-safe (p.4)

Single `Grain` instance shared across all integration points.
- `cpus=1` works (serial)
- `cpus=4` → segfault (race on Stress_v, tc, gamma_dot)
- Fix: create Grain per-call from STATEV, OR store state per-IP in a map

### 2. `A_to_V` ×2 for shears (p.5)

`DSTRAN` in Abaqus for shears is engineering (γ=2ε), same as our Voigt.
`A_to_V` should be pure permutation (no ×2). `V_to_A` for stress is correct.
Also `P_v(3,3)=c44` should be `2*c44` for consistency — the two errors cancel
in elastic regime but `computeDeviatoricDu()` overestimates du for shear loading.

### 3. No PNEWDT signal on fallback (p.6)

`solve_implicit` fallback returns without setting PNEWDT. Abaqus doesn't know
the step failed. Also `get_tangent_lab_into` throws if tangent_valid=false —
replace with elastic fallback under ABAQUS_UMAT.

### 4. Integ_mode validation, unsymmetric tangent, SSE/SPD (p.7)

Minor items — see full review in the original message.

## Stress uniformity note

Multi-element tests show **uniform stress** (49.78 MPa everywhere) because
all elements share identical orientation (identity), material (PROPS), and
loading (kinematic tension). Non-uniform fields require:
- Different initial orientations per element (texture)
- OR geometric imperfection/necking at large strains
- OR non-uniform loading (bending, indentation, etc.)
