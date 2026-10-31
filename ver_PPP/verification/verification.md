# Reference dataset for UMAT verification

This directory contains reference macro stress histories for five strain
trajectories, produced by an independent implementation of the two-level
statistical crystal-plasticity model used in this repository. The data are the
external reference for the Abaqus UMAT. Everything needed to use the dataset is
in this file; the files are self-contained.

## 1. Model formulation

**Macro level.** The macroscopic weighted Kirchhoff stress tensor `K` is the
volume average over all grains:

```
K = <k>
```

**Meso level (single grain).** The co-rotated rate of the weighted Kirchhoff
stress is the elasticity tensor applied to the elastic rate:

```
k^cor := dk/dt + k*w - w*k = n : (L - w - sum_k gdot^(k) b^(k) n^(k))
```

Slip rates (power-law viscosity):

```
gdot^(k) = g0 * (tau^(k) / tauc^(k))^m * H(tau^(k) - tauc^(k)),  k = 1..K
tau^(k)  = b^(k) n^(k) : k
```

**Taylor lattice rotation:**

```
w = 1/2 (L - L^T) - 1/2 sum_k gdot^(k) (b^(k) n^(k) - n^(k) b^(k))
```

**Anisotropic hardening:**

```
gdot_tauc^(k) = sum_j h^(kj) |gdot^(j)|
h^(kl) = [q_lat + (1 - q_lat) delta^(kl)] h^(l)
h^(l)  = h0 * |1 - tau_dis^(l) / tau_sat|^a
tau_dis^(l) = sum_j h^(lj) gamma^(j)
```

Here `k` is the weighted Kirchhoff stress tensor of the grain (essentially the
deviator; its trace stays below `6e-5` GPa), `n` is the elasticity tensor,
`L` is the velocity gradient, `b^(k)` and `n^(k)` are the unit slip direction
and slip-plane normal, and `K = 24` is the doubled number of FCC slip systems.
`tau_dis` starts from the initial critical stress `tau0`, so it coincides with
the critical resolved shear stress used by the UMAT (`tc0` in `PROPS`).

## 2. Material parameters (FCC aluminium alloy)

| Parameter | Symbol | Value |
|---|---|---|
| Elastic constants | c11, c12, c44 | 106.75, 60.41, 28.34 GPa |
| Rate sensitivity | m | 50 |
| Effective reference slip rate | g0 | 1e-3 s^-1 at the imposed strain-rate intensity 1e-3 s^-1 |
| Hardening exponent | a | 2.25 |
| Hardening multiplier | h0 | 0.115 GPa |
| Latent hardening | q_lat | 2 |
| Saturation shear stress | tau_sat | 0.030 GPa |
| Initial critical shear stress | tau0 | 0.006 GPa |

The effective reference slip rate is the value that reproduces the reference
response at the imposed strain-rate intensity of `1e-3 s^-1`. In the UMAT's
Hutchinson form

```
gdot = k * D_eq * (tau / tauc)^n
```

this corresponds to `k = 1.0` with `D_eq = 1e-3 s^-1` (i.e. an effective prefactor
of `1e-3`). The reference flow law uses a threshold; at this imposed rate the
response is effectively rate-independent, which is why the exact prefactor is
weakly observable (changing it by a factor of 100 changes the stress by about
10% at E = 0.2).

## 3. RVE

- 500 grains, FCC lattice, cubic symmetry, one polycrystal;
- equal grain volume fractions (1/500 each);
- initial orientations: uniform random over the cubic fundamental domain, fixed
  seed, identical for all five trajectories;
- the orientations are given in `quaternions_start.txt`; `quaternions_end.txt`
  contains the final orientations of the reference run.

## 4. Numerical scheme

- time step `dt = 5e-3 s`;
- explicit Euler integration of stress, strain, orientation and critical
  stresses within each step;
- Taylor spin; macroscopic stress is the volume-fraction-weighted mean;
- the accumulated strain intensity is

```
E = integral sqrt(2/3 D':D') dt
```

and every trajectory is scaled so that `dE/dt = 1e-3 s^-1`, hence `E = 1e-3 t`
and the run has 200 000 steps.

## 5. Trajectories

Each trajectory is run from a fresh state to `E = 1.0`.

| name | L (velocity gradient) | phases |
|---|---|---|
| `tension` | `1e-3 (k1k1 - 0.5 k2k2 - 0.5 k3k3)` | 1 |
| `compression` | `-1e-3 (k1k1 - 0.5 k2k2 - 0.5 k3k3)` | 1 |
| `simple_shear` | `sqrt(3) * 1e-3 k1k2` | 1 |
| `planar` | `(sqrt(3)/2) * 1e-3 (k1k1 - k2k2)` | 1 |
| `complex` | tension for E = 0..0.5, then simple shear for E = 0.5..1.0 | 2 |

For each `L`, `2/3 D':D' = (1e-3)^2` exactly, so the E grid is analytic:
after `step` steps, `E = 5e-6 * step`.

## 6. Output format

Per trajectory, directory `<name>/`.

### `stress_history.txt`

Macro stress tensor history vs accumulated strain intensity. Tab-separated,
12 significant digits, units **GPa**.

```
E       S11   S22   S33   S12   S13   S23   S_trace
0      ...
0.01   ...
...
1      ...
```

Header plus 101 data rows; `E = 0.01 k`, `k = 0..100`. `S_trace = tr(S)` is a
diagnostic (the material is nearly incompressible; the trace stays below
`6e-5` GPa). The stress is the weighted Kirchhoff stress; for the small elastic
strains here it is essentially the Cauchy stress.

### `quaternions_start.txt` / `quaternions_end.txt`

Grain orientations at the start and at `E = 1.0`. Tab-separated, 12 significant
digits.

```
grain_id   w   x   y   z
0          ...
...
499        ...
```

500 rows; quaternions are unit and sign-ambiguous (`q` and `-q` are the same
rotation). Start files are accurate to about `1e-12`; end files carry the
orientation drift accumulated over 200 000 steps (about `1e-6`).

## 7. Sanity values at E = 1.0 (S11, GPa)

| trajectory | S11 |
|---|---|
| tension | +0.066848179679 |
| compression | -0.0698296403372 |
| simple_shear | +0.0108104818491 |
| planar | +0.0485238270244 |
| complex | +0.00746942772162 |

## 8. How to use the dataset

- Per step the UMAT receives `dt = 5e-3 s` and the constant phase `L` from
  the table; the per-step increment of E is `5e-6`; total steps per trajectory
  is 200 000.
- Compare the macro stress at `E = 0, 0.01, ..., 1.0` (101 checkpoints).
- Initial grain orientations are `quaternions_start.txt`; rebuild rotation
  matrices from the quaternions (sign is irrelevant).
- Optional end-state check: compare `quaternions_end.txt` up to quaternion
  sign (or compare rotation matrices).
- The reference values in these files are the ground truth at this
  discretization. The reference uses explicit Euler; the UMAT uses an implicit
  scheme, so percent-level differences at large strain are expected and are a
  model/discretization difference, not an implementation error.
