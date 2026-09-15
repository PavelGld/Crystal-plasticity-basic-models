# Example 01 — single element tension

Smoke test: one `C3D8` element, one material, uniaxial tension to 10% strain.
Fastest way to verify that the UMAT builds, links and runs.

## Material
FCC copper, 12 slip systems {111}⟨110⟩, Hutchinson `n=50`, Kocks-Mecking
hardening (`h0=115, ts=34, alfa=2.25`), implicit integration, identity
orientation. See `tension.inp` PROPS block.

## Run
```bash
# 1. build the UMAT (once)
cd ../../cp_umat
build_cp_umat.bat

# 2. run
cd ../examples/01_single_element
run.bat
```
Expected: `tension.odb`, `tension.sta`, `tension.msg`.

## Check
```
abaqus python ../../scripts/extract_s33.py tension s33.txt
```
or in CAE: `abaqus cae database=tension.odb`.

## Notes
- `*UNSYMM=YES` is required (non-symmetric tangent).
- `*CONTROLS, PARAMETERS=LINE SEARCH` is important — Abaqus has it **off** by
  default and convergence suffers without it.
