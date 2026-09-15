# Example 02 — 343-element cube with random orientations

A 7×7×7 cube of `C3D8` elements. **Every element is its own grain**: it gets an
individual `*MATERIAL` with a random quaternion in `PROPS[10..13]`. Uniaxial
tension along **Z** to 30% strain (`NLGEOM=NO`).

This is the classic "material-per-element" polycrystal and the reference case
against which the standalone RVE was validated.

## Files
| File | Purpose |
|---|---|
| `gen.py` | generates `model.inp` (seed=42, 343 random quaternions) |
| `model.inp` | the deck (regenerate with `python gen.py`) |
| `run.bat` | builds the Abaqus command |

## Run
```bash
python gen.py          # (optional) regenerate the deck
run.bat
```

## Notes
- `NLGEOM=NO`: small-strain formulation, 30% strain.
- Individual `*ELSET` / `*SOLID SECTION` / `*MATERIAL` per element (343 of each).
- `*UNSYMM=YES` + `*CONTROLS, PARAMETERS=LINE SEARCH` (off by default in Abaqus).
- BCs: bottom face `U3=0`, top face `U3=0.3`, plus 2 anti-rigid-body constraints.
