@echo off
REM Run 03_neper_200grains_symmetry (deck already generated as model.inp)
cd /d "%~dp0"
call "C:\Program Files (x86)\Intel\OneAPI\setvars.bat" >nul 2>&1
abaqus job=model user=..\..\cp_umat\cp_umat.obj cpus=4 mp_mode=threads standard_parallel=all int
