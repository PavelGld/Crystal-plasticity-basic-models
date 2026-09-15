@echo off
REM ============================================================================
REM  Run example 01 -- single element uniaxial tension
REM  Builds nothing; expects ..\..\cp_umat\cp_umat.obj to exist
REM  (build it first: cd cp_umat && build_cp_umat.bat)
REM ============================================================================
cd /d "%~dp0"
call "C:\Program Files (x86)\Intel\OneAPI\setvars.bat" >nul 2>&1
abaqus job=tension user=..\..\cp_umat\cp_umat.obj cpus=4 mp_mode=threads standard_parallel=all int
