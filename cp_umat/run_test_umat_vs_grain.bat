@echo off
rem Regression: UMAT path must reproduce standalone Grain exactly.
rem Usage: run_test_umat_vs_grain.bat  (from anywhere)
setlocal
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "%~dp0"
if not exist build mkdir build

cl /GR- /EHsc- /std:c++17 /O2 /DNDEBUG /DABAQUS_UMAT /I ".." /I "D:\eigen-3.4.0" cp_umat.cpp /c /Fo:build\cp_umat_test.obj
if errorlevel 1 exit /b 1

cl /EHsc /std:c++17 /O2 /DABAQUS_UMAT /I ".." /I "D:\eigen-3.4.0" test_umat_vs_grain.cpp build\cp_umat_test.obj /Fe:build\test_umat_vs_grain.exe
if errorlevel 1 exit /b 1

cd /d "%~dp0\.."
cp_umat\build\test_umat_vs_grain.exe grain_implicit.txt
exit /b %ERRORLEVEL%
