@echo off
rem Interface tests: ordering, objectivity, energy, large strain.
setlocal
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "%~dp0"
if not exist build mkdir build

cl /GR- /EHsc- /std:c++17 /O2 /DNDEBUG /DABAQUS_UMAT /I ".." /I "D:\eigen-3.4.0" cp_umat.cpp /c /Fo:build\cp_umat_iface.obj
if errorlevel 1 exit /b 1

cl /EHsc /std:c++17 /O2 /DABAQUS_UMAT /I ".." /I "D:\eigen-3.4.0" test_umat_interface.cpp build\cp_umat_iface.obj /Fe:build\test_umat_interface.exe
if errorlevel 1 exit /b 1

cd /d "%~dp0\.."
cp_umat\build\test_umat_interface.exe grain_implicit.txt
exit /b %ERRORLEVEL%
