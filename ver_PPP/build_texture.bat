@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "D:\Pet_Projects\EigenProject"
if not exist ver_PPP\build mkdir ver_PPP\build

cl /GR- /EHsc- /std:c++17 /O2 /DNDEBUG /DABAQUS_UMAT /I . /I "D:\eigen-3.4.0" cp_umat\cp_umat.cpp /c /Fo:ver_PPP\build\cp_umat_tex.obj
if errorlevel 1 exit /b 1

cl /EHsc /O2 /openmp /std:c++17 /DABAQUS_UMAT /I . /I "D:\eigen-3.4.0" ver_PPP\run_umat_texture.cpp ver_PPP\build\cp_umat_tex.obj /Fe:ver_PPP\run_umat_texture.exe
if errorlevel 1 exit /b 1

cl /EHsc /O2 /openmp /std:c++17 /I . /I "D:\eigen-3.4.0" ver_PPP\export_initial.cpp /Fe:ver_PPP\export_initial.exe
if errorlevel 1 exit /b 1
echo BUILD_OK
