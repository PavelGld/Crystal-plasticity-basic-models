@echo off
cd /d "%~dp0"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cl /GR- /EHsc- /std:c++17 /O2 /DNDEBUG /DABAQUS_UMAT /I ".." /I "D:\eigen-3.4.0" cp_umat.cpp /c /Fo:cp_umat.obj
if %ERRORLEVEL% EQU 0 echo SUCCESS: cp_umat.obj built
