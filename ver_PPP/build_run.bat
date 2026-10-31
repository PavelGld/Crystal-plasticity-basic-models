@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "D:\Pet_Projects\EigenProject"
cl /EHsc /O2 /openmp /std:c++17 /I . /I "D:\eigen-3.4.0" /D_ENABLE_EXTENDED_ALIGNED_STORAGE ver_PPP\run_traj.cpp /Fe:ver_PPP\run_traj.exe
exit /b %ERRORLEVEL%
