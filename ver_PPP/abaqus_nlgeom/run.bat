@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul 2>&1
cd /d "%~dp0"
abaqus job=model input=model user=D:\Pet_Projects\EigenProject\cp_umat\cp_umat.obj cpus=1 interactive
exit /b %ERRORLEVEL%
