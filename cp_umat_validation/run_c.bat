@echo off
cd /d D:\Pet_Projects\EigenProject\cp_umat_validation
call "C:\Program Files (x86)\Intel\OneAPI\setvars.bat" >nul 2>&1
abaqus job=neper_c user=..\cp_umat\cp_umat.obj cpus=4 mp_mode=threads standard_parallel=all int
