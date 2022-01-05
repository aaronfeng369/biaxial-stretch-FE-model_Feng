@echo off
rem Batch file to run abaqus simulation and python script postprocessing
rem 
rem Record of Revision
rem
rem May-15-2014===YF===Orginal code based on previous versions
rem                    Use umat Struct2D_YF3.for to do parallel computing

::The directory of ABAQUS input file
rem SET INP_DIR=C:\Users\ices\Box Sync\BHV\biaxial testing\Abaqus simulation\Contact assembly\5 sample7mm_preconditioned\21 FxFy with K20MPa for paper\Fx300Fy1500\inp files

:: The input file without the extention .inp
SET INP_ABA=BiaxAssemStruct


:: Stores the current directory for use by POPD command,
:: then change to the specified directory
PUSHD ".\"
rem ABAQUS job=%INP_ABA% >%INP_ABA%.txt interactive


rem the parenthesis is needed for sequential execution of commands
(
ABAQUS user=Struct2D_YF3 job=%INP_ABA% >%INP_ABA%.txt interactive

ABAQUS python BiaxAssemStruct_odb_v8short.py

rem Creat "mat files" folder and move all the postprocessed .mat files to the folder
MKDIR ..\"mat files"
MOVE  .\*.mat ..\"mat files"\
)
