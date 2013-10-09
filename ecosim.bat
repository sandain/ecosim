@ECHO OFF

:: Find the Java binary.
set JAVA=java

:: Find the absolute path to this script, and change to that directory.
set ABS_PATH=%~d0%~p0
cd /D "%ABS_PATH%"

:: Run Ecotype Simulation.
cmd /c %JAVA% -jar bin\EcoSim.jar %*

