@ECHO OFF

:: Find the Java binary.
set JAVA=java

:: Run Ecotype Simulation.
"%JAVA%" -jar "%~dp0ecosim.jar" %*
