Ecotype Simulation README


INSTALLATION
================================================================================

Extract the contents of the archive to your hard drive:

  Windows: C:\EcoSim\
  POSIX (Linux, OSX, etc): ~/ecosim/

Create a Desktop shortcut to the shell script for your platform.

  Windows: C:\EcoSim\runEcoSim.bat or C:\EcoSim\runEcoSim.vbs
  POSIX (Linux, OSX, etc): ~/ecosim/runEcoSim.sh



REQUIREMENTS - EXECUTION
================================================================================

Java 6 JRE
------------------------------
A Java 6 Runtime Environment (JRE) is required to execute this program.  You
can download it here:

  http://java.sun.com/javase/downloads/index.jsp

NJplot
------------------------------
The NJplot executable must be accessible via your system's PATH, or in the
EcoSim directory.  You can download it here:

  http://pbil.univ-lyon1.fr/software/njplot.html


PHYLIP
------------------------------
The PHYLIP executable must be in your system's PATH, or in the EcoSim
directory.  You can download it here:

  http://www.phylip.com/


REQUIREMENTS - COMPILATION
================================================================================

To compile the Fortran programs, you will need to have a Fortran compiler
installed.  You can download and get installation instructions for the
GNU Fortran compiler here:

  http://gcc.gnu.org/install/

To compile the Java portion of the program, you will need to have the Java 6
Development Kit (JDK) installed.  You can download the JDK here:

  http://java.sun.com/javase/downloads/index.jsp

You will also need to have Apache Ant installed for the compilation of the
Java program.  You can get installation instructions and download binaries
here:

  http://ant.apache.org/manual/install.html

Windows
==============================

The Makefile was created to compile these programs in a POSIX environment with
access to the GNU Make system.  If you wish to compile these programs in a
Windows environment using the provided Makefile, you will need to install
GNU Make and GNU CoreUtils in addition to the other requirements.  You can
download the binaries here:

  http://gnuwin32.sourceforge.net/packages.html

Cygwin is also an acceptable POSIX environment for Windows.  You can get
installation instructions, and download the binaries here:

http://cygwin.com/install.html

You will also need to have the Pthreads-w32 library installed.  You can get
installation instructions, and download the binaries here:

  http://sourceware.org/pthreads-win32/


BUILDING THE SOURCE
================================================================================

Make
------------------------------

To build the binary files, and create the jar, issue the command:

  make

To build just the Fortran source files, issue the command:

  make fortran

To build just the Java source files, issue the command:

  make java

To clean the directory, issue the command:

  make clean


