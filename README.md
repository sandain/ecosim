Ecotype Simulation README
=========================


## INSTALLATION

Extract the contents of the archive to your hard drive:
* Windows: C:\EcoSim\
* POSIX (Linux, OSX, etc): ~/ecosim/

Create a Desktop shortcut to the shell script for your platform.
* Windows: C:\EcoSim\ecosim.bat
* POSIX (Linux, OSX, etc): ~/ecosim/ecosim.sh


## REQUIREMENTS - EXECUTION

### Binary Files

If you downloaded the Windows distribution, the necessary binary files have
already been compiled for you and can be found in the `bin` directory. Both
32-bit and 64-bit binary distributions are available, if you are running
32-bit Windows make sure you are using the 32-bit binary distribution.

If you are using Linux, OSX, BSD, or you downloaded the source distribution,
read the [REQUIREMENTS - COMPILATION](#requirements---execution) section
below.

### Java 8 JRE

A Java 8 Runtime Environment (JRE) is required to execute this program.  You
can download it here:

http://java.sun.com/javase/downloads/index.jsp


## REQUIREMENTS - COMPILATION

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

## Windows

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


## BUILDING THE SOURCE

### Make

To build the binary files, and create the jar, issue the command:

        make install

To clean the directory, issue the command:

        make clean


