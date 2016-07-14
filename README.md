Ecotype Simulation README
=========================


## DOWNLOAD

Download the latest stable release:

https://github.com/sandain/es2/releases

Windows users: you can download 32- or 64-bit versions of Ecotype Simulation.
Make sure that you download the version appropriate for your computer.

Linux and OSX users: you must compile your own copy of Ecotype Simulation. You
will need to download one of the source distributions. See the section on
[REQUIREMENTS - COMPILATION](#requirements---compilation) for more
information.

Download the development version as a zip file:

https://github.com/sandain/es2/archive/master.zip


## INSTALLATION

Extract the contents of the archive to your hard drive:
* Windows: `C:\ecosim\`
* POSIX (Linux, OSX, etc): `~/ecosim/`

Create a Desktop shortcut to the shell script for your platform.
* Windows: `C:\ecosim\ecosim.bat`
* POSIX (Linux, OSX, etc): `~/ecosim/ecosim.sh`


## REQUIREMENTS - EXECUTION

### Binary Files

If you downloaded the Windows distribution, the necessary binary files have
already been compiled for you and can be found in the `bin` directory. Both
32-bit and 64-bit binary distributions are available, if you are running
32-bit Windows make sure you are using the 32-bit binary distribution.

If you are using Linux, OSX, BSD, or you downloaded the source distribution,
read the [REQUIREMENTS - COMPILATION](#requirements---compilation) section
below.

### Java 8 JRE

A Java 8 Runtime Environment (JRE) is required to execute this program.  You
can download it here:

http://java.sun.com/javase/downloads/index.jsp


## USAGE

You can use the shell script provided for your platform: `ecosim.bat` for
Windows and `ecosim.sh` for Linux or OSX, or you can call the `ecosim.jar`
file directly with one of the following commands:

        ecosim.bat [OPTIONS]
        ./ecosim.sh [OPTIONS]
        java -jar ecosim.jar [OPTIONS]

The following command line options are available:

        -i, --input        : A XML formated save file for input.
        -o, --output       : A XML formated save for for output.
        -s, --sequences    : A Fasta formated file for input.
        -p, --phylogeny    : A Newick formatted file for input.
        -d, --debug        : Display debugging output.
        -h, --help         : Display helpful information.
        -n, --nogui        : Hide the default GUI.  Implies --runall.
        -r, --runall       : Run everything, including demarcation.
        -t=n, --threads=n  : Set the number of threads (n) to start, default
                             to system maximum.
        -v, --version      : Display the version number.

Sequences should be aligned and in a Fasta formated file, with the outgroup
listed first.

The phylogeny should be in Newick format and must include the same leaf node
names and number as the sequences in the Fasta file.

Output is saved in XML format, and can be used to save results for later
retrieval using the input option.

Debugging messages will be printed to the console when the debug flag is used.

## REQUIREMENTS - COMPILATION

To compile the Fortran programs, you will need to have a Fortran compiler
installed.  You can download and get installation instructions for the
GNU Fortran compiler here:

http://gcc.gnu.org/install/

To compile the Java portion of the program, you will need to have the Java 8
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

If you wish to use a different build system, just make sure that the compiled
Java `ecosim.jar` file is placed in the installation folder (e.g. 
`c:\ecosim\ecosim.jar`), and the Fortran programs end up in the `bin`
directory (e.g. `c:\ecosim\bin\`).


## BUILDING THE SOURCE

### Make

To build the binary files, and create the jar, issue the command:

        make install

To clean the directory, issue the command:

        make clean


