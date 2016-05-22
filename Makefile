# C compiler options.
CC := gcc
#CC := i686-w64-mingw32-gcc-win32
#CC := x86_64-w64-mingw32-gcc-win32
CCFLAGS := -DUSE_DOUBLE -DOPENMP -fopenmp -lm -O3 -finline-functions -funroll-loops -Wall

# Fortran compiler options.
FC := gfortran
#FC := i686-w64-mingw32-gfortran-win32
#FC := x86_64-w64-mingw32-gfortran-win32
FCFLAGS := -cpp -fbounds-check -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -g -lgfortran -lm -O0 -Wall
FLDFLAGS := -cpp -fbounds-check -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -g -O0 -Wall

# Apache Ant, used to compile the Java source.
ANT := ant

# Doxygen, used to compile the documentation for the Fortran and Java code.
DOXYGEN := doxygen

# Determine the operating system.
OS := $(shell uname -s)

# Set operating system specific variables.
ifneq (,$(findstring mingw, $(FC)))
  # Cross-compile for Windows.
  CCFLAGS := $(CCFLAGS) -static -lpthread
  FCFLAGS := $(FCFLAGS) -static -lpthread
  MKDIR_P := mkdir -p
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := /
else ifneq (,$(filter CYGWIN windows, $(OS)))
  # CygWin and GnuWin.
  CCFLAGS := $(CCFLAGS) -static -lpthread
  FCFLAGS := $(FCFLAGS) -static -lpthread
  MKDIR_P := mkdir
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := \\
else ifneq (,$(filter Linux, $(OS)))
  # Linux.
  MKDIR_P := mkdir -p
  BINARY_EXT :=
  DIRECTORY_SEPARATOR := /
else ifneq (,$(filter Darwin, $(OS)))
  # OS X and Darwin.
  MKDIR_P := mkdir -p
  BINARY_EXT := .app
  DIRECTORY_SEPARATOR := /
endif

BUILD_DIR := build$(DIRECTORY_SEPARATOR)
SOURCE_DIR := src$(DIRECTORY_SEPARATOR)
BIN_DIR := bin$(DIRECTORY_SEPARATOR)
DOCS_DIR := docs$(DIRECTORY_SEPARATOR)

# C configuration.
C_SOURCE_FILES := fasttree.c
C_BUILD_DIR := $(BUILD_DIR)c$(DIRECTORY_SEPARATOR)
C_SOURCE_DIR := $(SOURCE_DIR)c$(DIRECTORY_SEPARATOR)

# Fortran configuration.
FORTRAN_SOURCE_FILES := fredmethod.f90 hillclimb.f90 npopCI.f90 omegaCI.f90 sigmaCI.f90 demarcation.f90
FORTRAN_INCLUDE_FILES := darray.f90 ziggurat.f90 methods.f90 simplexmethod.f90
FORTRAN_BUILD_DIR := $(BUILD_DIR)fortran$(DIRECTORY_SEPARATOR)
FORTRAN_SOURCE_DIR := $(SOURCE_DIR)fortran$(DIRECTORY_SEPARATOR)

# Files to be created.
C_INSTALL_FILES := $(patsubst %.c, $(BIN_DIR)%$(BINARY_EXT), $(C_SOURCE_FILES))
C_BINARY_FILES  := $(patsubst %.c, $(C_BUILD_DIR)%$(BINARY_EXT), $(C_SOURCE_FILES))
FORTRAN_INSTALL_FILES := $(patsubst %.f90, $(BIN_DIR)%$(BINARY_EXT), $(FORTRAN_SOURCE_FILES))
FORTRAN_BINARY_FILES  := $(patsubst %.f90, $(FORTRAN_BUILD_DIR)%$(BINARY_EXT), $(FORTRAN_SOURCE_FILES))
FORTRAN_OBJECT_FILES  := $(patsubst %.f90, $(FORTRAN_BUILD_DIR)%.o, $(FORTRAN_INCLUDE_FILES))
FORTRAN_MOD_FILES     := $(patsubst %.f90, %.mod, $(FORTRAN_INCLUDE_FILES))

# List of phony build targets.
.PHONY: all clean install uninstall docs check

# The main entry point for building.
all: $(C_BINARY_FILES) $(FORTRAN_BINARY_FILES)
	$(ANT)

# Remove the old build files.
clean:
	rm -f Doxyfile.log
	rm -f $(FORTRAN_MOD_FILES)
	rm -Rf $(BUILD_DIR)
	rm -Rf $(DOCS_DIR)

# Install the binary files to the appropriate location.
install: $(C_BINARY_FILES) $(FORTRAN_BINARY_FILES)
	@$(MKDIR_P) $(BIN_DIR)
	cp -f $(C_BINARY_FILES) $(BIN_DIR)
	cp -f $(FORTRAN_BINARY_FILES) $(BIN_DIR)
	$(ANT) install

# Remove the binary files from their installed location.
uninstall:
	rm -f $(C_INSTALL_FILES)
	rm -f $(FORTRAN_INSTALL_FILES)
	$(ANT) uninstall

# Build the documentation.
docs:
	@$(MKDIR_P) $(DOCS_DIR)
	$(DOXYGEN) ecosim.doxy

# Run the unit tests.
check:
	$(ANT) check

# Build the c binary files.
$(C_BUILD_DIR)%$(BINARY_EXT): $(C_SOURCE_DIR)%.c
	@$(MKDIR_P) $(C_BUILD_DIR)
	$(CC) $^ $(CCFLAGS) -o $@

# Build the fortran binary files.
$(FORTRAN_BUILD_DIR)%$(BINARY_EXT): $(FORTRAN_SOURCE_DIR)%.f90 $(FORTRAN_OBJECT_FILES)
	@$(MKDIR_P) $(FORTRAN_BUILD_DIR)
	$(FC) $^ $(FCFLAGS) -o $@

# Build the fortran object files.
$(FORTRAN_BUILD_DIR)%.o: $(FORTRAN_SOURCE_DIR)%.f90
	@$(MKDIR_P) $(FORTRAN_BUILD_DIR)
	$(FC) -c $^ $(FLDFLAGS) -o $@

