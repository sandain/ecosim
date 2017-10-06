# C compiler options.
CC ?= gcc
#CC := i686-w64-mingw32-gcc-win32
#CC := x86_64-w64-mingw32-gcc-win32
CCFLAGS := -DUSE_DOUBLE -DOPENMP -fopenmp -lm -O3 -finline-functions -funroll-loops

# Fortran compiler options.
FC ?= gfortran
#FC := i686-w64-mingw32-gfortran-win32
#FC := x86_64-w64-mingw32-gfortran-win32
FCFLAGS := -cpp -fopenmp -lgfortran -lm -O3 -finline-functions -funroll-loops
FLDFLAGS := -cpp -fopenmp -O3 -finline-functions -funroll-loops

# Apache Ant, used to compile the Java source.
ANT := ant

# Grab the current version.
VERSION :=$(shell cat VERSION)

# The name of the zip file to build for the dist target.
DIST_ZIP := ~/ecosim-$(VERSION).zip

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
FORTRAN_SOURCE_FILES := hillclimb.f90 npopCI.f90 omegaCI.f90 sigmaCI.f90 demarcation.f90
FORTRAN_INCLUDE_FILES := darray.f90 ziggurat.f90 methods.f90 simplexmethod.f90
FORTRAN_BUILD_DIR := $(BUILD_DIR)fortran$(DIRECTORY_SEPARATOR)
FORTRAN_SOURCE_DIR := $(SOURCE_DIR)fortran$(DIRECTORY_SEPARATOR)

# Files to be created.
C_BINARY_FILES  := $(patsubst %.c, $(C_BUILD_DIR)%$(BINARY_EXT), $(C_SOURCE_FILES))
FORTRAN_BINARY_FILES  := $(patsubst %.f90, $(FORTRAN_BUILD_DIR)%$(BINARY_EXT), $(FORTRAN_SOURCE_FILES))
FORTRAN_OBJECT_FILES  := $(patsubst %.f90, $(FORTRAN_BUILD_DIR)%.o, $(FORTRAN_INCLUDE_FILES))
FORTRAN_MOD_FILES     := $(patsubst %.f90, %.mod, $(FORTRAN_INCLUDE_FILES))

# List of phony build targets.
.PHONY: all clean install uninstall docs check dist

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
	rm -Rf $(BIN_DIR)
	$(ANT) uninstall

# Build the documentation.
docs:
	@$(MKDIR_P) $(DOCS_DIR)
	$(DOXYGEN) ecosim.doxy

# Run the unit tests.
check:
	$(ANT) check

# Build the distribution zip file.
dist: uninstall install docs clean
	rm -Rf dist $(DIST_ZIP)
	@$(MKDIR_P) dist
	rsync -a --exclude='*.git*' --exclude dist . dist
	perl -i -ne 's/\%ECOSIM_VERSION\%/'$(VERSION)'/g; print;' dist/help/about.xhtml
	cd dist && zip -9 -r $(DIST_ZIP) .

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

