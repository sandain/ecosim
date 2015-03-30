# Fortran compiler options.
CC := gfortran
CFLAGS := -cpp -fbounds-check -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -g -lgfortran -lm -O0 -Wall -Wextra
LDFLAGS := -cpp -fbounds-check -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -g -O0 -Wall -Wextra

# Apache Ant, used to compile the Java source.
ANT := ant

# Doxygen, used to compile the documentation for the Fortran and Java code.
DOXYGEN := doxygen

# Source files used.
SOURCE_FILES := bruteforce.f90 hillclimb.f90 \
  npopCI.f90 omegaCI.f90 sigmaCI.f90 demarcation.f90

# Source files to include when building the above source files.
INCLUDE_FILES := darray.f90 ziggurat.f90 methods.f90 simplexmethod.f90

# Determine the operating system.
OS := $(shell uname -s)

# Set operating system specific variables.
ifeq ($(OS), Linux)
  MKDIR_P := mkdir -p
  ifeq ($(shell getconf LONG_BIT), 64)
    BINARY_EXT := .amd64
  else
    BINARY_EXT := .i386
  endif
  DIRECTORY_SEPARATOR := /
else ifeq ($(OS), Darwin)
  MKDIR_P := mkdir -p
  BINARY_EXT := .app
  DIRECTORY_SEPARATOR := /
else ifneq (,$(findstring windows, $(OS)))
  CFLAGS := $(CFLAGS) -static -lpthread
  MKDIR_P := mkdir
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := \\
else ifneq (,$(findstring CYGWIN, $(OS)))
  CFLAGS := $(CFLAGS) -static -lpthread
  MKDIR_P := mkdir
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := \\
endif

BUILD_DIR := build$(DIRECTORY_SEPARATOR)
SOURCE_DIR := src$(DIRECTORY_SEPARATOR)
BIN_DIR := bin$(DIRECTORY_SEPARATOR)
DOCS_DIR := docs$(DIRECTORY_SEPARATOR)

FORTRAN_BUILD_DIR := $(BUILD_DIR)fortran$(DIRECTORY_SEPARATOR)
FORTRAN_SOURCE_DIR := $(SOURCE_DIR)fortran$(DIRECTORY_SEPARATOR)

# Files to be created.
INSTALL_FILES := $(patsubst %.f90, $(BIN_DIR)%$(BINARY_EXT), $(SOURCE_FILES))
BINARY_FILES  := $(patsubst %.f90, $(FORTRAN_BUILD_DIR)%$(BINARY_EXT), $(SOURCE_FILES))
OBJECT_FILES  := $(patsubst %.f90, $(FORTRAN_BUILD_DIR)%.o, $(INCLUDE_FILES))
MOD_FILES     := $(patsubst %.f90, %.mod, $(INCLUDE_FILES))

# List of phony build targets.
.PHONY: all clean install uninstall docs check

# The main entry point for building.
all: $(BINARY_FILES)
	$(ANT)

# Remove the old build files.
clean:
	rm -f Doxyfile.log
	rm -f $(MOD_FILES)
	rm -Rf $(BUILD_DIR)
	rm -Rf $(DOCS_DIR)

# Install the binary files to the appropriate location.
install: $(BINARY_FILES)
	$(MKDIR_P) $(BIN_DIR)
	cp -f $(BINARY_FILES) $(BIN_DIR)
	$(ANT) install

# Remove the binary files from their installed location.
uninstall:
	rm -f $(INSTALL_FILES)
	$(ANT) uninstall

# Build the documentation.
docs:
	$(MKDIR_P) $(DOCS_DIR)
	$(DOXYGEN) ecosim.doxy

# Run the unit tests.
check:
	$(ANT) check

# Build the binary files.
$(FORTRAN_BUILD_DIR)%$(BINARY_EXT) : $(FORTRAN_SOURCE_DIR)%.f90 $(OBJECT_FILES)
	$(MKDIR_P) $(FORTRAN_BUILD_DIR)
	$(CC) $^ $(CFLAGS) -o $@

# Build the object files.
$(FORTRAN_BUILD_DIR)%.o : $(FORTRAN_SOURCE_DIR)%.f90
	$(MKDIR_P) $(FORTRAN_BUILD_DIR)
	$(CC) -c $^ $(LDFLAGS) -o $@

