# Fortran compiler options.
CC := gfortran
#CC := i686-w64-mingw32-gfortran-win32
#CC := x86_64-w64-mingw32-gfortran-win32

CFLAGS := -cpp -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -ffpe-summary=invalid,zero,overflow -g -lgfortran -lm -O0 -Wall -Wextra
LDFLAGS := -cpp -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -ffpe-summary=invalid,zero,overflow -g -O0 -Wall -Wextra
# OpenMP is currently broken with ES1.  Only enable it if you want to help debug.
#CFLAGS := -cpp -fbounds-check -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -ffpe-summary=invalid,zero,overflow -g -lgfortran -lm -O0 -Wall -Wextra
#LDFLAGS := -cpp -fbounds-check -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -ffpe-summary=invalid,zero,overflow -g -O0 -Wall -Wextra

# Apache Ant, used to compile the Java source.
ANT := ant

# Doxygen, used to compile the documentation for the Fortran and Java code.
DOXYGEN := doxygen

# Determine the operating system.
OS := $(shell uname -s)

# Set operating system specific variables.
ifneq (,$(findstring mingw, $(CC)))
  # Cross-compile for Windows.
  CFLAGS := $(CFLAGS) -static -lpthread
  MKDIR_P := mkdir -p
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := /
else ifneq (,$(findstring windows, $(OS)))
  # GnuWin.
  CFLAGS := $(CFLAGS) -static -lpthread
  MKDIR_P := mkdir
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := \\
else ifneq (,$(findstring CYGWIN, $(OS)))
  # Cygwin.
  CFLAGS := $(CFLAGS) -static -lpthread
  MKDIR_P := mkdir
  BINARY_EXT := .exe
  DIRECTORY_SEPARATOR := \\
else ifeq ($(OS), Linux)
  # Linux.
  MKDIR_P := mkdir -p
  ifeq ($(shell getconf LONG_BIT), 64)
    BINARY_EXT := .amd64
  else
    BINARY_EXT := .i386
  endif
  DIRECTORY_SEPARATOR := /
else ifeq ($(OS), Darwin)
  # OS X and Darwin.
  MKDIR_P := mkdir -p
  BINARY_EXT := .app
  DIRECTORY_SEPARATOR := /
endif

BUILD_DIR := build$(DIRECTORY_SEPARATOR)fortran$(DIRECTORY_SEPARATOR)
SOURCE_DIR := src$(DIRECTORY_SEPARATOR)fortran$(DIRECTORY_SEPARATOR)
BIN_DIR := bin$(DIRECTORY_SEPARATOR)
DOCS_DIR := docs$(DIRECTORY_SEPARATOR)

# Source files used.
SOURCE_FILES := binningdanny.f90 bruteforce.f90 correctpcr.f90 \
  demarcation.f90 divergencematrix.f90 driftCI.f90 hillclimb.f90 \
  npopCI.f90 omegaCI.f90 readsynec.f90 removegaps.f90 sigmaCI.f90

# Source files to include when building the above source files.
INCLUDE_FILES := darray.f90 tmatrix.f90 ziggurat.f90 methods.f90 simplexmethod.f90

# Files to be created.
INSTALL_FILES := $(patsubst %.f90, $(BIN_DIR)%$(BINARY_EXT), $(SOURCE_FILES))
BINARY_FILES  := $(patsubst %.f90, $(BUILD_DIR)%$(BINARY_EXT),$(SOURCE_FILES))
OBJECT_FILES  := $(patsubst %.f90,$(BUILD_DIR)%.o,$(INCLUDE_FILES))
MOD_FILES     := $(patsubst %.f90, %.mod, $(INCLUDE_FILES))

# List of phony build targets.
.PHONY: all clean install uninstall docs

# The main entry point for building.
all: $(BINARY_FILES)
	$(ANT)

# Remove the old build files.
clean:
	rm -f Doxyfile.log
	rm -f $(MOD_FILES)
	rm -Rf $(BUILD_DIR)
	rm -Rf $(DOCS_DIR)
	$(ANT) clean

# Install the binary files to the appropriate location.
install: $(BINARY_FILES)
	@$(MKDIR_P) $(BIN_DIR)
	cp -f $(BINARY_FILES) $(BIN_DIR)
	$(ANT) install

# Remove the binary files from their installed location.
uninstall:
	rm -f $(INSTALL_FILES)
	$(ANT) uninstall

# Build the documentation.
docs:
	@$(MKDIR_P) $(DOCS_DIR)
	$(DOXYGEN) ecosim.doxy

# Build the binary files.
$(BUILD_DIR)%$(BINARY_EXT): $(SOURCE_DIR)%.f90 $(OBJECT_FILES)
	@$(MKDIR_P) $(BUILD_DIR)
	$(CC) $^ $(CFLAGS) -o $@

# Build the object files.
$(BUILD_DIR)%.o: $(SOURCE_DIR)%.f90
	@$(MKDIR_P) $(BUILD_DIR)
	$(CC) -c $^ $(LDFLAGS) -o $@

