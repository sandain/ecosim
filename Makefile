# Fortran compiler options.
CC := gfortran
CFLAGS := -cpp -fbounds-check -fopenmp -g -lgfortran -lm -O3 -Wall -Wextra
LDFLAGS := -cpp -fbounds-check -fopenmp -g -O3 -Wall -Wextra

# Apache Ant, used to compile the Java source.
ANT := ant

# Doxygen, used to compile the documentation for the Fortran and Java code.
DOXYGEN := doxygen

# Directory in which to store the resulting binary files.
BIN_DIR := bin

# Directory in which to store the resulting documentation files.
DOCS_DIR := docs

# Source files used.
SOURCE_FILES := bruteforce.f90 hillclimb.f90 \
  npopCI.f90 demarcationCI.f90 omegaCI.f90 sigmaCI.f90

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

BUILD_DIR := build$(DIRECTORY_SEPARATOR)fortran
SOURCE_DIR := src$(DIRECTORY_SEPARATOR)fortran

# Files to be created.
INSTALL_FILES := $(patsubst %.f90, \
                   $(BIN_DIR)$(DIRECTORY_SEPARATOR)%$(BINARY_EXT), \
                   $(SOURCE_FILES) \
                 )
BINARY_FILES  := $(patsubst %.f90, \
                   $(BUILD_DIR)$(DIRECTORY_SEPARATOR)%$(BINARY_EXT), \
                   $(SOURCE_FILES) \
                 )
OBJECT_FILES  := $(patsubst %.f90, \
                   $(BUILD_DIR)$(DIRECTORY_SEPARATOR)%.o, \
                   $(INCLUDE_FILES) \
                 )
MOD_FILES     := $(patsubst %.f90, %.mod, $(INCLUDE_FILES))

# List of phony build targets.
.PHONY: all clean install uninstall docs check

# The main entry point for building.
all: $(BUILD_DIR) $(BINARY_FILES)
	$(ANT)

# Remove the old build files.
clean:
	rm -f $(BINARY_FILES) $(OBJECT_FILES) $(MOD_FILES)
	rm -Rf $(BUILD_DIR)
	rm -f Doxyfile.log
	rm -Rf $(DOCS_DIR)
	$(ANT) clean

# Install the binary files to the appropriate location.
install: $(BUILD_DIR) $(BINARY_FILES) $(BIN_DIR)
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
$(BINARY_FILES) : $(BUILD_DIR)$(DIRECTORY_SEPARATOR)%$(BINARY_EXT) : \
  $(SOURCE_DIR)$(DIRECTORY_SEPARATOR)%.f90 $(OBJECT_FILES)
	$(CC) $^ $(CFLAGS) -o $@

# Build the object files.
$(OBJECT_FILES) : $(BUILD_DIR)$(DIRECTORY_SEPARATOR)%.o : \
  $(SOURCE_DIR)$(DIRECTORY_SEPARATOR)%.f90
	$(CC) -c $^ $(LDFLAGS) -o $@

# Make sure the build directory is there.
$(BUILD_DIR):
	$(MKDIR_P) $(BUILD_DIR)

# Make sure the bin directory is there.
$(BIN_DIR):
	$(MKDIR_P) $(BIN_DIR)

