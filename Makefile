# Makefile for 2D Vortex FMM Simulation
# Compiles Fortran code for superconductor vortex dynamics

# Compiler settings
FC = gfortran
FFLAGS = -O3 -march=native -ffast-math -funroll-loops -Wall -Wextra -std=f2008
FFLAGS_DEBUG = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -std=f2008

# Target executable
TARGET = vortex_fmm
TARGET_DEBUG = vortex_fmm_debug

# Source files
SOURCES = multipole_2d_log.f90 main_vortex.f90
OBJECTS = multipole_2d_log.o main_vortex.o

# Default target
all: $(TARGET)

# Production build
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)
	@echo "==================================="
	@echo "Build successful: $(TARGET)"
	@echo "Run with: ./$(TARGET)"
	@echo "==================================="

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(TARGET_DEBUG)

$(TARGET_DEBUG): $(OBJECTS)
	$(FC) $(FFLAGS_DEBUG) -o $@ $(OBJECTS)
	@echo "==================================="
	@echo "Debug build successful: $(TARGET_DEBUG)"
	@echo "Run with: ./$(TARGET_DEBUG)"
	@echo "==================================="

# Compile module
multipole_2d_log.o: multipole_2d_log.f90
	$(FC) $(FFLAGS) -c $<

# Compile main program
main_vortex.o: main_vortex.f90 multipole_2d_log.o
	$(FC) $(FFLAGS) -c $<

# Clean build files
clean:
	rm -f *.o *.mod $(TARGET) $(TARGET_DEBUG)
	@echo "Cleaned build files"

# Clean all (including output)
cleanall: clean
	rm -rf output_*
	@echo "Cleaned all output directories"

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Help target
help:
	@echo "Available targets:"
	@echo "  make          - Build production version"
	@echo "  make debug    - Build debug version with checks"
	@echo "  make run      - Build and run simulation"
	@echo "  make clean    - Remove build files"
	@echo "  make cleanall - Remove build files and output"
	@echo "  make help     - Show this help"

.PHONY: all debug clean cleanall run help
