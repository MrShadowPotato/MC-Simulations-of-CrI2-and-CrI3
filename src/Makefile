# Makefile for CrI2 and CrI3 Monte Carlo simulations

# Compiler settings
FC = gfortran
FFLAGS = -O3 -fcheck=all -Wall
LDFLAGS = -O3

# Directories
CORE_DIR = core
IO_DIR = io
SIMULATION_DIR = simulation
TOOLS_DIR = tools
BUILD_DIR = ../build
BIN_DIR = ../bin

# Targets
TARGETS = tLoop hLoop sLoop generator

# Source files
CORE_SRCS = $(wildcard $(CORE_DIR)/*.f90)
IO_SRCS = $(wildcard $(IO_DIR)/*.f90)
SIMULATION_SRCS = $(wildcard $(SIMULATION_DIR)/*.f90)
TOOLS_SRCS = $(wildcard $(TOOLS_DIR)/*.f90)

# Object files
CORE_OBJS = $(patsubst $(CORE_DIR)/%.f90,$(BUILD_DIR)/%.o,$(CORE_SRCS))
IO_OBJS = $(patsubst $(IO_DIR)/%.f90,$(BUILD_DIR)/%.o,$(IO_SRCS))
SIMULATION_OBJS = $(patsubst $(SIMULATION_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SIMULATION_SRCS))
TOOLS_OBJS = $(patsubst $(TOOLS_DIR)/%.f90,$(BUILD_DIR)/%.o,$(TOOLS_SRCS))

# Default target
all: setup $(TARGETS)

# Setup directories
setup:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Build rules
tLoop: $(CORE_OBJS) $(IO_OBJS) $(BUILD_DIR)/tLoop.o
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$@ $^

hLoop: $(CORE_OBJS) $(IO_OBJS) $(BUILD_DIR)/hLoop.o
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$@ $^

sLoop: $(CORE_OBJS) $(IO_OBJS) $(BUILD_DIR)/sLoop.o
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$@ $^

generator: $(CORE_OBJS) $(IO_OBJS) $(BUILD_DIR)/generator.o
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$@ $^

# Compilation rules
$(BUILD_DIR)/%.o: $(CORE_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(IO_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SIMULATION_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TOOLS_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/tLoop.o: $(SIMULATION_DIR)/tLoop.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/hLoop.o: $(SIMULATION_DIR)/hLoop.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/sLoop.o: $(SIMULATION_DIR)/sLoop.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/generator.o: $(TOOLS_DIR)/generator.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(BIN_DIR)/* 

# Remove all generated files
distclean: clean
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all setup clean distclean 