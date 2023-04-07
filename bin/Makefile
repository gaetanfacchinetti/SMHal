################################
# Customized Makefile
# Create the executable exe
#
# author : Gaetan Facchinetti
# gaetan.facchinetti@gmail.com
#
################################
SHELL=/bin/bash

TARGET_EXEC = exe

# Define a build directory where to build the code
BUILD_DIR = ../build
DEP_DIR = $(BUILD_DIR)/deps
OBJ_DIR = $(BUILD_DIR)/obj
SRC_DIRS = ../src
CUR_DIR = $(shell pwd)
PROJECT_DIR = $(shell dirname $(CUR_DIR))
MYLIB_DIR = ../lib

DDCalc_DIR = ../local/DDCalc/


# Define the compiler, cpp version and the standard library suited for this version
# Insert your own favorite options here. C++11 minimum is preferable
# -Wall turns on most warning
CXX = clang++
CXXFLAGS += -std=c++20 -stdlib=libc++ -Xpreprocessor -fopenmp #-Wall

# Path to the libraries (for compilation and linking) for gsl, Eigen and DDCalc
# Take into account the fact that the library can either be in /usr/local or /opt/local
LIBFLAGS += -I/opt/local/include/eigen3 -I/usr/local/include -I/usr/local/include/eigen3 -I$(PROJECT_DIR)/local/DDCalc/include -I/opt/local/include/libomp
LINKFLAGS += -L/usr/local/lib -L$(PROJECT_DIR)/local/DDCalc/lib  -L/opt/local/lib/libomp  -lgsl -lgslcblas -lDDCalc

################################################################
# Define standard colors to help visualise compilation
ifneq (,$(findstring xterm,${TERM}))
	BLACK        := $(shell tput -Txterm setaf 0)
	RED          := $(shell tput -Txterm setaf 1)
	GREEN        := $(shell tput -Txterm setaf 2)
	YELLOW       := $(shell tput -Txterm setaf 3)
	LIGHTPURPLE  := $(shell tput -Txterm setaf 4)
	PURPLE       := $(shell tput -Txterm setaf 5)
	BLUE         := $(shell tput -Txterm setaf 6)
	WHITE        := $(shell tput -Txterm setaf 7)
	RESET := $(shell tput -Txterm sgr0)
else
	BLACK        := ""
	RED          := ""
	GREEN        := ""
	YELLOW       := ""
	LIGHTPURPLE  := ""
	PURPLE       := ""
	BLUE         := ""
	WHITE        := ""
	RESET        := ""
endif
################################################################

# Find all the C++ files we want to compile in the source 
SRC = $(shell find $(SRC_DIRS) -name *.cpp)

# String substitution for every C++ file.
# As an example, hello.cpp turns into ../build/hello.o
OBJ = $(SRC:../%.cpp=$(OBJ_DIR)/%.o)

# Flag to build dependencies
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.d


# Compilation of the executable only
exec:
	@echo "----- Compilation of SMHal -----"
	@$(MAKE) $(TARGET_EXEC)
#@echo "----- Creation of the library -----"
#@$(MAKE) $(MYLIB_DIR)/libSMHal.a

.PHONY: exec

# Recompilation of everything, even DDCalc
all: 
	@echo "----- Compilation of DDCalc library -----"
	@cd $(DDCalc_DIR) && $(MAKE) --silent bin lib 
	@echo "----- Compilation of SMHal -----"
	@$(MAKE) $(TARGET_EXEC)
#@echo "----- Creation of the library -----"
#@$(MAKE) $(MYLIB_DIR)/libSMHal.a
	

.PHONY: all


$(MYLIB_DIR)/libSMHal.a: $(OBJ)
	@mkdir -p $(dir $@)
	@ar ru $@ $^
	@ranlib $@

# The final build step
$(TARGET_EXEC): $(OBJ)
	@$(CXX) $(LINKFLAGS) -lomp -o $@ $^
	@echo "$(GREEN)Compilation completed$(RESET)"


# Build step for C++ source
$(OBJ_DIR)/%.o: ../%.cpp $(DEP_DIR)/%.d | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	@mkdir -p $(dir $(DEP_DIR)/$*.d)
	@echo "$(LIGHTPURPLE)> $(BLUE)Building$(subst $(BUILD_DIR),.., $@)$(RESET)"
	@$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(LIBFLAGS) -c $< -o $@

clean:
	@$(RM) -r $(TARGET_EXEC)
	@$(RM) -r $(BUILD_DIR)/*	
	@echo "$(RED)CLEANED: build folder content and executable deleted$(RESET)"

.PHONY: clean

cleanall:
	@cd $(DDCalc_DIR) && $(MAKE) clean
	@$(MAKE) clean

.PHONY: cleanall

## show all the colors
colors: 
	@echo "$(BLACK)BLACK${RESET}"
	@echo "${RED}RED${RESET}"
	@echo "${GREEN}GREEN${RESET}"
	@echo "${YELLOW}YELLOW${RESET}"
	@echo "${LIGHTPURPLE}LIGHTPURPLE${RESET}"
	@echo "${PURPLE}PURPLE${RESET}"
	@echo "${BLUE}BLUE${RESET}"
	@echo "${WHITE}WHITE${RESET}"

# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up.
DEPFILES := $(SRC:../%.cpp=$(DEP_DIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))
