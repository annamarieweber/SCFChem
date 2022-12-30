# Creator: Anna Weber
# Date Created: November 5, 2022
#
# This makefile creates the main executable and runs it

# GNU C++ Compiler
CXX         = g++-12

# compiler flags,
CXXFLAGS    = $(LOCAL_FLAGS) $(PROFILE_PATHS) -std=c++17 -o

OPTIM       = -O3

DEBUG_FLAG   = -g

# armadillo flags
ARMAFLAGS = -larmadillo

# The directory where the include files needed to create the library objects are
INC_DIR = ../../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib

# The directory for the executables
BIN_DIR     = ../../bin

# The directory for program output
OUTPUT_DIR  = ../../output

# The directory containing the source code
SRC_DIR = ../../src

# The directory for molecule input files
MOL_FILE_DIR = ../../molecule_files

# The directory where profileing data is stored
PROFILE_DIR = ../../profile_results

# The main src directory
MAIN_DIR = $(SRC_DIR)/main

# The util src directory
UTIL_DIR = $(SRC_DIR)/util

# The algorithm src directory
ALGO_DIR = $(SRC_DIR)/algorithms

# The exceptions src directory
EXCEPT_DIR = $(SRC_DIR)/exceptions

# The functions src directory
FUNC_DIR = $(SRC_DIR)/functions

# The molecules src directory
MOL_DIR = $(SRC_DIR)/molecules

# The util header files
UTIL_HEADERS = $(INC_DIR)/filereader.h $(INC_DIR)/combinatorics.h $(INC_DIR)/timedfunctional.h

# The algorithm header files
ALGO_HEADERS = $(INC_DIR)/algoutils.h $(INC_DIR)/cndo2.h $(INC_DIR)/scf.h

# The exceptions header files
EXCEPT_HEADERS = $(INC_DIR)/clusterconfigexception.h $(INC_DIR)/baseexception.h 

# The functions header files
FUNC_HEADERS = $(INC_DIR)/shell.h $(INC_DIR)/shelloverlapintegral.h

# The molecules header files
MOL_HEADERS = $(INC_DIR)/cluster.h $(INC_DIR)/molutils.h

# The constants header files
CONST_HEADERS = $(INC_DIR)/constants.h

# All of the header files
INC_HEADERS = $(ALGO_HEADERS) $(EXCEPT_HEADERS) $(FUNC_HEADERS) $(MOL_HEADERS) $(CONST_HEADERS) $(UTIL_HEADERS) 

# The main objects
MAIN_OBJ = $(SRC_DIR)/main

# The util objects
UTIL_OBJ = $(LIB_DIR)/filereader.o $(LIB_DIR)/combinatorics.o 
DEBUG_UTIL_OBJ = $(LIB_DIR)/debug_filereader.o $(LIB_DIR)/debug_combinatorics.o 

# The algorithm objects
ALGO_OBJ = $(LIB_DIR)/algoutils.o $(LIB_DIR)/cndo2.o $(LIB_DIR)/scf.o 
DEBUG_ALGO_OBJ =  $(LIB_DIR)/debug_cndo2.o $(LIB_DIR)/debug_scf.o $(LIB_DIR)/debug_algoutils.o

# The exceptions objects
EXCEPT_OBJ = $(LIB_DIR)/baseexception.o $(LIB_DIR)/clusterconfigexception.o
DEBUG_EXCEPT_OBJ = $(LIB_DIR)/debug_baseexception.o $(LIB_DIR)/debug_clusterconfigexception.o

# The functions objects
FUNC_OBJ =  $(LIB_DIR)/shell.o $(LIB_DIR)/shelloverlapintegral.o
DEBUG_FUNC_OBJ = $(LIB_DIR)/debug_shell.o $(LIB_DIR)/debug_shelloverlapintegral.o

# The molecules objects
MOL_OBJ = $(LIB_DIR)/cluster.o $(LIB_DIR)/molutils.o 
DEBUG_MOL_OBJ = $(LIB_DIR)/debug_cluster.o $(LIB_DIR)/debug_molutils.o 

# the lib objects to include
INC_OBJS = $(INC_DIR) $(UTIL_OBJ) $(EXCEPT_OBJ) $(FUNC_OBJ) $(MOL_OBJ) $(ALGO_OBJ)

#the list of debug opjects to include
DEBUG_OBJS = $(INC_DIR) $(DEBUG_UTIL_OBJ) $(DEBUG_EXCEPT_OBJ) $(DEBUG_FUNC_OBJ) $(DEBUG_MOL_OBJ) $(DEBUG_ALGO_OBJ)

# The next line contains the list of utility objects created by this Makefile.
EXECUTABLES = main_program

main_program: $(INC_HEADERS)  
	$(CXX) $(CXXFLAGS) main_program $(OPTIM) main.cxx -I$(INC_OBJS) $(ARMAFLAGS) $(PROFILE)
	mv main_program $(BIN_DIR)

all: $(EXECUTABLES)

debug: $(INC_HEADERS) 
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) main_program main.cxx -I$(DEBUG_OBJS) $(ARMAFLAGS) $(PROFILE)
	rm -rf $(BIN_DIR)/*
	mv main_program* $(BIN_DIR)

run:
	make all
	$(BIN_DIR)/main_program $(MOL_FILE_DIR)/$(mol).txt | tee $(OUTPUT_DIR)/$(mol).out

runprofile:
	make debug
	$(BIN_DIR)/main_program $(MOL_FILE_DIR)/$(mol).txt | tee $(OUTPUT_DIR)/$(mol).out

cleanall:
	cd $(BIN_DIR); rm $(EXECUTABLES)

