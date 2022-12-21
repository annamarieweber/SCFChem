# Creator:  Anna Weber (anna@scfchem.com)
# Date Created: December 21, 2022
#
# This makefile creates algoritms functions in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++14 -c  

OPTIM       = -O3

DEBUG_FLAG   = -g

# The directory where the include files needed to create the library objects are
INC_DIR = ../../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib

# targets this makefile builds
ALGO_OBJS = algoutils.o cndo2.o scf.o 

DEBUG_OBJS = debug_cndo2.o debug_scf.o debug_algoutils.o

cndo2.o: $(INC_DIR)/cndo2.h
	$(CXX) $(CXXFLAGS) $(OPTIM) cndo2.cxx -I$(INC_DIR)
	mv cndo2.o $(LIB_DIR)

debug_cndo2.o: $(INC_DIR)/cndo2.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) cndo2.cxx -I$(INC_DIR)
	mv cndo2.o $(LIB_DIR)/debug_cndo2.o

scf.o: $(INC_DIR)/scf.h
	$(CXX) $(CXXFLAGS) $(OPTIM) scf.cxx -I$(INC_DIR)
	mv scf.o $(LIB_DIR)

debug_scf.o: $(INC_DIR)/scf.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) scf.cxx -I$(INC_DIR)
	mv scf.o $(LIB_DIR)/debug_scf.o

algoutils.o: $(INC_DIR)/algoutils.h
	$(CXX) $(CXXFLAGS) $(OPTIM) algoutils.cxx -I$(INC_DIR)
	mv algoutils.o $(LIB_DIR)

debug_algoutils.o: $(INC_DIR)/algoutils.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) algoutils.cxx -I$(INC_DIR)
	mv algoutils.o $(LIB_DIR)/debug_algoutils.o


debug: $(DEBUG_OBJS)

all: $(ALGO_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(ALGO_OBJS); rm $(DEBUG_OBJS)
