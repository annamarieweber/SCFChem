# Creator:  Anna Weber (anna@scfchem.com)
# Date Created: Sept 10, 2022
#
# This makefile creates molecules objects in C++
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
CLUSTER_OBJS = cluster.o molutils.o

DEBUG_OBJS = debug_cluster.o debug_molutils.o

cluster.o: $(INC_DIR)/cluster.h
	$(CXX) $(CXXFLAGS) $(OPTIM) cluster.cxx -I$(INC_DIR)
	mv cluster.o $(LIB_DIR)

debug_cluster.o: $(INC_DIR)/cluster.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) cluster.cxx -I$(INC_DIR)
	mv cluster.o $(LIB_DIR)/debug_cluster.o

molutils.o: $(INC_DIR)/molutils.h
	$(CXX) $(CXXFLAGS) $(OPTIM) molutils.cxx -I$(INC_DIR)
	mv molutils.o $(LIB_DIR)

debug_molutils.o: $(INC_DIR)/molutils.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) molutils.cxx -I$(INC_DIR)
	mv molutils.o $(LIB_DIR)/debug_molutils.o

debug: $(DEBUG_OBJS)

all: $(CLUSTER_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(CLUSTER_OBJS); rm $(DEBUG_OBJS)

