# This Makefile generates the Components Needed for SCFChem++
# Author: Anna Weber - University of California, Berkeley
# Date Created: Dec 17, 2022

# This is the main build makefile

MOL_DIR = src/molecules
ALGO_DIR = src/algorithms
FUNC_DIR = src/functions
EXCEPTIONS = src/exceptions
MAIN_DIR = src/main
UTIL_DIR = src/util
DOXYFILE_PATH = doc/Doxyfile

ifneq ($(origin LOCAL_PATHS), environment)
  LOCAL_PATHS := "-I/usr/local/include -L/usr/local/lib"
endif

ifneq ($(origin PROFILE_LIB), environment)
  PROFILE_LIB := "-lprofiler"
endif

ifneq ($(origin PROFILE_PATHS), environment)
  PROFILE_PATHS := "-I/usr/local/Cellar/gperftools/2.10/include -L/usr/local/Cellar/gperftools/2.10/lib"
endif

mol:
	cd $(MOL_DIR); make all

func:
	cd $(FUNC_DIR); make all

algo:
	cd $(ALGO_DIR); make all

all:	
	cd $(EXCEPTIONS); make all LOCAL_FLAGS=$(LOCAL_PATHS) 
	cd $(MOL_DIR); make all LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(UTIL_DIR); make all LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(FUNC_DIR); make all LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(ALGO_DIR); make all LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(MAIN_DIR); make all LOCAL_FLAGS=$(LOCAL_PATHS)

profile:	
	cd $(EXCEPTIONS); make debug LOCAL_FLAGS=$(LOCAL_PATHS) 
	cd $(MOL_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS) 
	cd $(UTIL_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(FUNC_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS) 
	cd $(ALGO_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS) 
	cd $(MAIN_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS) PROFILE=$(PROFILE_LIB) PROFILE_PATHS=$(PROFILE_PATHS)

docs:
	doxygen $(DOXYFILE_PATH)

debug:
	cd $(EXCEPTIONS); make debug LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(MOL_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(UTIL_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(FUNC_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(ALGO_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS)
	cd $(MAIN_DIR); make debug LOCAL_FLAGS=$(LOCAL_PATHS)

buildrun:
	make all;
	make run;

run:
	cd $(MAIN_DIR); make run LOCAL_FLAGS=$(LOCAL_PATHS) PROFILE=$(PROFILE_LIB) PROFILE_PATHS=$(PROFILE_PATHS)

runprofile:
	cd $(MAIN_DIR); make runprofile LOCAL_FLAGS=$(LOCAL_PATHS) PROFILE=$(PROFILE_LIB) PROFILE_PATHS=$(PROFILE_PATHS)

cleanall:
	cd $(EXCEPTIONS); make cleanall
	cd $(MOL_DIR); make cleanall
	cd $(UTIL_DIR); make cleanall
	cd $(MAIN_DIR); make cleanall
	cd $(FUNC_DIR); make cleanall;
