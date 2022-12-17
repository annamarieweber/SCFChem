# This Makefile generates the Components Needed for SCFChem++
# Author: Anna Weber - University of California, Berkeley
# Date Created: Dec 17, 2022

# This is the main build makefile

MOL_DIR = src/molecules
FUNC_DIR = src/functions
EXCEPTIONS = src/exceptions
MAIN_DIR = src/main
UTIL_DIR = src/util
DOXYFILE_PATH = doc/Doxyfile

mol:
	cd $(MOL_DIR); make all

func:
	cd $(FUNC_DIR); make debug

all:	
	cd $(EXCEPTIONS); make all
	cd $(UTIL_DIR); make all
	cd $(FUNC_DIR); make all
	cd $(MOL_DIR); make all
	cd $(MAIN_DIR); make all

docs:
	doxygen $(DOXYFILE_PATH)

debug:
	cd $(EXCEPTIONS); make debug
	cd $(UTIL_DIR); make debug

	cd $(MOL_DIR); make debug
	cd $(MAIN_DIR); make debug

run:
	cd $(MAIN_DIR); make run

cleanall:
	cd $(EXCEPTIONS); make cleanall
	cd $(MOL_DIR); make cleanall
	cd $(UTIL_DIR); make cleanall
	cd $(MAIN_DIR); make cleanall
	cd $(FUNC_DIR); make cleanall;
