## C. Grandin, Nov 2015. Tested for g++/clang on Linux and g++ on Windows (mingw).
##
## **Make sure to set up the ADMB_HOME and ADMB_HOME_DEBUG
##  environment variables for debug and dist builds of ADMB on your system
##
## The makefiles for iscam assume you have the ADMB_HOME and ADMB_HOME_DEBUG
## environment variables set. They can be the same, but if you make debug
## using the non-debug (dist) build of ADMB, you will not be able to step into the
## ADMB source code during a debug session. You will however, still be able to
## debug the iscam.cpp file as it will be compiled with debug symbols.
##
## make-dist will make build/dist only (make sure ADMB_HOME is set)
## make-debug will make build/debug only (make sure ADMB_HOME_DEBUG is set)
## make will make both build/dist and build/debug (make sure both ADMB_HOME and ADMB_HOME_DEBUG are set)
##
## make clean-dist will clean up and remove build/dist
## make clean-debug will clean up and remove build/debug
## make clean will clean everything up and remove the build directory
##
## The examples are copied into both build/dist and build/debug.
## The data files are copied into the build/dist/bin and build/debug/bin
##  directories so that they reside with the iscam executable.
##
## The iscam.dat file is just a copy of the iscam_arrowtooth.dat file,
## so that iscam can be tested immediately after building.

COMPILER := g++

.PHONY: clean clean-debug rules testcompiler

all: dist debug

dist: testcompiler
	$(MAKE) dist --directory=src COMPILER=$(COMPILER)

debug: testcompiler
	$(MAKE) debug --directory=src COMPILER=$(COMPILER)

clean-dist:
	$(MAKE) clean-dist --directory=src

clean-debug:
	$(MAKE) clean-debug --directory=src

clean: clean-dist clean-debug
	-rm -rf build

rules: testcompiler
	$(MAKE) rules --directory=src COMPILER=$(COMPILER)

testcompiler:
ifeq (, $(shell which $(COMPILER)))
	$(error "Using shell: $(SHELL), compiler $(COMPILER) not found in path.")
	exit 1
else
	@echo "Using shell: $(SHELL), compiler: $(COMPILER) found in path."
endif
