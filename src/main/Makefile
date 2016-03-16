## C. Grandin, Nov 2015. Tested for g++/clang on Linux and g++ on Windows (mingw).
##
## **Make sure to set up the ADMB_HOME and ADMB_HOME_DEBUG
##  environment variables for debug and dist builds of ADMB on your system

ifeq ($(DEBUG),yes)
ADMB          := $(ADMB_HOME_DEBUG)
DESTDIR       := ../../build/debug/
OBJDIR        := ../../build/debug/objects/
BINDIR        := ../../build/debug/bin/
LINKERFLAGS   :=
LINKERLIBS    :=  $(ADMB)/lib/libadmb.a $(ADMB)/lib/libadmb-contrib.a
COMPILERFLAGS := -g -D__GNUDOS__ -Dlinux -DUSE_LAPLACE -std=c++11 -I. -I$(ADMB)/include -I$(ADMB)/contrib/include
LIBOBJS       := $(wildcard ../../build/debug/objects/*.o)
else
ADMB          := $(ADMB_HOME)
DESTDIR       := ../../build/dist/
OBJDIR        := ../../build/dist/objects/
BINDIR        := ../../build/dist/bin/
LINKERFLAGS   :=
LINKERLIBS    := -ladmb-contrib
COMPILERFLAGS := -c -O3 -D_FILE_OFFSET_BITS=64 -Wall -DSAFE_ALL -D__GNUDOS__ -Dlinux -DUSE_LAPLACE  -std=c++11 -I. -I$(ADMB)/include -I$(ADMB)/contrib/include
LIBOBJS       := $(wildcard ../../build/dist/objects/*.o)
endif
# LIBOBJS may contain scale.o, if so remove it from the list so linker call does not have two scale.o's.
LIBOBJS := $(filter-out ../../build/dist/objects/scale.o, $(LIBOBJS))

# Must come after the if statement so debug or dist dirs are correctly prepended
EXEC := scale
SRCS := $(EXEC).cpp
OBJS := $(EXEC).o

LINKER := g++
OBJECTS := $(addprefix $(OBJDIR), $(OBJS))

all: $(EXEC)
	@echo scale build completed, DEBUG = $(DEBUG)

$(EXEC): $(OBJECTS)
	$(LINKER) $(LINKERFLAGS) -o $@ $^ $(LIBOBJS) -L$(ADMB)/lib $(LINKERLIBS)
	@if [ ! -d "$(BINDIR)" ]; then mkdir -p $(BINDIR); fi
	-cp $(EXEC) $(BINDIR)
	-cp $(EXEC).exe $(BINDIR)

$(EXEC).cpp:
	$(ADMB)/bin/tpl2cpp $(EXEC)

$(OBJDIR)%.o: %.cpp
	@if [ ! -d "$(DESTDIR)" ]; then mkdir -p $(DESTDIR); fi
	@if [ ! -d "$(OBJDIR)" ]; then mkdir -p $(OBJDIR); fi
	$(COMPILER) $(COMPILERFLAGS) -c $< -o $@

# Show variable values and compiler information
rules:
	@echo EXEC = $(EXEC)
	@echo DESTDIR = $(DESTDIR)
	@echo OBJDIR = $(OBJDIR)
	@echo OBJS = $(OBJS)
	@echo OBJECTS = $(OBJECTS)
	@echo Library objects = $(LIBOBJS)
	@echo SRCS = $(SRCS)
	@echo ADMB = $(ADMB)
	@echo MAKE = $(MAKE)
	@echo Debug = $(DEBUG)
	@echo Compiler version = $(COMPILER) --version
	@echo Compiler flags = $(COMPILERFLAGS)
	@echo Linker = $(LINKER)
	@echo Linker flags = $(LINKERFLAGS)
	@echo Library object files = $(LIBOBJS)
