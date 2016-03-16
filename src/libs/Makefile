## C. Grandin, Nov 2015. Tested for g++/clang on Linux and g++ on Windows (mingw).
##
## **Make sure to set up the ADMB_HOME and ADMB_HOME_DEBUG
##  environment variables for debug and dist builds of ADMB on your system

ifeq ($(DEBUG),yes)
ADMB          = $(ADMB_HOME_DEBUG)
DESTDIR       = ../../build/debug/
OBJDIR        = ../../build/debug/objects/
BINDIR        = ../../build/debug/bin/
COMPILERFLAGS = -g -D__GNUDOS__ -Dlinux -DUSE_LAPLACE -std=c++11 -I. -I$(ADMB)/include -I$(ADMB)/contrib/include
else
ADMB          = $(ADMB_HOME)
DESTDIR       = ../../build/dist/
OBJDIR        = ../../build/dist/objects/
BINDIR        = ../../build/dist/bin/
COMPILERFLAGS = -c -O2 -D_FILE_OFFSET_BITS=64 -Wall -DSAFE_ALL -D__GNUDOS__ -Dlinux -DUSE_LAPLACE  -std=c++11 -I. -I$(ADMB)/include -I$(ADMB)/contrib/include
endif

# Must come after the if statement so debug or dist dirs are correctly prepended
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

OBJECTS := $(addprefix $(OBJDIR), $(OBJS))

all: $(OBJECTS)
	@echo scale libraries build completed, DEBUG = $(DEBUG)

$(OBJDIR)%.o: %.cpp
	@if [ ! -d "$(DESTDIR)" ]; then mkdir -p $(DESTDIR); fi
	@if [ ! -d "$(OBJDIR)" ]; then mkdir -p $(OBJDIR); fi
	$(COMPILER) $(COMPILERFLAGS) -c $< -o $@

# Show variable values and compiler information
rules:
	@echo DESTDIR = $(DESTDIR)
	@echo OBJDIR = $(OBJDIR)
	@echo OBJS = $(OBJS)
	@echo OBJECTS = $(OBJECTS)
	@echo SRCS = $(SRCS)
	@echo ADMB = $(ADMB)
	@echo MAKE = $(MAKE)
	@echo Debug = $(DEBUG)
	@echo Compiler version = $(COMPILER) --version
	@echo Compiler flags = $(COMPILERFLAGS)
