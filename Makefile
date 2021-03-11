#TARGET	= project
#LIBSRCS	= MersenneTwister.c RandomProcesses.c CSTR.c ImplicitEulerSolver.c
#LIBOBJS	= MersenneTwister.o RandomProcesses.o CSTR.o ImplicitEulerSolver.o
#
#OPT	= -g -O3 -funroll-all-loops -march=native
#PIC	= -fPIC
#DEFS	= -fopenmp -llapack -lblas
#CC	= gcc
#WARN	= -Wall
#CFLAGS= $(WARN) $(DEFS) $(OPT) $(PIC) $(XOPTS)
#
#LDFLAGS=-lm -fopenmp -llapack -lblas
#
#LINK.o=$(CC) $(LDFLAGS) 
#
#SOFLAGS = #-shared 
#XLIBS	= #-L /usr/lib64/atlas -lsatlas -L /opt/intel/lapack-3.2.1 #-llapack -lblas
#
#$(TARGET): $(LIBOBJS)
#	$(CC) -o $@ $(SOFLAGS) $(LIBOBJS) $(XLIBS) $(LDFLAGS)
#.PHONY:
#clean:
#	@/bin/rm -f core core.* $(LIBOBJS) 
	
CXX=gcc #Chooses c compiler
CFLAGS= -Werror -fsanitize=address -g
CFLAGS= -std=c11 -Wall #compiler flags
LDLIBS=  MersenneTwister.c RandomProcesses.c CSTR.c ImplicitEulerSolver.c -lblas -lm -llapack -fopenmp#library flags
LINK.o=$(CXX) $(LDFLAGS) -L /usr/lib64/atlas -lsatlas #Makes shure to compile with the object library

### Insert targets and prerequisites below
# target: prerequisites

all: project

.PHONY: clean
clean:
		-$(RM) *.o
		-$(RM) project
