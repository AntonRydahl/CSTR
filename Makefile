# Makefile for the CSTR example
#
default: project

CC = gcc

INCLUDES = -I.
OBJS = MersenneTwister.o RandomProcesses.o CSTR.o ImplicitEulerSolver.o

# Gather the full list of libraries for the linker line
DEFS = -fopenmp -llapack -lblas
OPT = -O3 -ftree-vectorize -march=native -funroll-all-loops
WARN = -Wall
CFLAGS = $(DEFS) -std=c11 $(WARN) $(OPT)

SOFLAGS = -shared 
XLIBS	= -L /usr/lib64/atlas -lsatlas

_DIST_HEADERS = MersenneTwister.h RandomProcesses.h CSTR.h ImplicitEulerSolver.h

MersenneTwister.o: MersenneTwister.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

RandomProcesses.o: RandomProcesses.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

CSTR.o: CSTR.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

ImplicitEulerSolver.o: ImplicitEulerSolver.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

project: project.o $(OBJS) $(XLIBS)
	$(CC) -o $@ $^ $(SOFLAGS) $(XLIBS) 

clean:
	rm -f *.o project
