# Makefile for the CSTR example
#
default: project

CC = gcc
CFLAGS = -std=c11 -Wall -O3 -ftree-vectorize -march=native -funroll-all-loops

INCLUDES = -I.
OBJS = MersenneTwister.o RandomProcesses.o CSTR.o ImplicitEulerSolver.o

LAPACKBLAS_PATH = /zhome/0e/2/36189/teaching/02616/2019/libraries/lib
LAPACKBLAS_LIBS = -L$(LAPACKBLAS_PATH) -llapack -lblas

# Gather the full list of libraries for the linker line
LIBS = $(LAPACKBLAS_LIBS) -fopenmp

_DIST_HEADERS = MersenneTwister.h RandomProcesses.h CSTR.h ImplicitEulerSolver.h

MersenneTwister.o: MersenneTwister.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

RandomProcesses.o: RandomProcesses.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

CSTR.o: CSTR.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

ImplicitEulerSolver.o: ImplicitEulerSolver.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

project: project.o $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o project
