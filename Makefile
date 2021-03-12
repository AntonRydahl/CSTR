TARGET = project
CC = gcc
DEFS = -std=c11
WARN = -Wall
SOURCES = project.c  MersenneTwister.c RandomProcesses.c CSTR.c ImplicitEulerSolver.c
LIBS = -lm -fopenmp -L/usr/lib64/atlas -lsatlas

$(CC) $(DEFS) $(WARN) -c $(SOURCES) -o $(TARGET)
