TARGET = project
default: $(TARGET)

CC = gcc

CFLAGS = -O3 -ftree-vectorize -march=native -funroll-all-loops

INCLUDES = -I.
OBJS = MersenneTwister.o RandomProcesses.o
OBJS += ImplicitEulerSolver.o CSTR.o

LIBS = -lm -fopenmp -llapack -L/usr/lib64/atlas -lsatlas

_DIST_HEADERS = MersenneTwister.h RandomProcesses.h ImplicitEulerSolver.h CSTR.h

MersenneTwister.o: MersenneTwister.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

RandomProcesses.o: RandomProcesses.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

ImplicitEulerSolver.o: ImplicitEulerSolver.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

CSTR.o: CSTR.c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

$(TARGET).o: $(TARGET).c $(_DIST_HEADERS)
	$(CC) -c $< $(INCLUDES) $(CFLAGS)

project: project.o $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(TARGET)
