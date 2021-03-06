TARGET = project
default: $(TARGET)

CC = gcc

CFLAGS = -O3 -march=native

WARN = -Wall

DEFS = -std=c11 -fopenmp

OBJS = MersenneTwister.o RandomProcesses.o
OBJS += ImplicitEulerSolver.o CSTR.o

LIBS = -lm -L/usr/lib64/atlas -lsatlas

_DIST_HEADERS = MersenneTwister.h RandomProcesses.h ImplicitEulerSolver.h CSTR.h

MersenneTwister.o: MersenneTwister.c $(_DIST_HEADERS)
	$(CC) $(WARN) $(DEFS) -c $< $(INCLUDES) $(CFLAGS)

RandomProcesses.o: RandomProcesses.c $(_DIST_HEADERS)
	$(CC) $(WARN) $(DEFS) -c $< $(INCLUDES) $(CFLAGS)

ImplicitEulerSolver.o: ImplicitEulerSolver.c $(_DIST_HEADERS)
	$(CC) $(WARN) $(DEFS) -c $< $(INCLUDES) $(CFLAGS)

CSTR.o: CSTR.c $(_DIST_HEADERS)
	$(CC) $(WARN) $(DEFS) -c $< $(INCLUDES) $(CFLAGS)

$(TARGET).o: $(TARGET).c $(_DIST_HEADERS)
	$(CC) $(WARN) $(DEFS) -c $< $(INCLUDES) $(CFLAGS)

project: project.o $(OBJS)
	$(CC) $(WARN) $(DEFS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(TARGET)
