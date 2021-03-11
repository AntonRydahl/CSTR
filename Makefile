TARGET	= project
LIBSRCS	= MersenneTwister.c RandomProcesses.c CSTR.c ImplicitEulerSolver.c
LIBOBJS	= MersenneTwister.o RandomProcesses.o CSTR.o ImplicitEulerSolver.o

OPT	= -g -O3 -funroll-all-loops -march=native
PIC	= -fPIC
DEFS	= -fopenmp
CC	= gcc
WARN	= -Wall
CFLAGS= $(DEFS) $(WARN) $(OPT) $(PIC) $(XOPTS)

LDFLAGS = -lm 

SOFLAGS = -shared 
XLIBS	= -L /usr/lib64/atlas -lsatlas -L /usr/lib64/blas -lblas -L /usr/lib64/lapack -llapack

$(TARGET): $(LIBOBJS)
	$(CC) -o $@ $(SOFLAGS) $(LIBOBJS) $(XLIBS)
.PHONY:
clean:
	@/bin/rm -f core core.* $(LIBOBJS) 
