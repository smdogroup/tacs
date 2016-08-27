include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = tutorial.o

default: ${OBJS}
	${CXX} -o tutorial tutorial.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o tutorial

test: default
	mpirun -np 1 ./tutorial
	mpirun -np 3 ./tutorial
	mpirun -np 4 ./tutorial
	mpirun -np 8 ./tutorial

test_complex: complex
	mpirun -np 1 ./tutorial fd=1e-30
	mpirun -np 3 ./tutorial fd=1e-30
	mpirun -np 4 ./tutorial fd=1e-30
	mpirun -np 8 ./tutorial fd=1e-30
