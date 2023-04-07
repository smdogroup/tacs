include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = spectral.o

default: ${OBJS}
	${CXX} -o spectral spectral.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o spectral

test: default
	mpirun -np 1 ./spectral
	mpirun -np 3 ./spectral
	mpirun -np 4 ./spectral
	mpirun -np 8 ./spectral

test_complex: complex
	mpirun -np 1 ./spectral
	mpirun -np 3 ./spectral
	mpirun -np 4 ./spectral
	mpirun -np 8 ./spectral
