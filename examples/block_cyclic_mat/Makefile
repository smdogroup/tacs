include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = block_cyclic_test.o

default: ${OBJS}
	${CXX} -o block_cyclic_test block_cyclic_test.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o block_cyclic_test

test: default
	mpirun -np 1 ./block_cyclic_test
	mpirun -np 3 ./block_cyclic_test
	mpirun -np 4 ./block_cyclic_test
	mpirun -np 8 ./block_cyclic_test

test_complex: complex
	mpirun -np 1 ./block_cyclic_test
	mpirun -np 3 ./block_cyclic_test
	mpirun -np 4 ./block_cyclic_test
	mpirun -np 8 ./block_cyclic_test
