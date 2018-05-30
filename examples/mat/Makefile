include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = pdmat_test.o

default: ${OBJS}
	${CXX} -o pdmat_test pdmat_test.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o pdmat_test

test: default
	mpirun -np 1 ./pdmat_test
	mpirun -np 3 ./pdmat_test
	mpirun -np 4 ./pdmat_test
	mpirun -np 8 ./pdmat_test

test_complex: complex
	mpirun -np 1 ./pdmat_test
	mpirun -np 3 ./pdmat_test
	mpirun -np 4 ./pdmat_test
	mpirun -np 8 ./pdmat_test
