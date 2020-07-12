include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = matrix_free.o

default: ${OBJS}
	${CXX} -o matrix_free matrix_free.o ${TACS_LD_FLAGS} -lprofiler

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o mg

test: default
	mpirun -np 1 ./matrix_free
	mpirun -np 3 ./matrix_free
	mpirun -np 4 ./matrix_free
	mpirun -np 8 ./matrix_free

test_complex: complex
	mpirun -np 1 ./matrix_free
	mpirun -np 3 ./matrix_free
	mpirun -np 4 ./matrix_free
	mpirun -np 8 ./matrix_free
