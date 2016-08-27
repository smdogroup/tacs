include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = triangle.o

default: ${OBJS}
	${CXX} -o triangle triangle.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o triangle

test: default
	mpirun -np 1 ./triangle
	mpirun -np 2 ./triangle
	mpirun -np 4 ./triangle

test_complex: complex
	mpirun -np 1 ./triangle
	mpirun -np 2 ./triangle
	mpirun -np 4 ./triangle
