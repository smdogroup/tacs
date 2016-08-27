include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = cylinder_verify.o

default: ${OBJS}
	${CXX} -o cylinder_verify cylinder_verify.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o cylinder_verify

test: default
	mpirun -np 1 ./cylinder_verify

test_complex: complex
	mpirun -np 1 ./cylinder_verify
