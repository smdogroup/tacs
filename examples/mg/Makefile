include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = mg.o

default: ${OBJS}
	${CXX} -o mg mg.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o mg

test: default
	mpirun -np 1 ./mg
	mpirun -np 3 ./mg
	mpirun -np 4 ./mg
	mpirun -np 8 ./mg

test_complex: complex
	mpirun -np 1 ./mg
	mpirun -np 3 ./mg
	mpirun -np 4 ./mg
	mpirun -np 8 ./mg
