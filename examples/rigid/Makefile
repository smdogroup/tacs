include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = double_pendulum.o \
	flexible_pendulum.o

default: ${OBJS}
	${CXX} -o double_pendulum double_pendulum.o ${TACS_LD_FLAGS}
	${CXX} -o flexible_pendulum flexible_pendulum.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o double_pendulum flexible_pendulum

test: default
	mpirun -np 1 ./double_pendulum

test_complex: complex
	mpirun -np 1 ./double_pendulum
