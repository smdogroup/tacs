include ../../Makefile.in
include ../../TACS_Common.mk

# OBJS = double_pendulum.o \
#	flexible_pendulum.o \

OBJS = two_link_mechanism.o

default: ${OBJS}
#	${CXX} -o double_pendulum double_pendulum.o ${TACS_LD_FLAGS}
#	${CXX} -o flexible_pendulum flexible_pendulum.o ${TACS_LD_FLAGS}
	${CXX} -o two_link_mechanism two_link_mechanism.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o double_pendulum flexible_pendulum two_link_mechanism

test: default
	mpirun -np 1 ./double_pendulum

test_complex: complex
	mpirun -np 1 ./double_pendulum
