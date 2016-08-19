include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = smd.o oscillator.o

default: ${OBJS}
	${CXX} -o smd smd.o ${TACS_LD_FLAGS}
	${CXX} -o oscillator oscillator.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o smd oscillator
