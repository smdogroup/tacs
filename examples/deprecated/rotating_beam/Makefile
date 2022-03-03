include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = rotating_beam.o

default: ${OBJS}
	${CXX} -o rotating_beam rotating_beam.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o beam

test: default
	./rotating_beam

test_complex: complex
	./rotating_beam
