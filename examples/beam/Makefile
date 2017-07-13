include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = beam.o

default: ${OBJS}
	${CXX} -o beam beam.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o beam

test: default
	./beam

test_complex: complex
	./beam
