include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = profile_elements.o

default: ${OBJS}
	${CXX} -o profile_elements profile_elements.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o profile_elements

test: default
	./profile_elements

test_complex: complex
	./profile_elements
