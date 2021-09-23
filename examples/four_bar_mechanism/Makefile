include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = four_bar_mechanism.o

default: ${OBJS}
	${CXX} -o four_bar_mechanism four_bar_mechanism.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
complex_debug: debug

clean:
	rm -f *.o four_bar_mechanism

test: default
	./four_bar_mechanism

test_complex: complex
	./four_bar_mechanism
