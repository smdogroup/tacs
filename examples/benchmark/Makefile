include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = 	benchmark.o

default: ${OBJS}
	${CXX} -o benchmark ${OBJS} ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o benchmark

test: default
	./benchmark

test_complex: complex
	./benchmark
