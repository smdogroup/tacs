include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = tetrahedral.o

default: ${OBJS}
	${CXX} -o tetrahedral tetrahedral.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o tetrahedral
