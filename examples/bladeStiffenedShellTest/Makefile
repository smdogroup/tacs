include ../../Makefile.in
include ../../TACS_Common.mk

NAME = bladeStiffenedShellTest-Composite

OBJS = ${NAME}.o

default: ${OBJS}
	${CXX} ${TACS_CC_FLAGS} -o ${NAME} ${NAME}.o ${TACS_LD_FLAGS} -lopen-pal

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o ${NAME}

test: default
	./${NAME}

test_complex: complex
	./${NAME}
