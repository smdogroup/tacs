include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = jacobian.o plate.o cylinder.o beam.o

default: jacobian plate cylinder beam

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

plate: ${OBJS}
	${CXX} -o plate plate.o ${TACS_LD_FLAGS}

cylinder: ${OBJS}
	${CXX} -o cylinder cylinder.o ${TACS_LD_FLAGS}

jacobian: ${OBJS}
	${CXX} -o jacobian jacobian.o ${TACS_LD_FLAGS}

beam: ${OBJS}
	${CXX} -o beam beam.o ${TACS_LD_FLAGS}

clean:
	${RM} *.o plate cylinder jacobian
