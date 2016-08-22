include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = tutorial.o

default: ${OBJS}
	${CXX} -o tutorial tutorial.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -f *.o tutorial
