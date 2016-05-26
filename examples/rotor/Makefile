include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = rotor.o

default: ${OBJS}
	${CXX} -o rotor rotor.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -f *.o rotor
