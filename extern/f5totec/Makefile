include ../../Makefile.in
include ../../TACS_Common.mk

default: TACS_CC_FLAGS=${TACS_OPT_CC_FLAGS}
default: f5totec

f5totec: f5totec.o
	${CXX} f5totec.o ${TACS_LD_FLAGS} -o f5totec

clean:
	${RM} f5totec *.o