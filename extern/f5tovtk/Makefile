include ../../Makefile.in
include ../../TACS_Common.mk

default: TACS_CC_FLAGS=${TACS_OPT_CC_FLAGS}
default: f5tovtk

f5tovtk: f5tovtk.o
	${CXX} f5tovtk.o ${TACS_LD_FLAGS} -o f5tovtk

clean:
	${RM} f5tovtk *.o
