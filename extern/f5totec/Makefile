include ../../Makefile.in
include ../../TACS_Common.mk

default: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS} -I../tecio/tecsrc/ -I/usr/include
default: f5totec

f5totec: f5totec.o
	${CXX} -Wl,--no-as-needed f5totec.o ../tecio/tecio.a ${TACS_LD_FLAGS} -o f5totec

clean:
	${RM} f5totec 
