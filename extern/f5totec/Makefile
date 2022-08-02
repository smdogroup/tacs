include ../../Makefile.in
include ../../TACS_Common.mk

default: TACS_CC_FLAGS=${TACS_OPT_CC_FLAGS} -I../tecio/teciosrc/ -I/usr/include
default: LD_FLAGS=${TACS_LD_FLAGS} -lopen-pal
default: f5totec

f5totec: f5totec.o
	${CXX} f5totec.o ../tecio/teciompisrc/libtecio.a ${LD_FLAGS} -o f5totec

clean:
	${RM} f5totec *.o
