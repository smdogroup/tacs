include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = crm.o

default: ${OBJS}
	${CXX} -o crm crm.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -f *.o crm

test: default
	mpirun -np 2 ./crm
	mpirun -np 2 python crm.py 
