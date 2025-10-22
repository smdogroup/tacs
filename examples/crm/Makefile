include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = crm.o \
	crm_frequency.o

default: ${OBJS}
	${CXX} -o crm crm.o ${TACS_LD_FLAGS}
	${CXX} -o crm_frequency crm_frequency.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o crm

test: default
	mpirun -np 2 ./crm
	mpirun -np 2 python crm.py

test_complex: complex
	mpirun -np 2 ./crm
	mpirun -np 2 python crm.py
