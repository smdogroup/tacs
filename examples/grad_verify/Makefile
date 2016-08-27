include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = grad_verify.o

default: ${OBJS}
	${CXX} -o grad_verify grad_verify.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o grad_verify

test: default
	mpirun -np 4 ./grad_verify CRM_box_2nd.bdf

test_complex: complex
	mpirun -np 4 ./grad_verify CRM_box_2nd.bdf dh=1e-30
