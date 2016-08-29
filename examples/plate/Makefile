include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = plate.o

default: ${OBJS}
	${CXX} -o plate plate.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o plate

test: default
	python generate_plate.py
	mpirun -np 1 ./plate NBG
	mpirun -np 1 ./plate BDF1
	mpirun -np 1 ./plate BDF2
	mpirun -np 1 ./plate BDF3
	mpirun -np 1 ./plate DIRK2
	mpirun -np 1 ./plate DIRK3
	mpirun -np 1 ./plate DIRK4
	mpirun -np 1 ./plate ABM1
