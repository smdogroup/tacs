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
	mpirun -np 1 ./plate BDF2  &> bdf2.log

test_all: default
	python generate_plate.py
	mpirun -np 1 ./plate NBG   &> nbg.log
	mpirun -np 1 ./plate BDF1  &> bdf1.log
	mpirun -np 1 ./plate BDF2  &> bdf2.log
	mpirun -np 1 ./plate BDF3  &> bdf3.log
	mpirun -np 1 ./plate DIRK2 &> dirk2.log
	mpirun -np 1 ./plate DIRK3 &> dirk3.log
	mpirun -np 1 ./plate DIRK4 &> dirk4.log
	mpirun -np 1 ./plate ABM1  &> abm1.log
