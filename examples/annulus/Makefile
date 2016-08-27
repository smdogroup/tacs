include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = annulus.o

default: ${OBJS}
	${CXX} -o annulus annulus.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o annulus

test: default
	python generate_bdf.py --nyelems 25 --order 3
	mpirun -np 2 ./annulus circular_annulus91x51_order3.bdf
	mpirun -np 4 ./annulus circular_annulus91x51_order3.bdf

test_complex: complex
	python generate_bdf.py --nyelems 25 --order 3
	mpirun -np 2 ./annulus circular_annulus91x51_order3.bdf
	mpirun -np 4 ./annulus circular_annulus91x51_order3.bdf
