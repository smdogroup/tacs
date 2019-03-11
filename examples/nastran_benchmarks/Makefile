include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = nastran_benchmarks.o

default: ${OBJS}
	${CXX} -o ../bin/nastran_benchmarks nastran_benchmarks.o ${TACS_LD_FLAGS} 

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o nastran_benchmarks

test: default
	mpirun -np 1 ./nastran_benchmarks
	mpirun -np 3 ./nastran_benchmarks
	mpirun -np 4 ./nastran_benchmarks
	mpirun -np 8 ./nastran_benchmarks

test_complex: complex
	mpirun -np 1 ./nastran_benchmarks fd=1e-30
	mpirun -np 3 ./nastran_benchmarks fd=1e-30
	mpirun -np 4 ./nastran_benchmarks fd=1e-30
	mpirun -np 8 ./nastran_benchmarks fd=1e-30