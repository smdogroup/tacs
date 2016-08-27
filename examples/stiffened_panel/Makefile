include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = create_stiff_panel.o \
	create_skewed_plate.o \
	stiffened_panel.o \
	skewed_plate.o

default: ${OBJS}
	${CXX} -o create_stiff_panel create_stiff_panel.o ${TACS_LD_FLAGS}
	${CXX} -o create_skewed_plate create_skewed_plate.o ${TACS_LD_FLAGS}
	${CXX} -o stiffened_panel stiffened_panel.o ${TACS_LD_FLAGS}
	${CXX} -o skewed_plate skewed_plate.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

clean:
	rm -f *.o create_stiff_panel create_skewed_plate stiffened_panel skewed_plate

test: default
	./create_stiff_panel
	./create_skewed_plate
	mpirun -np 4 ./stiffened_panel
	mpirun -np 4 ./skewed_plate

test_complex: default
	./create_stiff_panel
	./create_skewed_plate
	mpirun -np 4 ./stiffened_panel
	mpirun -np 4 ./skewed_plate
