include ../../Makefile.in
include ../../TACS_Common.mk

OBJS =  two_blade.o \
	two_blade_noswash.o \
	four_blade.o \
	four_blade_noswash.o

default: ${OBJS}
	${CXX} -o two_blade two_blade.o ${TACS_LD_FLAGS}
	${CXX} -o two_blade_noswash two_blade_noswash.o ${TACS_LD_FLAGS}
	${CXX} -o four_blade four_blade.o ${TACS_LD_FLAGS}
	${CXX} -o four_blade_noswash four_blade_noswash.o ${TACS_LD_FLAGS}	

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -f *.o two_blade two_blade_noswash four_blade four_blade_noswash
