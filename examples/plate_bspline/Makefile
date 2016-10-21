TMR_DIR=${HOME}/hg/tmr
include ../../Makefile.in
include ../../TACS_Common.mk
include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = plate_bspline.o \
	plate_simple.o \
	plate_tmr.o \
	plate_curvilinear.o
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} ${TACS_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} -o plate plate_bspline.o ${TACS_LD_FLAGS}
	${CXX} -o test_plate plate_simple.o ${TACS_LD_FLAGS} 
	${CXX} -o tmr plate_tmr.o ${TACS_LD_FLAGS} ${TMR_LD_FLAGS}
	${CXX} -o curvilinear plate_curvilinear.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug

clean:
	rm -f *.o plate
	rm -f *.o test_plate
