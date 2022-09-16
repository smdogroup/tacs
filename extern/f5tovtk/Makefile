include ../../Makefile.in
include ../../TACS_Common.mk

default: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
default: LD_FLAGS=${TACS_LD_FLAGS}
default: f5tovtk f5tovtk_element

f5tovtk: f5tovtk.o
	${CXX} f5tovtk.o ${LD_FLAGS} -o f5tovtk

f5tovtk_element: f5tovtk_element.o
	${CXX} f5tovtk_element.o ${LD_FLAGS} -o f5tovtk_element

clean:
	${RM} f5tovtk f5tovtk_element *.o
