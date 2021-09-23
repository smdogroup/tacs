TMR_DIR=${HOME}/git/tmr
TACS_DIR=${HOME}/git/tacs

include ${TACS_DIR}/Makefile.in
include ${TACS_DIR}/TACS_Common.mk
include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = JD.o

%.o: %.cpp
	${CXX} ${TMR_CC_FLAGS} -I/home/yfu97/git/nlopt/install/include  -c $< -o $*.o

default: ${OBJS}
	${CXX} -o JD JD.o ${TMR_LD_FLAGS} -L/home/yfu97/git/nlopt/install/lib -lnlopt -lm -fsanitize=address -static-libasan -Wall

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS} 
debug: default


clean:
	rm -f *.o JD

