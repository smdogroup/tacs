include ../Makefile.in
include ../TACS_Common.mk

#Flags for python compilation
PY_INCLUDE = ${PYTHON_INCLUDE} -I${NUMPY_DIR} -I${MPI4PY_DIR}
CYTHON_INCLUDE = -I${NUMPY_DIR} -I${MPI4PY_DIR}

# ParOpt shared objects
CYTHON_SO = TACS.${SO_EXT} \
	constitutive.${SO_EXT} \
	elements.${SO_EXT} \
	functions.${SO_EXT}

default: ${CYTHON_SO}

complex: default

complex_debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
complex_debug: complex

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

%.cpp:%.pyx
	${CYTHON} --cplus $*.pyx

%.o:%.cpp
	${CXX} ${TACS_CC_FLAGS} ${PY_INCLUDE} -c $< -o $@

%.${SO_EXT}:%.o
	${CXX} ${SO_LINK_FLAGS} $<  ${TACS_LD_FLAGS} -o $@

clean:
	${RM} *.cpp *.${SO_EXT} *.o
