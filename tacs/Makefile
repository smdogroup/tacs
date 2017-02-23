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

default: real_defs ${CYTHON_SO}

complex: complex_defs ${CYTHON_SO}

complex_debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
complex_debug: complex

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex_defs:
	@echo "ctypedef complex TacsScalar" > TacsTypedefs.pxi;
	@echo "TACS_NPY_SCALAR = np.NPY_CDOUBLE" > TacsDefs.pxi;
	@echo "dtype = np.complex" >> TacsDefs.pxi;

real_defs:
	@echo "ctypedef double TacsScalar" > TacsTypedefs.pxi;
	@echo "TACS_NPY_SCALAR = np.NPY_DOUBLE" > TacsDefs.pxi;
	@echo "dtype = np.double" >> TacsDefs.pxi;

%.cpp:%.pyx
	${CYTHON} --cplus $*.pyx

%.o:%.cpp
	${CXX} ${TACS_CC_FLAGS} ${PY_INCLUDE} -c $< -o $@

%.${SO_EXT}:%.o
	${CXX} ${SO_LINK_FLAGS} $<  ${TACS_LD_FLAGS} -o $@

clean:
	${RM} *.cpp *.${SO_EXT} *.o
