include ../../Makefile.in
include ../../TACS_Common.mk

#Flags for python compilation
PY_INCLUDE = ${PYTHON_INCLUDE} -I${NUMPY_DIR} -I${MPI4PY_DIR}
CYTHON_INCLUDE = -I${NUMPY_DIR} -I${MPI4PY_DIR}

# ParOpt shared objects
CYTHON_SO = shell_element.so
OBJS = shell.o

default: ${CYTHON_SO} shell

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

shell: ${OBJS}
	${CXX} -o shell shell.o ${TACS_LD_FLAGS}

%.cpp: %.pyx
	cython -I${TACS_DIR}/tacs/ --cplus $*.pyx --include-dir ${TACS_DIR}/tacs/

%.o: %.cpp
	${CXX} -I${TACS_DIR}/tacs/ ${TACS_CC_FLAGS} ${PY_INCLUDE} -c $< -o $@

shell_element.so: shell_element.o
	${CXX} -undefined dynamic_lookup ${SO_LINK_FLAGS} shell_element.o ${TACS_LD_FLAGS} -o $@

clean:
	${RM} *.so *.o
