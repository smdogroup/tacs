# ============================================
#
# Make file for TACS_DIR/src/elements/dynamics
#
# ============================================

include ../../../Makefile.in
include ../../../TACS_Common.mk

CXX_OBJS = MITC3.o \
	TACSRigidBody.o \
  	TACSKinematicConstraints.o
#  	MITC9.o \

DIR=${TACS_DIR}/src/elements/dynamics/

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
