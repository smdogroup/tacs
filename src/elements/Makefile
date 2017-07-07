# ============================================
#
# Make file for TACS_DIR/src/elements
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = TACSElement.o \
	FElibrary.o \
	TACSShell.o \
	ShellUtils.o \
	LargeRotUtils.o \
	PlaneStressTri6.o \
	MITC3.o \
	MITC9.o \
	RigidBody.o \
	KinematicConstraints.o	

DIR=${TACS_DIR}/src/elements

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
