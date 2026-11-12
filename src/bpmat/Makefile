# ============================================
#
# Make file for TACS_DIR/src/bpmat
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = AMDInterface.o \
	BCSRMat.o \
	BCSRMatFact.o \
	BCSRMatMult.o \
	BCSCMatPivot.o \
	TACSNodeMap.o \
	TACSBVec.o \
	TACSBVecDistribute.o \
	TACSBVecInterp.o \
	TACSMatDistribute.o \
	TACSParallelMat.o \
	TACSBlockCyclicMat.o \
	TACSSerialPivotMat.o \
	TACSSchurMat.o \
	KSM.o \
	GSEP.o \
	JacobiDavidson.o

DIR=${TACS_DIR}/src/bpmat

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o

