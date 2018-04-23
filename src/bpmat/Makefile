# ============================================
#
# Make file for TACS_DIR/src/bpmat
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = BCSRMat.o \
	BCSRMatFact.o \
	BCSRMatMult.o \
	BCSRMatFact1.o \
	BCSRMatMult1.o \
	BCSRMatFact2.o \
	BCSRMatMult2.o \
	BCSRMatFact3.o \
	BCSRMatMult3.o \
	BCSRMatFact4.o \
	BCSRMatMult4.o \
	BCSRMatFact5.o \
	BCSRMatMult5.o \
	BCSRMatFact6.o \
	BCSRMatMult6.o \
	BCSRMatFact8.o \
	BCSRMatMult8.o \
	BVec.o \
	BVecDist.o \
	BVecInterp.o \
	PMat.o \
	KSM.o \
	MatUtils.o \
	DistMat.o \
	ScMat.o \
	FEMat.o \
	GSEP.o \
	JacobiDavidson.o \
	PDMat.o \
	AMDInterface.o \
	BCSCMatPivot.o \
	SerialBCSCMat.o

DIR=${TACS_DIR}/src/bpmat

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o

