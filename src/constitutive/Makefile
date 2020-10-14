# ============================================
#
# Make file for TACS_DIR/src/constitutive
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = TACSConstitutive.o \
	YSlibrary.o \
	MaterialProperties.o \
	SolidStiffness.o \
	PlaneStressStiffness.o \
	FSDTStiffness.o \
	isoFSDTStiffness.o \
	compFSDTStiffness.o \
	EBStiffness.o \
	TimoshenkoStiffness.o \
	TACSPanelAnalysis.o \
	CoupledThermoPlaneStressStiffness.o \
	CoupledThermoSolidStiffness.o \
	bladeFSDTStiffness.o

DIR=${TACS_DIR}/src/constitutive

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
