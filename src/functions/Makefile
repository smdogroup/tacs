# ============================================
#
# Make file for TACS_DIR/src/functions
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = TACSFunction.o \
	Compliance.o \
	MaxFailure.o \
	KSFailure.o \
	AverageMaxFailure.o \
	AverageKSFailure.o \
	PNormFailure.o \
	InducedFailure.o \
	InducedBuckling.o \
	MaxBuckling.o \
	KSBuckling.o \
	AverageKSBuckling.o \
	PNormBuckling.o \
	StructuralMass.o \
	StructuralSize.o \
	StructuralMassMoment.o \
	AggregateDisplacement.o \
	MinimumDistance.o \
	TACSElementFunction.o \
	KSElementFailure.o

DIR=${TACS_DIR}/src/functions

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} ${CXX_OBJS}
