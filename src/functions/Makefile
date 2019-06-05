# ============================================
#
# Make file for TACS_DIR/src/functions
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = TACSFunction.o \
	KSFailure.o \
	KSDisplacement.o \
	DisplacementIntegral.o \
	Compliance.o \
	InducedFailure.o \
	StructuralMass.o \
	ThermalKSFailure.o \
	HeatFlux.o \
	KSTemperature.o \
	KSMatTemperature.o

DIR=${TACS_DIR}/src/functions

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
