# ============================================
#
# Make file for TACS_DIR/src/elements
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = FElibrary.o \
	TACSElement.o \
	TACSElement2D.o \
	TACSElement3D.o \
	TACSElementTypes.o \
	TACSLinearElasticity.o \
	TACSThermoelasticity.o \
	TACSHeatConduction.o \
	TACSPlateModel.o \
	TACSThermoelasticPlateModel.o \
	TACSNeohookean.o \
	TACSTraction2D.o \
	TACSTraction3D.o \
	TACSPressure2D.o \
	TACSPressure3D.o \
	TACSInertialForce2D.o \
	TACSInertialForce3D.o \
	TACSCentrifugalForce3D.o \
	TACSConvectiveTraction2D.o \
	TACSConvectiveTraction3D.o \
	TACSElementVerification.o \
	TACSPCMHeatConduction.o \

DIR=${TACS_DIR}/src/elements

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
