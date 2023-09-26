# ============================================
#
# Make file for TACS_DIR/src/constitutive
#
# ============================================

include ../../Makefile.in
include ../../TACS_Common.mk

CXX_OBJS = TACSMaterialProperties.o \
	TACSConstitutive.o \
	TACSConstitutiveVerification.o \
	TACSPlaneStressConstitutive.o \
	TACSSolidConstitutive.o \
	TACSShellConstitutive.o \
	TACSIsoShellConstitutive.o \
	TACSLamParamShellConstitutive.o \
	TACSGeneralMassConstitutive.o \
	TACSPointMassConstitutive.o \
	TACSGeneralSpringConstitutive.o \
	TACSDOFSpringConstitutive.o \
	TACSCompositeShellConstitutive.o \
	TACSBeamConstitutive.o \
	TACSBasicBeamConstitutive.o \
	TACSIsoTubeBeamConstitutive.o \
	TACSIsoRectangleBeamConstitutive.o \
	TACSPhaseChangeMaterialConstitutive.o \
	TACSBladeStiffenedShellConstitutive.o \
	TACSSmearedCompositeShellConstitutive.o \

DIR=${TACS_DIR}/src/constitutive

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
