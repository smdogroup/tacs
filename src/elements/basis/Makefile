# ============================================
#
# Make file for TACS_DIR/src/elements/basis
#
# ============================================

include ../../../Makefile.in
include ../../../TACS_Common.mk

CXX_OBJS = TACSElementBasis.o \
	TACSLineBasis.o \
	TACSQuadBasis.o \
	TACSHexaBasis.o \
	TACSTensorProductBasisImpl.o \
	TACSTriangularBasis.o \
	TACSTetrahedralBasis.o \
	TACSQuadBernsteinBasis.o \
	TACSHexaBernsteinBasis.o \
	TACSElementMappedBasis.o \
	TACSBeamBasis.o \
	TACSMITCBasis.o \
	TACSMITCQuadBasis.o

DIR=${TACS_DIR}/src/elements/basis

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
