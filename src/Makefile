# ============================================
#
# Make file for TACS_DIR/src/
#
# ============================================

include ../Makefile.in
include ../TACS_Common.mk

CXX_OBJS = TACSObject.o \
	TacsUtilities.o \
	TACSAssembler.o \
	TACSAuxElements.o \
	TACSCreator.o \
	TACSMg.o \
	TACSBuckling.o \
	TACSAssembler_thread.o \
	TACSIntegrator.o \
	TACSMatrixFreeMat.o \
	TACSContinuation.o \
	TACSSpectralIntegrator.o

DIR=${TACS_DIR}/src

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
