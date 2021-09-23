# ============================================
#
# Make file for TACS_DIR/src/elements/shell
#
# ============================================

include ../../../Makefile.in
include ../../../TACS_Common.mk

CXX_OBJS = RBE2.o \
	       RBE3.o

DIR=${TACS_DIR}/src/elements/shell

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
