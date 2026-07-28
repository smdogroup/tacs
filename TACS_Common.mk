
TACS_LIB = ${TACS_DIR}/lib/libtacs.a

# A2D is a header-only automatic-differentiation library used by newer TACS
# elements, vendored as a git submodule at extern/a2d. These defaults can be
# overridden in Makefile.in (which is always included before this file).
A2D_DIR ?= ${TACS_DIR}/extern/a2d
A2D_INCLUDE ?= -I${A2D_DIR}/include

# Warn if the A2D headers are missing (typically because the submodule has
# not been initialized). Skip the check when cleaning.
ifeq ($(filter clean,$(MAKECMDGOALS)),)
A2D_INC_DIRS = $(patsubst -I%,%,$(filter -I%,$(A2D_INCLUDE)))
ifeq ($(wildcard $(addsuffix /a2dcore.h,$(A2D_INC_DIRS))),)
$(warning ==========================================================================)
$(warning WARNING: A2D headers not found (no a2dcore.h in: $(A2D_INC_DIRS)))
$(warning TACS source files that use A2D ("a2dcore.h") will fail to compile.)
$(warning A2D is provided as a git submodule. To fetch it, run:)
$(warning >   git submodule update --init)
$(warning from the root TACS directory, or set A2D_DIR (or A2D_INCLUDE) in your)
$(warning Makefile.in to point to a checkout of https://github.com/smdogroup/a2d)
$(warning ==========================================================================)
endif
endif

TACS_INCLUDE = -I${TACS_DIR}/src \
	-I${TACS_DIR}/src/bpmat \
	-I${TACS_DIR}/src/elements \
	-I${TACS_DIR}/src/elements/a2d \
	-I${TACS_DIR}/src/elements/dynamics \
	-I${TACS_DIR}/src/elements/basis \
	-I${TACS_DIR}/src/elements/shell \
	-I${TACS_DIR}/src/constitutive \
	-I${TACS_DIR}/src/functions \
	-I${TACS_DIR}/src/io \
	${A2D_INCLUDE}

# Set the command line flags to use for compilation
TACS_OPT_CC_FLAGS = ${TACS_DEF} ${EXTRA_CC_FLAGS} ${METIS_INCLUDE} ${AMD_INCLUDE} ${TECIO_INCLUDE} ${TACS_INCLUDE}
TACS_DEBUG_CC_FLAGS = ${TACS_DEF} ${EXTRA_DEBUG_CC_FLAGS} ${METIS_INCLUDE} ${AMD_INCLUDE} ${TECIO_INCLUDE} ${TACS_INCLUDE}

# By default, use the optimized flags
TACS_CC_FLAGS = ${TACS_OPT_CC_FLAGS}

# Set the linking flags to use
TACS_EXTERN_LIBS = ${AMD_LIBS} ${METIS_LIB} ${LAPACK_LIBS} ${TECIO_LIBS}
TACS_LD_FLAGS = ${EXTRA_LD_FLAGS} ${TACS_LD_CMD} ${TACS_EXTERN_LIBS}

# This is the one rule that is used to compile all the
# source code in TACS
%.o: %.cpp
	${CXX} ${TACS_CC_FLAGS} -std=c++17 -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.cpp successfully ---"
	@echo
