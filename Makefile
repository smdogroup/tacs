# ============================================
#
# Make file for TACS_DIR/
#
# ============================================

include Makefile.in
include TACS_Common.mk

# When A2D is taken from its default submodule location, fetch the submodule
# automatically before building so a plain (non-recursive) clone still works,
# and warn if the checkout no longer matches the commit pinned by TACS (e.g.
# after a git pull). The stale case is deliberately not synced automatically
# so that a local A2D checkout is never clobbered.
ifeq (${A2D_DIR},${TACS_DIR}/extern/a2d)
default debug: a2d-submodule

.PHONY: a2d-submodule
a2d-submodule: ${A2D_DIR}/include/a2dcore.h
	@case "$$(cd ${TACS_DIR} && git submodule status extern/a2d 2>/dev/null)" in \
	   +*) echo "WARNING: extern/a2d does not match the A2D commit pinned by TACS."; \
	       echo "Run 'git submodule update' from ${TACS_DIR} to sync it." ;; \
	esac

${A2D_DIR}/include/a2dcore.h:
	@if [ -e ${TACS_DIR}/.git ]; then \
	   echo "A2D submodule not initialized; running git submodule update --init extern/a2d"; \
	   cd ${TACS_DIR} && git submodule update --init extern/a2d; \
	else \
	   echo "ERROR: A2D headers not found at ${A2D_DIR}/include and this is not a"; \
	   echo "git checkout, so the extern/a2d submodule cannot be fetched. Download"; \
	   echo "https://github.com/smdogroup/a2d and set A2D_DIR in Makefile.in."; \
	   exit 1; \
	fi
endif

TACS_SUBDIRS = src \
	src/bpmat \
	src/elements \
	src/elements/dynamics \
	src/elements/basis \
	src/elements/shell \
	src/constitutive \
	src/functions \
	src/io

TACS_OBJS := $(addsuffix /*.o, ${TACS_SUBDIRS})

default:
	@if [ "${TACS_IS_COMPLEX}" = "true" ]; then \
	   echo "Building Complex TACS"; \
	   for subdir in $(TACS_SUBDIRS) ; do \
	      echo "making $@ in $$subdir"; \
	      echo; (cd $$subdir && $(MAKE) TACS_DIR=${TACS_DIR} TACS_DEF="${TACS_DEF} -DTACS_USE_COMPLEX") || exit 1; \
            done \
	else \
	   echo "Building Real TACS"; \
	   for subdir in $(TACS_SUBDIRS) ; do \
	      echo "making $@ in $$subdir"; \
	      echo; (cd $$subdir && $(MAKE) TACS_DIR=${TACS_DIR} TACS_DEF=${TACS_DEF}) || exit 1; \
            done \
	fi
	${CXX} ${SO_LINK_FLAGS} ${TACS_OBJS} ${TACS_EXTERN_LIBS} -o ${TACS_DIR}/lib/libtacs.${SO_EXT}
	@if [ "${TACS_IS_COMPLEX}" = "true" ]; then \
		echo "ctypedef complex TacsScalar" > tacs/cpp_headers/TacsTypedefs.pxi; \
		echo "TACS_NPY_SCALAR = np.NPY_CDOUBLE" > tacs/TacsDefs.pxi; \
		echo "dtype = complex" >> tacs/TacsDefs.pxi; \
	else \
		echo "ctypedef double TacsScalar" > tacs/cpp_headers/TacsTypedefs.pxi; \
		echo "TACS_NPY_SCALAR = np.NPY_DOUBLE" > tacs/TacsDefs.pxi; \
		echo "dtype = np.double" >> tacs/TacsDefs.pxi; \
	fi

debug:
	@if [ "${TACS_IS_COMPLEX}" = "true" ]; then \
	   echo "Building Complex TACS"; \
	   for subdir in $(TACS_SUBDIRS) ; do \
	      echo "making $@ in $$subdir"; \
	      echo; (cd $$subdir && $(MAKE) debug TACS_DIR=${TACS_DIR} TACS_DEF="${TACS_DEF} -DTACS_USE_COMPLEX") || exit 1; \
            done \
	else \
	   echo "Building Real TACS"; \
	   for subdir in $(TACS_SUBDIRS) ; do \
	      echo "making $@ in $$subdir"; \
	      echo; (cd $$subdir && $(MAKE) debug TACS_DIR=${TACS_DIR} TACS_DEF=${TACS_DEF}) || exit 1; \
            done \
	fi
	${CXX} ${SO_LINK_FLAGS} ${TACS_OBJS} ${TACS_EXTERN_LIBS} -o ${TACS_DIR}/lib/libtacs.${SO_EXT}
	@if [ "${TACS_IS_COMPLEX}" = "true" ]; then \
		echo "ctypedef complex TacsScalar" > tacs/cpp_headers/TacsTypedefs.pxi; \
		echo "TACS_NPY_SCALAR = np.NPY_CDOUBLE" > tacs/TacsDefs.pxi; \
		echo "dtype = complex" >> tacs/TacsDefs.pxi; \
	else \
		echo "ctypedef double TacsScalar" > tacs/cpp_headers/TacsTypedefs.pxi; \
		echo "TACS_NPY_SCALAR = np.NPY_DOUBLE" > tacs/TacsDefs.pxi; \
		echo "dtype = np.double" >> tacs/TacsDefs.pxi; \
	fi

interface:
	${PYTHON} -m pip install -e .\[all\];

complex_interface:
	CFLAGS=-DTACS_USE_COMPLEX ${PYTHON} -m pip install -e .\[all\];

complex: TACS_IS_COMPLEX=true
complex: default

complex_debug: TACS_IS_COMPLEX=true
complex_debug: debug

clean:
	${RM} lib/libtacs.a lib/libtacs.s${SO_EXT}
	${RM} tacs/*.${SO_EXT} tacs/*.cpp
	@for subdir in $(TACS_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; \
	     (cd $$subdir && $(MAKE) $@ TACS_DIR=${TACS_DIR}) || exit 1; \
	done
