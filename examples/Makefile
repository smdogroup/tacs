# ============================================
#
# Make file for TACS_DIR/examples
#
# ============================================

include ../Makefile.in
include ../TACS_Common.mk

EXAMPLE_SUBDIRS = \
	annulus \
	benchmark \
	block_cyclic_mat \
	crm \
	cylinder \
	grad_verify \
	matrix_free \
	mg \
	plate \
	profile_elements \
	shell \
	stiffened_panel \
	tetrahedral \
	triangle \
	tutorial

# Examples that should work:
# four_bar_mechanism
# pyelement
# rigid
# rotating_beam
# plate_dynamics
#

default:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE)); \
	done

debug:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) debug); \
	done

complex:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) complex); \
	done

complex_debug:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) complex_debug); \
	done

test:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) test); \
	done

test_complex:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) test_complex); \
	done

clean:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) clean); \
	done
