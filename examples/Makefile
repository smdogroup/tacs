# ============================================
#
# Make file for TACS_DIR/examples
#
# ============================================

include ../Makefile.in
include ../TACS_Common.mk

EXAMPLE_SUBDIRS = \
	annulus \
	mg \
        tetrahedral \
	triangle \
	tutorial \
	four_bar_mechanism \
	block_cyclic_mat

#	plate \
#	profile_elements \
#	rigid \
#	stiffened_panel \
#	benchmark \
#	crm \
#	cylinder \
#	grad_verify \

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
