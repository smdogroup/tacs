# ============================================
#
# Make file for TACS_DIR/examples
#
# ============================================

include ../Makefile.in
include ../TACS_Common.mk

EXAMPLE_SUBDIRS = \
	tutorial \
	benchmark \
	annulus \
	cylinder \
	triangle \
	crm \
	mg \
	stiffened_panel \
	grad_verify \
	profile_elements \
	rigid

default:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE)); \
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
