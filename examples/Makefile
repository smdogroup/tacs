# ============================================
#
# Make file for TACS_DIR/examples
#
# ============================================

include ../Makefile.in
include ../TACS_Common.mk

EXAMPLE_SUBDIRS = \
	tutorial \
	annulus \
	cylinder \
	triangle \
	crm \
	rigid

default: test

test:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) test) || exit 1; \
	done

clean:
	@for subdir in $(EXAMPLE_SUBDIRS) ; do \
	  echo "making $@ in $$subdir"; \
	  echo; (cd $$subdir && $(MAKE) clean) || exit 1; \
	done
