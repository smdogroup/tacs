include ../Makefile.in

EXTERN_DIR=${shell pwd}

# --- CMake Version Check ---
# If you have CMake>=4 then you need to set CMAKE_POLICY_VERSION_MINIMUM=3.5 for METIS to configure correctly.
# This block checks the installed CMake version and sets CMAKE_POLICY_FLAG if the major version is 4 or greater. The
# flag will be empty otherwise.
# We redirect stderr to /dev/null to suppress errors if cmake is not found and then pass "version 0" to the rest of the
# command to ensure the variable is set even if cmake is not installed.
CMAKE_MAJOR_VERSION := $(shell (cmake --version 2>/dev/null || echo "version 0") | head -n 1 | awk '{print $$3}' | cut -d. -f1)

# Conditionally set the flag. The shell command inside the if will execute.
# If the condition [ ... ] is true, it echoes the flag; otherwise, it echoes nothing.
CMAKE_POLICY_FLAG := $(shell if [ $(CMAKE_MAJOR_VERSION) -ge 4 ]; then echo "CMAKE_POLICY_VERSION_MINIMUM=3.5"; fi)
# --- End CMake Version Check ---

default: suitesparse tecio metis

.phony: suitesparse
suitesparse: # Allowed to fail since SuiteSparse/AMD is optional
	tar -xzf SuiteSparse*.tar.gz && cd ${SUITESPARSE_CONFIG_DIR} && make CFLAGS="-O3 -fPIC -fexceptions" && cd ${AMD_DIR} && make CFLAGS="-O3 -fPIC -fexceptions" || exit 0;

.phony: metis
metis:
	tar -xzf metis*.tar.gz;
	cd `ls -d metis-*/` && make config prefix=${METIS_DIR} CFLAGS="-O3 -fPIC" ${CMAKE_POLICY_FLAG} && make install;

.phony: tecio
tecio: # Allowed to fail since TECIO is optional
	tar -xzf tecio*.tgz;
	cd ${TECIO_DIR} && make -f Makefile.linux || exit 0;

.phony: clean
clean:
	rm -rf SuiteSparse*/
	rm -rf metis*/
	rm -rf tecio*/
