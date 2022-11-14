include ../Makefile.in

EXTERN_DIR=${shell pwd}

default: suitesparse tecio metis

.phony: suitesparse
suitesparse: # Allowed to fail since SuiteSparse/AMD is optional
	tar -xzf suitesparse*.tar.gz && cd ${SUITESPARSE_CONFIG_DIR} && make CFLAGS="-O3 -fPIC -fexceptions" && cd ${AMD_DIR} && make CFLAGS="-O3 -fPIC -fexceptions" || exit 0;

.phony: metis
metis:
	tar -xzf metis*.tar.gz;
	cd `ls -d metis-*/` && make config prefix=${METIS_DIR} CFLAGS="-O3 -fPIC" && make install;

.phony: tecio
tecio: # Allowed to fail since TECIO is optional
	tar -xzf tecio*.tgz;
	cd ${TECIO_DIR} && make -f Makefile.linux || exit 0;

.phony: clean
clean:
	rm -rf SuiteSparse*/
	rm -rf metis*/
	rm -rf tecio*/
