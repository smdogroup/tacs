EXTERN_DIR=${shell pwd}

default: ufconfig amd metis tecio

ufconfig: # Allowed to fail since UFConfig/AMD is optional
	tar -xzf UFconfig*.tar.gz && cd UFconfig && make CFLAGS="-O3 -fPIC -fexceptions" || exit 0;

amd: # Allowed to fail since UFConfig/AMD is optional
	tar -xzf AMD*.tar.gz && cd AMD && make CFLAGS="-O3 -fPIC -fexceptions" || exit 0;

metis:
	tar -xzf metis*.tar.gz;
	cd `ls -d metis-*/` && make config prefix=${EXTERN_DIR}/metis CFLAGS="-O3 -fPIC" && make install;

tecio:
	tar -xzf tecio*.tar.gz;
	cd tecio/teciosrc && make -f Makefile.linux;

clean:
	rm -rf metis*/
	rm -rf UFconfig*/
	rm -rf AMD*/
	rm -rf tecio/teciosrc/*.o
