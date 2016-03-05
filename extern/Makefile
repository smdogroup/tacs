default:
	tar -xzf TACS_extern.tar.gz
	tar -xzf metis*.tar.gz
	tar -xzf UFconfig*.tar.gz
	tar -xzf AMD*.tar.gz

	echo; (cd metis* && make COPTIONS="-fPIC") || exit 1;
	echo; (cd UFconfig && make CFLAGS="-O3 -fPIC -fexceptions") || exit 1;
	echo; (cd AMD && make CFLAGS="-O3 -fPIC -fexceptions") || exit 1;
