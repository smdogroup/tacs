Installing TACS
***************

Prerequisites
-------------
The following packages are required to use TACS:

* MPI and mpi4py
* BLAS
* LAPACK
* Metis 5.1
* numpy

Optional packages:

* AMD - TACS can use the approximate minimum degree ordering routines from AMD/UFConfig
* TecIO - to convert TACS FH5 output files to tecplot-compatible files


Steps to Compile
----------------
#. Clone the TACS git repository
#. In the base 'tacs' directory, copy the Makefile.in.info to Makefile.in. Edit
   the Makefile.in and follow the directions within to set the variables. Set
   the following:

	* TACS_DIR: the root director of TACS
	* CXX: the C++ compiler, must be MPI-enabled
	* LAPACK_LIBS: linking arguments for the LAPACK libraries
	* METIS_INCLUDE/METIS_LIB: the include/linking arguments for METIS
	* AMD_INCLUDE/AMD_LIBS: optional include/linking arguments for AMD

#. To compile, from the base directory, run *make* then *make interface*
#. To set up the Python interface, run *python setup.py*
