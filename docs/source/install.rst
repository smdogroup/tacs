Install
*******

Prerequisites
-------------

The following packages are required to use TACS:

* MPI
* BLAS
* LAPACK
* Metis 5.1

To use the python interface to TACS you will also require:

* Cython
* mpi4py
* numpy

Optional packages:

* AMD - TACS can use the approximate minimum degree ordering routines from AMD/UFConfig
* TecIO - to convert TACS FH5 output files to tecplot-compatible files

Basic steps to compile TACS
---------------------------

#. Clone the TACS git repository tacs_
#. In the base 'tacs' directory, copy the Makefile.in.info to Makefile.in. Edit
   the Makefile.in and follow the directions within to set the variables. Set
   the following:

    * ``TACS_DIR``: the root director of TACS
    * ``CXX``: the C++ compiler, must be MPI-enabled
    * ``LAPACK_LIBS``: linking arguments for the LAPACK libraries
    * ``METIS_INCLUDE`` and ``METIS_LIB``: set the include/linking arguments for METIS
    * ``AMD_INCLUDE`` and ``AMD_LIBS``: *optional* set include/linking arguments for AMD

#. To compile, from the base directory, run ``make`` then ``make interface``
#. To set up the Python interface, run ``python setup.py develop --user``

.. _tacs: https://github.com/smdogroup/tacs

Detailed installation instructions
----------------------------------

The Toolbox for the Analysis of Composite Structures (TACS) is a
parallel finite-element code written in C++, with an optional python
interface. TACS has been implemented from the start with gradient-based
optimization in mind. There are built-in routines for evaluating
functions of interest and their derivatives. TACS is object oriented
and can be extended to include new elements, constitutive properties
or functions of interest.

There are several software package dependencies required to install
TACS. The dependencies are divided into the following categories:

#. Dependencies for the C++ interface
#. Dependencies for the python interface

Checking out the code
---------------------

Using git checkout the source

::

    git clone https://github.com/smdogroup/tacs

After you have cloned TACS, copy the file ``Makefile.in.info`` to a file called ``Makefile.in``.
When compiling, TACS will look for the paths and settings in ``Makefile.in``.
Make sure to set the following:

#. ``TACS_DIR``: the root director of TACS
#. ``CXX``: the C++ compiler, must be MPI-enabled
#. ``LAPACK_LIBS``: linking arguments for the LAPACK libraries
#. ``METIS_INCLUDE`` and ``METIS_LIB``: set the include/linking arguments for METIS

In addition, it is recommended to copy the default python settings by copying ``setup.cfg.info`` to ``setup.cfg``.
The settings in ``setup.cfg.info`` are intended for development, if you are just going to use the code as-is,
you may wish to modify these settings.

Install dependencies
--------------------

These instructions direct you to install METIS and other dependencies in the directory ``tacs/extern``.
This location for the dependencies is not required, and indeed may not be best.
If you already have these libraries installed, simply adjust the variables in ``tacs/Makefile.in`` accordingly.

Go to the directory ``tacs/extern``. Download ``metis-5.1.0`` from ``http://glaros.dtc.umn.edu/gkhome/metis/metis/download`` and place the file ``metis-5.1.0.tar.gz`` there.
Note that METIS needs CMake to build and install.

::

    cd metis-5.1.0
    make config prefix=$PWD CFLAGS="-fPIC"
    make

Make the C++ TACS library
-------------------------

Return to the root TACS directory.
Ensure that all appropriate variables are set in ``Makefile.in``.
Make the TACS libraries by running ``make`` from the root directory.

f5tovtk
-------

``f5tovtk`` is an executable that converts ``.f5`` files to VTK format compatible with Paraview.
After compiling the C++ libraries, go to the subdirectory ``tacs/extern/f5tovtk`` and make the executable there.
It is useful to put this utility on your path if possible.
I add the directory ``$HOME/bin`` to my ``PATH`` and then from the directory ``$HOME/bin`` execute

::

    ln -s $HOME/git/tacs/extern/f5tovtk

Installing the python interface
-------------------------------

The python interface is generated in the ``tacs/tacs`` sub-directory.
The interface is generated using Cython.

The python interface requires the following packages:

#. ``Cython``: Python interface generator
#. ``numpy``: Numerical python packages
#. ``mpi4py``: Python interface for MPI

Use ``pip`` to install these packages if they are not already installed.
TACS works with python 3.

To build the python interface to ``tacs``, you can use the Makefile
or you can type the following command in the root directory:

::

    python setup.py build_ext --inplace

The ``--inplace`` option places the shared objects direclty in the source directories.
There are several options to make the python ``tacs`` package visible to python.
I recommend using the user install option which places a link to the library in a local directory that python looks in.
This can be performed by typing:

::

    python setup.py develop --user

Once you have installed TACS in this way, you can use the shortcut in the ``Makefile`` and type:

::

    make interface