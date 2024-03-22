Install
*******
TACS can be readily installed on both Linux and MacOS systems.
At present, Windows installation is not supported.
Windows users are recommended to try the following alternatives for accessing TACS:

1. Install a Linux-based Virtual Machine (VM)

2. Install `Windows Subsystem for Linux <https://learn.microsoft.com/en-us/windows/wsl/install>`_

3. Install `Docker <https://www.docker.com/>`_ and then use one of the `MDOLab's public docker images <https://hub.docker.com/r/mdolab/public>`_, which include TACS.

For either 1 or 2 any of the installation methods described below should get you setup with TACS.

From Anaconda
=============

For users, TACS can be easily installed through `Anaconda <https://www.anaconda.com/>`_.
`Conda packages <https://anaconda.org/smdogroup/tacs>`_ are available for MacOS and Linux platforms.
To get started, run the following in a conda terminal:

::

    conda create -n TACS -c conda-forge python=3.9 mamba
    conda activate TACS
    mamba install -c conda-forge -c smdogroup tacs

This will create an environment named "TACS" and install the `tacs` package and all
necessary dependencies. Once installed, the user will have access to all
TACS C++/python libraries, f5tovtk, f5totec, etc through their conda environment.

From source
===========

For developers, TACS can also be installed from the source code.

Prerequisites
-------------

The following packages are required to use TACS:

* MPI
* BLAS
* LAPACK
* Metis 5.1.0

To use the python interface to TACS you will also require:

* Cython
* mpi4py
* numpy

Optional packages:

* SuiteSparse - TACS can use the approximate minimum degree (AMD) ordering routines from SuiteSparse
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
    * ``METIS_DIR``: set the location of METIS
    * ``SUITESPARSE_DIR``: *optional* set location of SuiteSparse
    * ``TECIO_DIR``: *optional* set location of TecIO

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
#. ``METIS_DIR``: set the location of METIS
#. ``SUITESPARSE_DIR``: *optional* set location of SuiteSparse
#. ``TECIO_DIR``: *optional* set location of TecIO (note you can use either the ``teciosrc`` or ``teciompisrc`` implementations)

In addition, it is recommended to copy the default python settings by copying ``setup.cfg.info`` to ``setup.cfg``.
The settings in ``setup.cfg.info`` are intended for development, if you are just going to use the code as-is,
you may wish to modify these settings.

Install dependencies
--------------------

These instructions direct you to install METIS and other dependencies in the directory ``tacs/extern``.
This location for the dependencies is not required, and indeed may not be best.
If you already have these libraries installed, simply adjust the variables in ``tacs/Makefile.in`` accordingly.

Go to the directory ``tacs/extern``. Download ``metis-5.1.0`` from `<https://src.fedoraproject.org/lookaside/pkgs/metis/metis-5.1.0.tar.gz/5465e67079419a69e0116de24fce58fe/>`_ and place the file ``metis-5.1.0.tar.gz`` there.
Note that METIS needs CMake to build and install.

Optionally, you can also place ``SuiteSparse-5.13.0.tar.gz`` (available from `<https://github.com/DrTimothyAldenDavis/SuiteSparse/releases>`_) in the same directory if you want to use the approximate minimum degree ordering routines from SuiteSparse.

Also optionally, place ``tecio.tgz`` (available from `<https://www.tecplot.com/products/tecio-library/>`_) in the same directory if you want to build ``f5totec``.
Note that TecIO requires the boost library, which can be install with ``sudo apt-get install libboost-dev`` on debian systems.

Then, to build the dependencies, simply run ``make``. If the build process ends with something like:

::

    make[2]: *** No rule to make target 'w'.  Stop.
    make[2]: Leaving directory 'SomeDirectory/tacs/extern/metis-5.1.0/build/Linux-x86_64'
    make[1]: *** [Makefile:64: install] Error 2
    make[1]: Leaving directory 'SomeDirectory/tacs/extern/metis-5.1.0'
    make: *** [Makefile:11: default] Error 1

Then try manually running ``make install`` within the ``metis-5.1.0`` directory.

Make the C++ TACS library
-------------------------

Return to the root TACS directory.
Ensure that all appropriate variables are set in ``Makefile.in``.
Make the TACS libraries by running ``make`` from the root directory.

Install postprocessing tools
----------------------------

``f5tovtk`` and ``f5totec`` are executables that convert ``.f5`` files to Paraview ``.vtk`` and ``.plt`` formats compatible with Paraview and Tecplot respectively.
After compiling the C++ TACS library, go to the subdirectory ``tacs/extern/f5tovtk`` and run ``make`` there.

``f5totec`` requires Tecplot's ``tecio`` library, the installation of which is described above.

The ``extern`` directory also contains two bash scripts, ``f5convert`` and ``f5clean``, that can be used to convert and clean ``.f5`` files.
``f5convert`` converts any ``.f5`` files that don't have an up-to-date ``.vtk`` or ``.plt`` file, and ``f5clean`` removes the ``.vtk`` or ``.plt`` file corresponding to each ``.f5`` file.
Both scripts accept a ``-s`` flag that will also convert or clean the ``.f5`` files in any subdirectories that contain ``.f5`` files.
Run ``f5convert -h`` or ``f5clean -h`` for more information.

Add the following lines to your ``.bashrc`` file to add the executables to your path:

::

    export PATH="<path to the tacs directory>/extern/f5totec:$PATH"
    export PATH="<path to the tacs directory>/extern/f5tovtk:$PATH"
    export PATH="<path to the tacs directory>/extern:$PATH"


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

To build the python interface to ``tacs``, and install all dependencies, you can use the Makefile
or you can type the following command in the root directory:

::

    pip install -e .\[all\]

or alternatively, you can use the shortcut in the ``Makefile`` and type:

::

    make interface

.. note::
  If the user is using an older version of pip (<21.3) and runs into a missing ``libtacs.so`` error when importing
  tacs in python, they may need to add the following to their pip install command ``pip install -e .\[all\] --use-feature=in-tree-build``.
  This option is on by default in newer pip versions and therefore should not be necessary.

Once this process is complete the python interface install should be complete and tacs should be importable from python.



