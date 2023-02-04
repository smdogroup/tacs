[![Build, unit tests, and docs](https://github.com/smdogroup/tacs/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/smdogroup/tacs/actions/workflows/unit_tests.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

[![Anaconda-Server Badge](https://anaconda.org/smdogroup/tacs/badges/version.svg)](https://anaconda.org/smdogroup/tacs)
[![Anaconda-Server Badge](https://anaconda.org/smdogroup/tacs/badges/platforms.svg)](https://anaconda.org/smdogroup/tacs)
[![Anaconda-Server Badge](https://anaconda.org/smdogroup/tacs/badges/downloads.svg)](https://anaconda.org/smdogroup/tacs)

# TACS Overview #

The Toolkit for the Analysis of Composite Structures (TACS) is a parallel finite-element code for analysis and gradient-based design of advanced structures. Originally, TACS was primarily designed as a tool for the analysis of shell structures, such as wing-boxes. More recently it has been extended to perform topology optimization of large three-dimensional structures using gradient-based methods.

TACS has been under continuous development since 2010 by the [Structural and Multidisciplinary Design Optimization group at Georgia Tech](http://gkennedy.gatech.edu) and by the [Multidisciplinary Design Optimization Lab at the University of Michigan](http://mdolab.engin.umich.edu/).

Online documentation and examples are located at [https://smdogroup.github.io/tacs/](https://smdogroup.github.io/tacs/)

# How to cite TACS #

If you use TACS, please cite one or more of the following papers.

This paper describes the time-dependent flexible multibody dynamics and adjoint capabilities implemented in TACS:

K. Boopathy and G. J. Kennedy.  "Parallel Finite Element Framework for Rotorcraft Multibody Dynamics and Discrete Adjoint Sensitivities", 2019, https://doi.org/10.2514/1.J056585 

This paper describes the core functionality of TACS, including the adjoint-based gradient evaluation techniques it implements:

Kennedy, G. J. and Martins, J. R. R. A, "A parallel finite-element framework for large-scale gradient-based design optimization of high-performance structures", Finite Elements in Analysis and Design, 2014, doi:http://dx.doi.org/10.1016/j.finel.2014.04.011

These papers describe in detail the aggregation functional implementation in TACS:

Kennedy, G. J. and Hicken, J. E., "Improved Constraint-Aggregation Methods", Computer Methods in Applied Mechanics and Engineering, 2015, doi:http://dx.doi.org/10.1016/j.cma.2015.02.017

Kennedy, G. J., "Strategies for adaptive optimization with aggregation constraints using interior-point methods", 2015, doi:http://dx.doi.org/10.1016/j.compstruc.2015.02.024

# Setting up and installing TACS through anaconda #
The easiest way to get started with TACS is through a conda install in an [Anaconda](https://www.anaconda.com/) environment. [Conda packages](https://anaconda.org/smdogroup/tacs) are
available for MacOS and Linux platforms. To get started, run the following in a conda terminal:

    conda create -n TACS python=3.8
    conda activate TACS
    conda install -c conda-forge -c smdogroup tacs
    
This will create an environment named "TACS" and install the `tacs` package and all
necessary dependencies. Once installed, the user will have access to all TACS C++ and python libraries. 

The best way to get started is to check out and run the files in the `examples/`
folder. For instance, running the script under the `examples/crm` directory:

    mpirun python analysis.py
    
The conda install also sets up the `f5totec` and `f5tovtk` executables in the user's conda environment.
So TACS-generated .f5 solution files can be converted into a format for viewing in Tecplot or Paraview using the following respective commands:
    
    f5totec cruise_000.f5

or 

    f5tovtk cruise_000.f5

This will generate a .plt or .vtk file respectively.

# Setting up and installing TACS from source #

In addition to a working implementation of MPI, BLAS and LAPACK, TACS requires Metis 5.1 for mesh partitioning. The latest version of Metis can be obtained [here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download). TACS can optionally use the approximate minimum degree ordering routines from AMD/UFConfig. These were distributed separately, but can now be obtained from SuiteSparse package. If you use AMD, be sure to define the TACS_HAS_AMD_LIBRARY flag within the Makefile.in configuration file.

To convert TACS FH5 output files to tecplot-compatible files, you must install [TecIO](https://tecplot.azureedge.net/products/tecio/2021r2/tecio.tgz). This can be placed in the tacs/extern directory. There is also a FH5 to VTK converter as well that produces (large) ASCII files.

Once the external dependencies are installed, copy Makefile.in.info to Makefile.in. Open Makefile.in and follow the directions within to set the variables. In particular, set the following:

1. TACS_DIR: the root directory of TACS
2. CXX: the C++ compiler - must be MPI-enabled
3. LAPACK_LIBS: the linking arguments for the LAPACK libraries
4. METIS_INCLUDE/METIS_LIB: the include/linking arguments for METIS
5. AMD_INCLUDE/AMD_LIBS: the include/linking arguments for AMD

Note that the default values can often be used without modification. Of all these options, it is most important for performance reasons to link to an optimized version of LAPACK, if possible.

The C++ interface can then be compiled by running the `make` command from the tacs root directory:

   ```
   make
   ```

### Setting up the Python interface ###

The python interface can be installed after the C++ interface has been compiled. The python interface (and all dependencies) can be created with a call to setup.py. The setup.cfg.info contains the recommended defaults for the configuration script. For development, create a local development installation by executing

    pip install -e .\[all\]

This command is also executed by the command `make interface`.

If the user does not intend to modify the source code and wishes to install the interface to their python site-packages, they can instead run

    pip install .\[all\]

### Converting FH5 files ###

The utility f5totec can be used to convert the .f5 files generated by TACS to tecplot .plt files. This utility is located under tacs/extern/f5totec/. I find it convenient to create a symbolic link to f5totec in a bin directory.

### Building docs ###

In addition to the publicly hosted docs (available [here](https://smdogroup.github.io/tacs/)), users can build the TACS docs locally.
Source code for the docs are located under the `docs` directory.
To compile the docs in html format, run the following command from the `docs` directory:

  ```
  make html BUILDDIR=.
  ```

The docs can then be viewed by opening the `index.html` file:

  ```
  open ./html/index.html
  ```

# Contributing to TACS development #

We welcome users who wish to propose new bugfixes/feature additions for the TACS codebase.
All propositions for feature additions should be made as a pull request (PR) through TACS' GitHub [page](https://github.com/smdogroup/tacs).
PR's should be succinct in describing the fix/feature they intend to add, and the code changes should pertain to that specified goal (avoid feature creep).
Before finalizing their PR, we ask that developers run the following checks on their end:
1. Ensure that your modification works when TACS is compiled both in real (default) and complex mode. TACS can be compiled in complex mode using the following command: 
   ```
   make complex; make complex_interface
   ```
2. Ensure that none of TACS current capabilities have been broken by running all [unit tests](https://github.com/smdogroup/tacs/tree/master/tests) (in real and complex mode) using [Testflo](https://pypi.org/project/testflo/). 
Testflo can be run by calling the following command from TACS' root directory:
   ```
   testflo ./tests
   ```
3. Add unit/integration tests that include coverage for any features that have been added to TACS' test [library](https://github.com/smdogroup/tacs/tree/master/tests).
4. Run formatting checks on any modified C++ source/header code using [clang-format](https://clang.llvm.org/docs/ClangFormat.html).
   ```
   clang-format --style=Google -i filename.cpp
   ```
5. Run formatting checks on any modified Python code using [Black](https://black.readthedocs.io/en/stable/).
   ```
   python -m black filename.py
   ```
6. If the change is a new feature, make sure that its expected use is described through a docstring and that it is added in a relevant section of the docs
