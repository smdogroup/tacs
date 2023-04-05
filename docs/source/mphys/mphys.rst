MPhys
*****
`MPhys <https://pypi.org/project/mphys/>`_ is a package that standardizes high-fidelity multiphysics problems in OpenMDAO. Mphys eases the problem set up, provides straightforward extension to new disciplines, and has a library of OpenMDAO groups for multidisciplinary problems addressed by its standard.
The interface consists of two main groups of classes: an builder class called :class:`~tacs.mphys.mphys_tacs.TacsBuilder` and
a set of problem classes for analysis. The details of the interfaces will be discussed in the sections below.

Workflow
========

The most common usage of TACS is to evaluate the values and gradients of desired
structural functions with respect to specified design variables. Using the pyTACS interface, this workflow
proceeds as follows:

#. Load in a finite element model of the desired structure (in the form of a NASTRAN-style
   file) using an instance of the :class:`~tacs.pytacs.pyTACS` class.
#. Setup tacs element objects and design variables for the structure using the
   :meth:`pyTACS.initialize <tacs.pytacs.pyTACS.initialize>` method.
#. Create an instance of the :ref:`problem<pytacs/problems:Problem classes>` class and add loads and functions of interest.
#. Solve the :ref:`problem<pytacs/problems:Problem classes>` and evaluate the functions and their gradients with respect to the
   design variables.

These function values and gradients can then be passed to an optimizer (such as :mod:`~paropt.ParOpt`)
in order to minimize the value of a particular function subject to some constraints.
Improved design variable values are iteratively computed by the optimizer and Step 4 is
repeated until the optimization criteria are satisfied.

.. toctree::
  :maxdepth: 1

  builder