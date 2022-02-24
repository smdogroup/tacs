pyTACS
******
The pyTACS interface utilizes a number of Python classes that automate the setup and running of models and analyses in TACS.
This interface offers the following benefits over the :ref:`core/TACS:Direct` approach:
has a more robust NASTRAN BDF mesh parsing capablity, improved interface for applying loads to structural problems,
fewer Python objects to keep track of when running typical analyses.
The interface consists of two main groups of classes: an assembler class called :class:`~tacs.pytacs.pyTACS` and
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

  pytacs_module
  problems