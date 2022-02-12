Interface
*********

TACS has interfaces at the C++ level and the Python level. TACS is implemented in C++,
so the interface through C++ contains all publicly accessible class member functions.
The Python level interface wraps the most important classes and functions, of which the
most frequently used are discussed below.


Workflow
--------

The most common usage of TACS is to evaluate the values and gradients of desired
structural functions with respect to specified design variables. In general, this workflow
proceeds as follows:

#. Load in a finite element model of the desired structure (in the form of a NASTRAN-style
   file) using an instance of the :class:`~TACS.MeshLoader` class.
#. For each component of the loaded model, generate an element with the desired
   cosntitutive properties and design variables.
#. Create an instance of the :class:`~TACS.Assembler` class and apply boundary conditions.
#. Solve the system and evaluate the functions and their gradients with respect to the
   design variables.

These function values and gradients can then be passed to an optimizer (such as ParOpt)
in order to minimize the value of a particular function subject to some constraints.
Improved design variable values are iteratively computed by the optimizer and Step 4 is
repeated until the optimization criteria are satisfied.