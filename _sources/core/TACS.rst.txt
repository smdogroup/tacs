Direct
******

In the direct Python interface for TACS, the user will be responsible for setting up and bookkeeping of most TACS objects.
This includes objects like: the MeshLoader, Assembler, state vectors, FE matrices, etc.
This approach allows for more visibility into the workings of the C++ source code, but can be daunting for new users.
For a more simplified user interface where most of this setup and bookkeeping has been automated for the user, see
:ref:`pytacs/pytacs:pyTACS`.

Workflow
--------

The most common usage of TACS is to evaluate the values and gradients of desired
structural functions with respect to specified design variables. Using the direct interface, this workflow
proceeds as follows:

#. Load in a finite element model of the desired structure (in the form of a NASTRAN-style
   file) using an instance of the :class:`~TACS.MeshLoader` class.
#. For each component of the loaded model, generate an element with the desired
   constitutive properties and design variables.
#. Create an instance of the :class:`~TACS.Assembler` class and apply boundary conditions.
#. Solve the system and evaluate the functions and their gradients with respect to the
   design variables.

These function values and gradients can then be passed to an optimizer (such as :mod:`~paropt.ParOpt`)
in order to minimize the value of a particular function subject to some constraints.
Improved design variable values are iteratively computed by the optimizer and Step 4 is
repeated until the optimization criteria are satisfied.

Assembler
---------

The :class:`~TACS.Assembler` object can be created using the :class:`~TACS.MeshLoader` or :class:`~TACS.Creator` classes. It contains the methods that allow the user to solve the finite element problem and evaluate the desired structural functions and their gradients. Once an instance has been created, its typical usage is as follows:

#. Apply loads to the model with the :func:`~TACS.Assembler.applyBCs`.
#. Create the solver using the :func:`~TACS.Assembler.createVec` and
   :func:`~TACS.Assembler.createFEMat` functions.
#. Use the :func:`~TACS.Assembler.setDesignVars` function to set design variable
   vales.
#. Evaluate structural functions (e.g. Structural Mass, KSFailure) using the
   :func:`~TACS.Assembler.evalFunctions` call.
#. The gradients of the functions with respect to the design variables can
   be evaluated using the adjoint method with the :func:`~TACS.Assembler.addDVSens`,
   :func:`~TACS.Assembler.addSVSens`, and :func:`~TACS.Assembler.addAdjointResProducts` functions.

.. autoclass:: TACS.Assembler
	:members:

MeshLoader
----------
:class:`~TACS.MeshLoader` is an interface for reading in FEM data from
NASTRAN-style files (such as .bdf files).

The typical usage for a :class:`~TACS.MeshLoader` class is as follows:

#. Create the object and call :func:`~TACS.MeshLoader.scanBDFFile` on
   the desired NASTRAN .bdf file.
#. Retrieve the number of components using :func:`~TACS.MeshLoader.getNumComponents`
#. Iterate through each component, create elements with desired constitutive properties
   and design variables. Set elements into the :class:`~TACS.MeshLoader` object and
   create the :class:`~TACS.Assembler` object.

.. autoclass:: TACS.MeshLoader
	:members:

Creator
-------

The :class:`~TACS.Creator` object is similar to the :class:`~TACS.MeshLoader` object, but
sets the nodes, elements, boundary conditions, etc. manually rather than loading them
from a NASTRAN-style file. This involves the use of the :func:`~TACS.Creator.setNodes`,
:func:`~TACS.Creator.setElements`, and :func:`~TACS.Creator.setBoundaryConditions` functions,
and finally the :func:`~TACS.Creator.createTACS` function which creates the :class:`~TACS.Assembler`
object.

.. autoclass:: TACS.Creator
	:members:

FrequencyAnalysis
-----------------

The :class:`~TACS.FrequencyAnalysis` object solves the natural frequency eigenproblem
and extracts the eigenvalues and eigenvectors (natural frequencies and mode shapes).
This could be used to, for example, minimize the mass of a beam by varying element
thicknesses with a lower bound constraint on its lowest natural frequency and an upper
bound constraint on the KSFailure function.


Integrator
----------

The :class:`~TACS.Integrator` class contains functions for solving the adjoint equations and governing equations forward in time. Classes for BDF, DIRK, ABM, and NBG integration inherit from this class.

.. autoclass:: TACS.Integrator
	:members: