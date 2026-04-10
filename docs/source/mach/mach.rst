MACH
====

The MACH interface provides a specialized wrapper for coupling TACS with MDOLab's MACH framework codes, particularly for aerostructural analysis and multidisciplinary optimization. It is built around the :class:`~tacs.mach.struct_problem.StructProblem` class, which adapts a TACS :class:`~tacs.problems.StaticProblem` to conform to the ``baseclasses.StructProblem`` interface expected by MACH framework codes.

Overview
--------

The MACH interface is designed for scenarios where TACS is integrated with external computational libraries, such as:

- `ADflow <https://github.com/mdolab/adflow>`_ — aerodynamic solver for aerostructural coupling
- `pyGeo <https://github.com/mdolab/pygeo>`_ — geometric parametrization and design variables
- `pyoptsparse <https://github.com/mdolab/pyoptsparse>`_ — gradient-based optimization

The typical workflow is:

1. Create a :class:`~tacs.pytacs.pyTACS` assembler and a :class:`~tacs.problems.StaticProblem`.
2. Set up any external forces on the :class:`~tacs.problems.StaticProblem`.
3. Wrap the :class:`~tacs.problems.StaticProblem` with :class:`~tacs.mach.struct_problem.StructProblem`.
4. Optionally attach a ``DVGeometry`` object for geometric design variables.
5. At each optimization step: call :meth:`~tacs.mach.struct_problem.StructProblem.solve`, then call :meth:`~tacs.mach.struct_problem.StructProblem.evalFunctions` / :meth:`~tacs.mach.struct_problem.StructProblem.evalFunctionsSens`.

Key Features
------------

- **External Force Integration**: The :attr:`~tacs.mach.struct_problem.StructProblem.Fext` property holds the aerodynamic coupling force vector. It is automatically included in the residual during :meth:`~tacs.mach.struct_problem.StructProblem.solve`.
- **Aitken Acceleration**: :meth:`~tacs.mach.struct_problem.StructProblem.solve` optionally applies Aitken relaxation to stabilize fixed-point iteration in aerostructural coupling loops.
- **Adjoint Capabilities**: The class manages all adjoint state vectors needed for coupled sensitivity analysis, including :attr:`~tacs.mach.struct_problem.StructProblem.phi`, :attr:`~tacs.mach.struct_problem.StructProblem.pLoad`, :attr:`~tacs.mach.struct_problem.StructProblem.dSdu`, and :attr:`~tacs.mach.struct_problem.StructProblem.adjRHS`.
- **Geometric Design Variables**: When a ``DVGeometry`` object is provided (via :meth:`~tacs.mach.struct_problem.StructProblem.setDVGeo`), node coordinates are automatically updated before each analysis call and geometric sensitivities are folded in by :meth:`~tacs.mach.struct_problem.StructProblem.evalFunctionsSens`.
- **pyoptsparse Integration**: :meth:`~tacs.mach.struct_problem.StructProblem.addVariablesPyOpt` and :meth:`~tacs.mach.struct_problem.StructProblem.addConstraintsPyOpt` register structural design variables and constraints directly with a ``pyoptsparse`` optimization problem.
- **Force File I/O**: External coupling forces can be saved to and loaded from BDF-format files for post-processing or restarting aerostructural analyses.

Adjoint Vector Properties
--------------------------

During coupled adjoint computation the following state vectors are managed by ``StructProblem``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Property
     - Description
   * - ``Fext``
     - External aerodynamic coupling force. Set by the aero solver before each :meth:`solve` call.
   * - ``phi``
     - Structural adjoint vector :math:`\phi`. Updated by :meth:`solveAdjoint`.
   * - ``pLoad``
     - Aerodynamic adjoint contribution to the structural adjoint RHS: :math:`[dA/du]^T \psi`. Supplied by the aero solver.
   * - ``dSdu``
     - Product :math:`[dS/du]^T \phi` from a coupled structural block. Supplied by the coupling layer.
   * - ``adjRHS``
     - Assembled adjoint RHS: :math:`dI/du - [dA/du]^T\psi - [dS/du]^T\phi`. Built by :meth:`assembleAdjointRHS`.
   * - ``dIdu``
     - Partial derivative of the objective with respect to state variables. Set by :meth:`setAdjointRHS` and BCs are applied automatically.
   * - ``matVecRHS`` / ``matVecSolve``
     - Scratch vectors used for matrix–vector products in the coupled direct/adjoint system.

StructProblem Class
-------------------

.. autoclass:: tacs.mach.struct_problem.StructProblem
   :members:
   :undoc-members:

Related Documentation
---------------------

- `ADflow <https://mdolab-adflow.readthedocs-hosted.com/en/latest/>`_ — ADflow documentation
- `pyGeo <https://mdolab-pygeo.readthedocs-hosted.com/en/latest/>`_ — pyGeo documentation
- `pyoptsparse <https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/>`_ — pyoptsparse documentation
