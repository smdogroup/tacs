StaticProblem
-------------
.. automodule:: tacs.problems.static

Options
^^^^^^^
Options can be set for :class:`~tacs.problems.StaticProblem` at time of creation for the class in the
:meth:`pyTACS.createStaticProblem <tacs.pytacs.pyTACS.createStaticProblem>` method or using the
:meth:`StaticProblem.setOption <tacs.problems.StaticProblem.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`StaticProblem.printOption <tacs.problems.StaticProblem.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.problems import StaticProblem; StaticProblem.printDefaultOptions()"

Nonlinear solvers
^^^^^^^^^^^^^^^^^
When you create an assembler using a nonlinear element type or constitutive model, any static problems you create will automatically setup a nonlinear solver.
A continuation solver is used to control an adaptive load incrementation process, and a Newton solver is used to solve the nonlinear system of equations at each load increment.
The options for these solvers should be set directly to ``problem.nonlinearSolver`` and ``problem.nonlinearSolver.innerSolver``.
See :ref:`ContinuationSolver <pytacs/continuation_solver:continuationsolver>` and :ref:`NewtonSolver <pytacs/newton_solver:newtonsolver>` for more information.


API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.problems.StaticProblem
  :members:
  :inherited-members:
