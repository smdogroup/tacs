NewtonSolver
-------------
.. automodule:: tacs.solvers.NewtonSolver

Options
^^^^^^^
Options can be set for :class:`~tacs.solvers.NewtonSolver` at time of creation for the class in the
:meth:`pyTACS.createStaticProblem <tacs.solvers.NewtonSolver.__init__>` method or using the
:meth:`NewtonSolver.setOption <tacs.solvers.NewtonSolver.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`NewtonSolver.printOption <tacs.solvers.NewtonSolver.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.solvers import NewtonSolver; NewtonSolver.printDefaultOptions()"

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.solvers.NewtonSolver
  :members:
  :inherited-members:
