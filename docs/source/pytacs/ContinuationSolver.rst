ContinuationSolver
-------------
.. automodule:: tacs.solvers.ContinuationSolver

Options
^^^^^^^
Options can be set for :class:`~tacs.solvers.ContinuationSolver` at time of creation for the class in the
:meth:`pyTACS.createStaticProblem <tacs.solvers.ContinuationSolver.__init__>` method or using the
:meth:`ContinuationSolver.setOption <tacs.solvers.ContinuationSolver.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`ContinuationSolver.printOption <tacs.solvers.ContinuationSolver.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.solvers import ContinuationSolver; ContinuationSolver.printDefaultOptions()"

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.solvers.ContinuationSolver
  :members:
  :inherited-members:
