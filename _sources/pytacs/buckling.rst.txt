BucklingProblem
---------------
.. automodule:: tacs.problems.buckling

Options
^^^^^^^
Options can be set for :class:`~tacs.problems.BucklingProblem` at time of creation for the class in the
:meth:`pyTACS.createBucklingProblem <tacs.pytacs.pyTACS.createBucklingProblem>` method or using the
:meth:`BucklingProblem.setOption <tacs.problems.BucklingProblem.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`BucklingProblem.printOption <tacs.problems.BucklingProblem.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.problems import BucklingProblem; BucklingProblem.printDefaultOptions()"

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.problems.BucklingProblem
  :members:
  :inherited-members:
