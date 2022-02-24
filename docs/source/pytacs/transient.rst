TransientProblem
----------------
.. automodule:: tacs.problems.transient

Options
^^^^^^^
Options can be set for :class:`~tacs.problems.TransientProblem` at time of creation for the class in the
:meth:`pyTACS.createTransientProblem <tacs.pytacs.pyTACS.createTransientProblem>` method or using the
:meth:`TransientProblem.setOption <tacs.problems.TransientProblem.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`TransientProblem.printOption <tacs.problems.TransientProblem.printOptions>` method.
The following options their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.problems import TransientProblem; TransientProblem.printDefaultOptions()"

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.problems.TransientProblem
  :members:
  :inherited-members: