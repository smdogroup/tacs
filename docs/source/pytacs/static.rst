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

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.problems.StaticProblem
  :members:
  :inherited-members: