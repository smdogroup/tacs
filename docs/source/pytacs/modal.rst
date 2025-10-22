ModalProblem
------------
.. automodule:: tacs.problems.modal

Options
^^^^^^^
Options can be set for :class:`~tacs.problems.ModalProblem` at time of creation for the class in the
:meth:`pyTACS.createModalProblem <tacs.pytacs.pyTACS.createModalProblem>` method or using the
:meth:`ModalProblem.setOption <tacs.problems.ModalProblem.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`ModalProblem.printOption <tacs.problems.ModalProblem.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.problems import ModalProblem; ModalProblem.printDefaultOptions()"

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.problems.ModalProblem
  :members:
  :inherited-members: