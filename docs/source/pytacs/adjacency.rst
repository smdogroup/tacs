AdjacencyConstraint
-------------------
.. automodule:: tacs.constraints.adjacency

Options
^^^^^^^
Options can be set for :class:`~tacs.constraints.AdjacencyConstraint` at time of creation for the class in the
:meth:`pyTACS.createAdjacencyConstraint <tacs.pytacs.pyTACS.createAdjacencyConstraint>` method or using the
:meth:`AdjacencyConstraint.setOption <tacs.constraints.AdjacencyConstraint.setOption>` method. Current option values for a class
instance can be printed out using the :meth:`AdjacencyConstraint.printOption <tacs.constraints.AdjacencyConstraint.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs.constraints import AdjacencyConstraint; AdjacencyConstraint.printDefaultOptions()"

API Reference
^^^^^^^^^^^^^
.. autoclass:: tacs.constraints.AdjacencyConstraint
  :members:
  :inherited-members: