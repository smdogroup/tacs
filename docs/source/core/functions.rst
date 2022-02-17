functions module
****************

There are several structural functions built into TACS. All of them can be found in the `tacs.functions` module. These functions and their gradients with respect to design variables can be evaluated by the :class:`~TACS.Assembler` object. These functions inherit from the :class:`~TACS.Function` class.
Each class should only ever be passed to a single instance of TACS. If the function needs to be calculated for separate instances, this should be handled by separate instances of function.
The current available function in TACS:

.. automodule:: tacs.functions
  :members:
  :show-inheritance: