constitutive module
*******************

The `tacs.constitutive` contains classes responsible for defining constitutive behaviors for elements.

Material classes
----------------
Most constitutive classes require a material properties class to setup. The following classes are available in TACS:

.. automodule:: tacs.constitutive
  :members: MaterialProperties, OrthotropicPly
  :show-inheritance:

Constitutive classes
--------------------
Most :class:`~TACS.Element` classes require a :class:`~TACS.Constitutive` class in their setup procedure.
These objects typically defines mass, stiffness, failure and buckling behaviors for the elements that reference them.
The following :class:`~TACS.Constitutive` classes are available in TACS:

.. automodule:: tacs.constitutive
  :exclude-members: MaterialProperties, OrthotropicPly
  :members:
  :show-inheritance: