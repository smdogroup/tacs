elements module
***************

The `tacs.elements` module contains the full library elements supported by TACS.

Basis classes
-------------
Some :class:`~TACS.Element` classes are capable of running under different element parameterizations (number of nodes, connectivity, etc).
These elements may require an :class:`~TACS.ElementBasis` class at setup. The following :class:`~TACS.ElementBasis` classes are available in TACS:

.. automodule:: tacs.elements
  :members: LinearTetrahedralBasis, QuadraticTetrahedralBasis, CubicTetrahedralBasis,
    LinearHexaBasis, QuadraticHexaBasis, CubicHexaBasis,
    LinearQuadBasis, QuadraticQuadBasis, CubicQuadBasis, QuarticQuadBasis, QuinticQuadBasis,
    LinearTriangleBasis, QuadraticTriangleBasis, CubicTriangleBasis
  :show-inheritance:

Model classes
-------------
Some :class:`~TACS.Element` classes require :class:`~TACS.ElementModel` classes in their setup procedure.
The following :class:`~TACS.ElementModel` classes are available in TACS:

.. automodule:: tacs.elements
  :members: HeatConduction2D, HeatConduction3D,
    PCMHeatConduction2D,
    LinearElasticity2D, LinearElasticity3D,
    LinearThermoelasticity2D, LinearThermoelasticity3D
  :show-inheritance:

Transform classes
-----------------
Some :class:`~TACS.Element` classes require transform classes in their setup procedure.
The following transform classes are available in TACS:

.. automodule:: tacs.elements
  :members: ShellNaturalTransform, ShellRefAxisTransform, BeamRefAxisTransform,
    SpringIdentityTransform, SpringRefAxisTransform, SpringRefFrameTransform
  :show-inheritance:

Element classes
---------------
The following :class:`~TACS.Element` classes are available in TACS:

.. automodule:: tacs.elements
  :members: Element2D, Element3D,
    Quad4Shell, Quad9Shell, Quad16Shell, Tri3Shell,
    Quad4NonlinearShell, Quad9NonlinearShell, Quad16NonlinearShell, Tri3NonlinearShell,
    Quad4NonlinearThermalShell, Quad9NonlinearThermalShell, Quad16NonlinearThermalShell, Tri3NonlinearThermalShell,
    Quad4ThermalShell, Quad9ThermalShell, Quad16ThermalShell, Tri3ThermalShell,
    Beam2, Beam3, Beam2ModRot, Beam3ModRot,
    RBE2, RBE3, MassElement, SpringElement
  :show-inheritance: