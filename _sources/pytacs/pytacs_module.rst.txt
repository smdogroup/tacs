pyTACS class
============
.. automodule:: tacs.pytacs

Options
-------
Options can be set for :class:`~pyTACS` at time of creation for the class or using the
:meth:`pyTACS.setOption <tacs.pytacs.pyTACS.setOption>`. Current option values for a class
instance can be printed out using the :meth:`pyTACS.printOption <tacs.pytacs.pyTACS.printOptions>` method.
The following options, their default values and descriptions are listed below:

.. program-output:: python -c "from tacs import pyTACS; pyTACS.printDefaultOptions()"

Initializing
------------
Before the class can be used to create problems, it must first be initialized by calling the
:meth:`pyTACS.initialize <tacs.pytacs.pyTACS.initialize>` method. This method reads in all of the element cards
from the BDF file and sets up the equivalent TACS element objects necessary for analysis. The
class can be initialized through two ways: with or without ``elemCallBack`` function.

Initializing with elemCallBack
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`~tacs.pytacs.pyTACS` can be initialized py passing a user-defined ``elemCallBack`` function to
:meth:`pyTACS.initialize <tacs.pytacs.pyTACS.initialize>`, that will be used to setup the correct
TACS elements at runtime.

The ``elemCallBack`` function should have the following structure:

.. function:: elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs)

  User-defined function used by :class:`~tacs.pytacs.pyTACS` to setup a :ref:`Element<core/elements:Element classes>`
  for each element type in a given component (property) group.

  :param dvNum: The current counter which must be used when creating
           constitutive object with design
           variables.
  :type dvNum: int
  :param compID: The ID number used by TACS to reference this property group.
           Use kwargs['propID'] to get the corresponding NASTRAN property ID that
           is read in from the BDF.
  :type compID: int
  :param compDescript: The component description label read in from optional
           formatted comments in BDF file.
  :type compDescript: str
  :param elemDescripts: The name of the NASTRAN element cards belonging to this group
           (e.g. CQUAD4, CTRIA3, CTETRA, etc). This value will be a list since
           one component may contain multiple compatible element types.
           Example: ['CQUAD4', CTRIA3'].
  :type elemDescripts: list[str]
  :param globalDVs: Dictionary containing information about any
           global DVs that have been added.
  :type globalDVs: dict
  :returns: List containing as many TACS element
           objects as there are element types in `elemDescripts` (one for each).
  :rtype: list[:ref:`Element<core/elements:Element classes>`]

.. note::
  Some elements should not be setup through the ``elemCallBack`` function (RBE2, RBE3, CONM1, CONM2, etc.).
  These elements are setup automatically by :class:`~tacs.pytacs.pyTACS`. As a general rule of thumb,
  if the NASTRAN element card doesn't contain a property ID, it should not show up in the ``elemCallBack`` function.

Initializing without elemCallBack
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If NASTRAN material property card definitions for every element exist in the BDF file,
:class:`~tacs.pytacs.pyTACS` can be initialized without a :func:`~elemCallBack` function.
This can be done by calling :meth:`pyTACS.initialize <tacs.pytacs.pyTACS.initialize>` without any arguments.

Currently supported NASTRAN cards and their corresponding TACS-equivelant classes are listed below:

  Material cards:
    - MAT1 -> :class:`~tacs.constitutive.MaterialProperties`
    - MAT2 -> :class:`~tacs.constitutive.MaterialProperties`
    - MAT8 -> :class:`~tacs.constitutive.MaterialProperties`

  Property cards:
    - PSHELL -> :class:`~tacs.constitutive.IsoShellConstitutive`
    - PCOMP -> :class:`~tacs.constitutive.CompositeShellConstitutive`
    - PBAR, PROD -> :class:`~tacs.constitutive.BasicBeamConstitutive`
    - PSOLID -> :class:`~tacs.constitutive.SolidConstitutive`
    - PBUSH -> :class:`~tacs.constitutive.DOFSpringConstitutive`

  Elements cards:
    - CQUAD4, CQUADR -> :class:`~tacs.elements.Quad4Shell`
    - CQUAD9 -> :class:`~tacs.elements.Quad9Shell`
    - CTRIA3, CTRIAR -> :class:`~tacs.elements.Tri3Shell`
    - CBAR, CROD -> :class:`~tacs.elements.Beam2`
    - CHEXA -> :class:`~tacs.elements.Element3D` (:class:`~tacs.elements.LinearHexaBasis`, :class:`~tacs.elements.LinearElasticity3D`)
    - CTETRA -> :class:`~tacs.elements.Element3D` (:class:`~tacs.elements.LinearTetrahedralBasis`, :class:`~tacs.elements.LinearElasticity3D`)
    - RBE2 -> :class:`~tacs.elements.RBE2`
    - RBE3 -> :class:`~tacs.elements.RBE3`
    - CONM1 -> :class:`~tacs.elements.MassElement`
    - CONM2 -> :class:`~tacs.elements.MassElement`
    - CBUSH -> :class:`~tacs.elements.SpringElement`

  Design Variable cards:
    - DESVAR

.. note::
  Not every NASTRAN element card has one unique TACS element counterpart.
  For example a user may want a `CQUAD4` card interpreted as :class:`~tacs.elements.Quad4Shell` in the case of an elastic analysis,
  but might want it interpreted as a :class:`~tacs.elements.Quad4ThermalShell` in the case of a thermoelastic analysis.
  In the case where the default TACS elements above don't match the user's desired TACS element,
  the user is recommended to setup that element through an
  :ref:`elemCallBack procedure <pytacs/pytacs_module:Initializing with elemCallBack>` instead.

Tagging component groups in BDF
-------------------------------
Several pyTACS methods (:meth:`pyTACS.selectCompIDs <tacs.pytacs.pyTACS.selectCompIDs>`, :func:`~elemCallBack`, etc.)
allow for the use of user-defined component labels. These labels are read in through formatted comment in the BDF file.

There are currently two supported formats for labels: ICEM-format and FEMAP-format. These are described below.

ICEM component label format
^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the ICEM format, the elements are organized by property groups and
the component label is appended as a string comment above the first element in every property group.
The string comment should start with ``$       Shell element data for family`` followed by the component name.
No spaces are allowed within the component name. An example of this format is provided below:

.. code-block:: none

  $       Shell element data for family    WING_SPARS/LE_SPAR/SEG.00
  CQUADR         1       1    3600    3310    3797     731       1
  CQUADR         2       1     731    3797    3798     732       1
  CQUADR         3       1     732    3798    3799     733       1
  CQUADR         4       1     733    3799    3800     734       1
  CQUADR         5       1     734    3800    3801     735       1
  CQUADR         6       1     735    3801    3802     736       1
  $       Shell element data for family    WING_SPARS/LE_SPAR/SEG.01
  CQUADR        97       2    3262    3882     782    3601       2
  CQUADR        98       2     782    3882    3881     781       2
  CQUADR        99       2    3875    3888    3885    3874       2
  CQUADR       100       2    3885    3888    3887    3884       2
  CQUADR       101       2    3892    3899    3896    3891       2
  CQUADR       102       2    3896    3899    3898    3895       2

FEMAP component label format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the FEMAP format, the component label is appended as a string comment above each property card.
The string comment should start with ``$ Femap Property i :`` (where ``i`` is replaced with the property number the label is referencing)
followed by the component name. No spaces are allowed within the component name. An example of this format is provided below:

.. code-block:: none

  $ Femap Property 1 : Plate
  PSHELL         1       1      .1       1               1              0.
  $ Femap Property 2 : Stiffener
  PSHELL         2       1      .1       1               1              0.

API Reference
-------------
.. autoclass:: tacs.pytacs.pyTACS
  :members:
  :inherited-members:
