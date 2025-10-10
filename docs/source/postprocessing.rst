Postprocessing TACS Solutions
=============================

TACS provides comprehensive postprocessing capabilities through its native f5 file format and conversion utilities. This section explains how to generate f5 files from TACS analyses and convert them to standard visualization formats for detailed analysis and presentation.

Overview
--------

The TACS postprocessing workflow consists of three main steps:

1. **Generate f5 files** from your TACS analysis containing solution data
2. **Convert f5 files** to standard visualization formats (Tecplot or ParaView)
3. **Visualize and analyze** results using your preferred visualization software

F5 File Format
--------------

TACS uses the f5 file format (based on HDF5) to store solution data in a parallel, binary format. F5 files contain both continuous (nodal) and element-wise data, making them suitable for detailed postprocessing and visualization.

Creating F5 Files
~~~~~~~~~~~~~~~~~

**Method 1: Using the Direct (C++) Interface**

To create an f5 file from a TACS analysis, use the ``TACSToFH5`` class. This should be done after solving your analysis but before cleaning up the assembler object.

.. code-block:: cpp

   #include "TACSToFH5.h"

   // After solving your analysis
   // Create an TACSToFH5 object for writing output to files
   ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
   int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                     TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                     TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
   TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
   f5->incref();
   f5->writeToFile("solution.f5");
   f5->decref();

.. tip::
   The ``ElementType`` should match the type of elements in your model. Common types include:

   - ``TACS_BEAM_OR_SHELL_ELEMENT`` for shell and beam elements
   - ``TACS_SOLID_ELEMENT`` for 3D solid elements
   - ``TACS_PLANE_STRESS_ELEMENT`` for 2D plane stress elements

**Method 2: Using the pyTACS (Python) interface**

To create an f5 file from a TACS analysis, use the `writeSolution` method of any :doc:`problem class <pytacs/problems>`
after solving the problem.

.. code-block:: python

   from tacs import pyTACS, TACS
   bdfFile = "input.bdf"
   options = {# Specify the element type to output in f5 file
              # If not specified, pyTACS will choose this option automatically
              # based on first element type in model
              "outputElement": TACS.BEAM_OR_SHELL_ELEMENT,
              # Output flags
              "writeConnectivity": True,
              "writeNodes": True,
              "writeDisplacements": True,
              "writeStrains": True,
              "writeStresses": True,
              "writeExtras": True}
   FEAAssembler = pyTACS(bdfFile, options=options)
   FEAAssembler.initialize()
   # Creata a pytacs problem
   staticProb = FEAAssembler.createStaticProblem("grav")
   staticProb.addInertialLoad([0.0, 0.0, -9.81])
   # Solve
   staticProb.solve()
   # Write the solution to an f5 file
   FEAAssembler.writeSolution()

Output Flags
~~~~~~~~~~~~

The following output flags control what data is written to the f5 file. The table shows both the C++ interface flags and their equivalent pyTACS options. You can combine multiple C++ flags using the bitwise OR operator (``|``), while :ref:`pyTACS options<pytacs/pytacs_module:Options>` are set as boolean values in the options dictionary when the :class:`~tacs.pytacs.pyTACS` object is created:

.. list-table:: Output Flags
   :widths: 25 25 50
   :header-rows: 1

   * - Direct (C++) Flag
     - pyTACS (Python) Option
     - Description
   * - ``TACS_OUTPUT_CONNECTIVITY``
     - ``writeConnectivity``
     - Element connectivity information (required for visualization)
   * - ``TACS_OUTPUT_NODES``
     - ``writeNodes``
     - Nodal coordinates (X, Y, Z) - essential for geometry visualization
   * - ``TACS_OUTPUT_DISPLACEMENTS``
     - ``writeDisplacements``
     - Nodal displacements and rotations - needed for deformed shape visualization
   * - ``TACS_OUTPUT_STRAINS``
     - ``writeStrains``
     - Element strains - useful for strain analysis and contour plots
   * - ``TACS_OUTPUT_STRESSES``
     - ``writeStresses``
     - Element stresses - essential for stress analysis and failure assessment
   * - ``TACS_OUTPUT_EXTRAS``
     - ``writeExtras``
     - Additional quantities (failure indices, design variables) - useful for optimization
   * - ``TACS_OUTPUT_LOADS``
     - ``writeLoads``
     - Applied loads - helpful for load verification and visualization
   * - ``TACS_OUTPUT_COORDINATE_FRAME``
     - ``writeCoordinateFrame``
     - Element coordinate frames - useful for composite material analysis

.. note::
   For basic visualization, you typically need at least ``TACS_OUTPUT_CONNECTIVITY``, ``TACS_OUTPUT_NODES``, and ``TACS_OUTPUT_DISPLACEMENTS``. Add other flags based on your analysis requirements.

Converting F5 Files
-------------------

TACS provides two utilities for converting f5 files to standard visualization formats. These utilities are typically located in the ``extern/`` directory of your TACS installation.

f5totec: Convert to Tecplot Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``f5totec`` utility converts f5 files to Tecplot format (.plt files):

.. code-block:: bash

   # Basic conversion
   f5totec solution.f5

This creates a ``solution.plt`` file that can be opened in Tecplot.

f5tovtk: Convert to VTK Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``f5tovtk`` utility converts f5 files to VTK format (.vtk files) for use with ParaView:

.. code-block:: bash

   # Basic conversion
   f5tovtk solution.f5

This creates a ``solution.vtk`` file that can be opened in ParaView.

.. note::
   When a node is used by multiple elements, each element may have a different value for variables such as stress,
   strain, failure criteria, and design variables at that node. f5totec produces a single value for each node by
   averaging the values from each element. This can lead to unrealistic values of these variables in certain situations
   (e.g design variable values at the boundaries between different components and stress/strain/failure criteria values
   at points where shell elements meet at very different orientations.

**Troubleshooting Conversion Issues:**

- Ensure the f5 file was generated successfully and contains the expected data
- Check that the conversion utilities are compiled and accessible in your PATH

Output Variables by Element Type
--------------------------------

The following tables describe the output variables available for each element type in TACS.

.. tip::
   For elements that have a local coordinate system (e.g shells and beams), the stress and strain outputs are in the local coordinate system. For example, ``ex0``/``sx0`` are in the direction of the local reference axis, ``ey0``/``sy0`` are in the direction of the second reference frame vector. Other variables, such as forces ``fx``, ``fy``, ``fz`` are given in the global reference frame.

Beam/Shell Elements (TACS_BEAM_OR_SHELL_ELEMENT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Beam/Shell Element Output Variables
   :widths: 20 20 60
   :header-rows: 1

   * - Category
     - Variable
     - Description
   * - Displacements
     - u, v, w
     - Translational displacements
   * -
     - rotx, roty, rotz
     - Rotational displacements
   * - Strains
     - ex0, ey0, exy0
     - Membrane strains
   * -
     - ex1, ey1, exy1
     - Bending strains
   * -
     - eyz0, exz0
     - Transverse shear strains
   * -
     - erot
     - Drilling strain
   * - Stresses
     - sx0, sy0, sxy0
     - Membrane stress resultants
   * -
     - sx1, sy1, sxy1
     - Bending stress resultants
   * -
     - syz0, sxz0
     - Transverse shear stress resultants
   * -
     - srot
     - Drilling stress resultant
   * - Extras
     - failure0-failure6
     - Failure indices for different failure criteria
   * -
     - dv1-dv7
     - Design variables
   * - Loads
     - fx, fy, fz
     - Applied forces
   * -
     - mx, my, mz
     - Applied moments
   * - Coordinate Frame
     - t0x, t0y, t0z
     - First element reference frame vector (i.e. reference axis) components
   * -
     - t1x, t1y, t1z
     - Second element reference frame vector components
   * -
     - t2x, t2y, t2z
     - Third element reference frame vector (i.e. normal vector) components

Solid Elements (TACS_SOLID_ELEMENT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Solid Element Output Variables
   :widths: 20 20 60
   :header-rows: 1

   * - Category
     - Variable
     - Description
   * - Displacements
     - u, v, w
     - Translational displacements
   * - Strains
     - exx, eyy, ezz
     - Normal strains
   * -
     - gyz, gxz, gxy
     - Shear strains
   * - Stresses
     - sxx, syy, szz
     - Normal stresses
   * -
     - syz, sxz, sxy
     - Shear stresses
   * - Extras
     - failure
     - Failure index
   * -
     - dv1, dv2, dv3
     - Design variables
   * - Loads
     - fx, fy, fz
     - Applied forces

Plane Stress Elements (TACS_PLANE_STRESS_ELEMENT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Plane Stress Element Output Variables
   :widths: 20 20 60
   :header-rows: 1

   * - Category
     - Variable
     - Description
   * - Displacements
     - u, v
     - In-plane displacements
   * - Strains
     - exx, eyy, gxy
     - In-plane strains
   * - Stresses
     - sxx, syy, sxy
     - In-plane stresses
   * - Extras
     - failure
     - Failure index
   * -
     - dv1, dv2, dv3
     - Design variables
   * - Loads
     - fx, fy
     - Applied forces

Scalar Elements (TACS_SCALAR_2D_ELEMENT, TACS_SCALAR_3D_ELEMENT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Scalar Element Output Variables
   :widths: 20 20 60
   :header-rows: 1

   * - Category
     - Variable
     - Description
   * - Displacements
     - u
     - Scalar displacement
   * - Strains
     - ux, uy (2D) / ux, uy, uz (3D)
     - Gradient components
   * - Stresses
     - sx, sy (2D) / sx, sy, sz (3D)
     - Flux components
   * - Extras
     - failure
     - Failure index
   * -
     - dv1, dv2, dv3
     - Design variables
   * - Loads
     - f
     - Applied load

PCM Elements (TACS_PCM_ELEMENT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: PCM Element Output Variables
   :widths: 20 20 60
   :header-rows: 1

   * - Category
     - Variable
     - Description
   * - Displacements
     - dT
     - Temperature change
   * - Strains
     - gradx, grady
     - Temperature gradient components
   * - Stresses
     - fluxx, fluxy
     - Heat flux components
   * - Extras
     - rho
     - Density
   * -
     - dv1, dv2, dv3
     - Design variables
   * -
     - phase
     - Phase field
   * - Loads
     - Q
     - Applied heat source

Visualization Tips
------------------

1. **Element-wise vs. Nodal Data**: F5 files contain both element-wise and nodal data. The conversion utilities automatically perform nodal averaging for element-wise quantities.

2. **Higher-order Elements**: Higher-order elements are split into multiple lower order element for visualization (e.g., each quadratic triangle becomes 3 linear triangles).

3. **Component Separation**: In Tecplot, each component in the model can be written as a separate zone in the output files, making it easy to visualize different parts of the structure.

Visualizing Deformed Surfaces
-----------------------------

One of the most common postprocessing tasks is visualizing the deformed shape of structures. TACS provides both nodal coordinates (X, Y, Z) and displacements (u, v, w) that can be used to create deformed surface visualizations.

Creating Deformed Geometry in Tecplot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Tecplot, you can visualize deformed surfaces by creating new variables that represent the deformed coordinates:

1. **Open the converted .plt file** in Tecplot
2. **Create new variables** for deformed coordinates:

   - Go to ``Data > Alter > Specify Equations``
   - Under the Equations box enter:

     ::

      {XDEF} = {X} + {u}
      {YDEF} = {Y} + {v}
      {ZDEF} = {Z} + {w}

3. **Create the deformed plot**:

   - Go to ``Plot > Assign XYZ...``
   - Set ``X``, ``Y``, ``Z`` to ``XDEF``, ``YDEF``, ``ZDEF``
   - Choose appropriate surface rendering (``Surface``, ``Mesh``, or ``Contour``)

Tecplot Macro
^^^^^^^^^^^^^
You can also automate this process using a Tecplot macro. Place this code in your `tecplot.mcr` file to make it available as a quick macro in the Tecplot GUI:

.. code-block::

    $!MACROFUNCTION NAME = "TACS - Apply Deformations"
    $!PromptForTextString |DefFactor|
      Instructions = "Enter scaling factor for deformations"

    $!AlterData
      Equation = "{X}={X}+|DefFactor|*{u}"
    $!AlterData
      Equation = "{Y}={Y}+|DefFactor|*{v}"
    $!AlterData
      Equation = "{Z}={Z}+|DefFactor|*{w}"
    $!ENDMACROFUNCTION

Note that this macro modifies the original X, Y, Z variables, rather than creating new variables for the deformed coordinates as is done above.

Creating Deformed Geometry in ParaView
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParaView provides several methods to visualize deformed surfaces:

1. **Open the converted .vtk file** in ParaView

2. **Add Calculator filter**:
   - Select the dataset
   - Go to ``Filters > Alphabetical > Calculator``

3. **Create deformed coordinates**:

   - Set ``Result Array Name`` to ``def_vec``
   - Set ``Function`` to ``u*iHat + v*jHat + w*kHat``
   - Click ``Apply``

4. **Add Warp By Vector filter**:

   - Select the dataset
   - Go to ``Filters > Alphabetical > Warp By Vector``

5. **Configure the warp**:

   - Set ``Vector`` to ``def_vec``
   - Adjust ``Scale Factor`` to control deformation magnification
   - Click ``Apply``
