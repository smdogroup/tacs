Postprocessing TACS Solutions
=============================

TACS provides comprehensive postprocessing capabilities through its native f5 file format and conversion utilities. This section explains how to generate f5 files from TACS analyses and convert them to standard visualization formats for detailed analysis and presentation.

Overview
--------

The TACS postprocessing workflow consists of three main steps:

1. **Generate f5 files** from your TACS analysis containing solution data
2. **Convert f5 files** to standard visualization formats (Tecplot or ParaView)
3. **Visualize and analyze** results using your preferred visualization software

The f5 format stores both nodal (continuous) and element-wise data, making it suitable for comprehensive structural analysis visualization.

F5 File Format
--------------

TACS uses the f5 file format (based on HDF5) to store solution data in a parallel, binary format. F5 files contain both continuous (nodal) and element-wise data, making them suitable for detailed postprocessing and visualization.

.. important::
   F5 files are binary and platform-dependent. They should be generated and converted on the same system architecture for best compatibility.

**Key Features of F5 Files:**
- **Parallel I/O**: Efficient for large-scale parallel simulations
- **Binary format**: Fast read/write operations and compact storage
- **Hierarchical structure**: Organized data storage with metadata

Creating F5 Files
~~~~~~~~~~~~~~~~~

To create an f5 file from a TACS analysis, use the ``TACSToFH5`` class. This should be done after solving your analysis but before cleaning up the assembler object.

.. code-block:: cpp

   #include "TACSToFH5.h"
   
   // After solving your analysis (e.g., assembler->solve())
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

Output Flags
~~~~~~~~~~~~

The following output flags control what data is written to the f5 file. You can combine multiple flags using the bitwise OR operator (``|``):

.. list-table:: Output Flags
   :widths: 30 70
   :header-rows: 1

   * - Flag
     - Description
   * - ``TACS_OUTPUT_CONNECTIVITY``
     - Element connectivity information (required for visualization)
   * - ``TACS_OUTPUT_NODES``
     - Nodal coordinates (X, Y, Z) - essential for geometry visualization
   * - ``TACS_OUTPUT_DISPLACEMENTS``
     - Nodal displacements and rotations - needed for deformed shape visualization
   * - ``TACS_OUTPUT_STRAINS``
     - Element strains - useful for strain analysis and contour plots
   * - ``TACS_OUTPUT_STRESSES``
     - Element stresses - essential for stress analysis and failure assessment
   * - ``TACS_OUTPUT_EXTRAS``
     - Additional quantities (failure indices, design variables) - useful for optimization
   * - ``TACS_OUTPUT_LOADS``
     - Applied loads - helpful for load verification and visualization
   * - ``TACS_OUTPUT_COORDINATE_FRAME``
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
   
   # Convert with strands (for time-dependent data)
   f5totec --use_strands solution.f5
   
   # Specify output filename
   f5totec solution.f5 -o wing_analysis.plt

This creates a ``solution.plt`` file that can be opened in Tecplot.

.. tip::
   The ``--use_strands`` option is particularly useful for transient analyses or optimization histories where you want to animate the results over time.

f5tovtk: Convert to VTK Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``f5tovtk`` utility converts f5 files to VTK format (.vtk files) for use with ParaView:

.. code-block:: bash

   # Basic conversion
   f5tovtk solution.f5
   
   # Specify output filename
   f5tovtk solution.f5 -o wing_analysis.vtk

This creates a ``solution.vtk`` file that can be opened in ParaView.

.. note::
   Both utilities automatically handle the conversion of element-wise data to nodal data through averaging, making the results suitable for smooth visualization.

**Troubleshooting Conversion Issues:**
- Ensure the f5 file was generated successfully and contains the expected data
- Check that the conversion utilities are compiled and accessible in your PATH
- For large files, conversion may take several minutes - this is normal

Output Variables by Element Type
--------------------------------

The following tables describe the output variables available for each element type in TACS.

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
     - Rotational strain
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
     - Rotational stress
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
     - First tangent vector components
   * - 
     - t1x, t1y, t1z
     - Second tangent vector components
   * - 
     - t2x, t2y, t2z
     - Normal vector components

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

2. **Higher-order Elements**: Higher-order elements are automatically converted to basic element types for visualization (e.g., quadratic triangles become linear triangles).

3. **Component Separation**: In Tecplot, each component in the model can be written as a separate zone in the output files, making it easy to visualize different parts of the structure.

4. **Time-dependent Data**: Use the ``--use_strands`` option with ``f5totec`` for time-dependent analyses to create animated visualizations.

5. **Large Models**: For very large models, consider using only the essential output flags to reduce file size and processing time.

Visualizing Deformed Surfaces
-----------------------------

One of the most common postprocessing tasks is visualizing the deformed shape of structures. TACS provides both nodal coordinates (X, Y, Z) and displacements (u, v, w) that can be used to create deformed surface visualizations.

Creating Deformed Geometry in Tecplot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Tecplot, you can visualize deformed surfaces by creating new variables that represent the deformed coordinates:

1. **Open the converted .plt file** in Tecplot
2. **Create new variables** for deformed coordinates:
   - Go to ``Data > Alter > Specify Equations``
   - Create new variables:
     - ``{XDEF} = {X} + {u}`` (deformed X coordinate)
     - ``{YDEF} = {Y} + {v}`` (deformed Y coordinate)  
     - ``{ZDEF} = {Z} + {w}`` (deformed Z coordinate)
3. **Create the deformed plot**:
   - Go to ``Plot > 3D``
   - Set ``X``, ``Y``, ``Z`` to ``XDEF``, ``YDEF``, ``ZDEF``
   - Choose appropriate surface rendering (``Surface``, ``Mesh``, or ``Contour``)
4. **Add displacement magnitude** for color coding:
   - Create variable: ``{DISP_MAG} = sqrt({u}*{u} + {v}*{v} + {w}*{w})``
   - Set ``Contour`` to ``DISP_MAG`` to color-code by displacement magnitude

**Example Tecplot equations:**
::

   {XDEF} = {X} + {u}
   {YDEF} = {Y} + {v}
   {ZDEF} = {Z} + {w}
   {DISP_MAG} = sqrt({u}*{u} + {v}*{v} + {w}*{w})

Creating Deformed Geometry in ParaView
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParaView provides several methods to visualize deformed surfaces:

**Method 1: Using the Calculator Filter**
1. **Open the converted .vtk file** in ParaView
2. **Add Calculator filter**:
   - Select the dataset
   - Go to ``Filters > Alphabetical > Calculator``
3. **Create deformed coordinates**:
   - Set ``Result Array Name`` to ``XDEF``
   - Set ``Function`` to ``X + u``
   - Click ``Apply``
   - Repeat for Y and Z coordinates:
     - ``YDEF = Y + v``
     - ``ZDEF = Z + w``
4. **Create displacement magnitude**:
   - Add another Calculator filter
   - Set ``Function`` to ``sqrt(u*u + v*v + w*w)``
   - Set ``Result Array Name`` to ``DISP_MAG``
5. **Visualize the deformed surface**:
   - Set ``Representation`` to ``Surface`` or ``Surface With Edges``
   - Set ``Coloring`` to ``DISP_MAG`` for displacement-based coloring

**Method 2: Using the Warp By Vector Filter**
1. **Open the .vtk file** in ParaView
2. **Add Warp By Vector filter**:
   - Select the dataset
   - Go to ``Filters > Alphabetical > Warp By Vector``
3. **Configure the warp**:
   - Set ``Vector`` to ``[u, v, w]`` (displacement components)
   - Adjust ``Scale Factor`` to control deformation magnification
   - Click ``Apply``
4. **Set visualization properties**:
   - Choose appropriate representation
   - Color by displacement magnitude or stress/strain variables

**Method 3: Using the Transform Filter**
1. **Add Transform filter**:
   - Go to ``Filters > Alphabetical > Transform``
2. **Configure translation**:
   - Set ``Translation`` to ``[u, v, w]``
   - Enable ``Transform All Input Vectors``
   - Click ``Apply``

Displacement Scaling and Animation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For better visualization, consider these techniques:

**Displacement Scaling:**
- **Tecplot**: Multiply displacement components by a scale factor in the equations (e.g., ``XDEF = X + 10*u``)
- **ParaView**: Use the ``Scale Factor`` in Warp By Vector filter or multiply in Calculator expressions

**Animation for Time-Dependent Results:**
- **Tecplot**: Use the ``Animation`` panel to cycle through time steps
- **ParaView**: Use the ``Animation View`` to create smooth animations of the deformation

**Combined Visualization:**
- Show both undeformed (wireframe) and deformed (solid) geometry
- Use different colors for different displacement components
- Overlay stress/strain contours on the deformed geometry

Example Workflow
----------------

Here's a complete example of generating and visualizing TACS results with deformed geometry:

.. code-block:: cpp

   // 1. Create f5 file from TACS analysis
   ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
   int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                     TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                     TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
   TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
   f5->writeToFile("wing_analysis.f5");
   f5->decref();

.. code-block:: bash

   # 2. Convert to Tecplot format
   f5totec wing_analysis.f5
   
   # 3. Convert to VTK format for ParaView
   f5tovtk wing_analysis.f5

**Tecplot Visualization Steps:**
1. Open ``wing_analysis.plt``
2. Create 10x deformed coordinates: ``{XDEF} = {X} + 10*{u}``, ``{YDEF} = {Y} + 10*{v}``, ``{ZDEF} = {Z} + 10*{w}``
3. Create displacement magnitude: ``{DISP_MAG} = sqrt({u}*{u} + {v}*{v} + {w}*{w})``
4. Plot using ``XDEF``, ``YDEF``, ``ZDEF`` with ``DISP_MAG`` coloring

**ParaView Visualization Steps:**
1. Open ``wing_analysis.vtk``
2. Add ``Warp By Vector`` filter with vector ``[u, v, w]``
3. Set scale factor to 10 for better visibility
4. Color by displacement magnitude or stress components

The resulting visualizations will show the deformed structure with appropriate scaling and coloring for effective analysis.
