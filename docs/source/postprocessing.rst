Postprocessing TACS Solutions
=============================

TACS provides comprehensive postprocessing capabilities through its native f5 file format and conversion utilities. This section explains how to generate f5 files from TACS analyses and convert them to standard visualization formats.

F5 File Format
--------------

TACS uses the f5 file format (based on HDF5) to store solution data in a parallel, binary format. F5 files contain both continuous (nodal) and element-wise data, making them suitable for detailed postprocessing and visualization.

Creating F5 Files
~~~~~~~~~~~~~~~~~

To create an f5 file from a TACS analysis, use the ``TACSToFH5`` class:

.. code-block:: cpp

   #include "TACSToFH5.h"
   
   // Create an TACSToFH5 object for writing output to files
   ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
   int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                     TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                     TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
   TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
   f5->incref();
   f5->writeToFile("solution.f5");
   f5->decref();

Output Flags
~~~~~~~~~~~~

The following output flags control what data is written to the f5 file:

- ``TACS_OUTPUT_CONNECTIVITY``: Element connectivity information
- ``TACS_OUTPUT_NODES``: Nodal coordinates (X, Y, Z)
- ``TACS_OUTPUT_DISPLACEMENTS``: Nodal displacements and rotations
- ``TACS_OUTPUT_STRAINS``: Element strains
- ``TACS_OUTPUT_STRESSES``: Element stresses
- ``TACS_OUTPUT_EXTRAS``: Additional quantities (failure indices, design variables)
- ``TACS_OUTPUT_LOADS``: Applied loads
- ``TACS_OUTPUT_COORDINATE_FRAME``: Element coordinate frames

Converting F5 Files
-------------------

TACS provides two utilities for converting f5 files to standard visualization formats:

f5totec: Convert to Tecplot Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``f5totec`` utility converts f5 files to Tecplot format (.plt files):

.. code-block:: bash

   # Basic conversion
   f5totec solution.f5
   
   # Convert with strands (for time-dependent data)
   f5totec --use_strands solution.f5

This creates a ``solution.plt`` file that can be opened in Tecplot.

f5tovtk: Convert to VTK Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``f5tovtk`` utility converts f5 files to VTK format (.vtk files) for use with ParaView:

.. code-block:: bash

   f5tovtk solution.f5

This creates a ``solution.vtk`` file that can be opened in ParaView.

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

Example Workflow
----------------

Here's a complete example of generating and visualizing TACS results:

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

The resulting files (``wing_analysis.plt`` and ``wing_analysis.vtk``) can then be opened in Tecplot or ParaView for visualization and further analysis.
