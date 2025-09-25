TACS Interfaces
===============

TACS provides multiple interfaces to suit different user needs and workflows. The core functionality is implemented in C++, with Python interfaces created using Cython for easier integration and analysis.

Interface Overview
------------------

TACS offers several interface options:

.. list-table:: TACS Interface Options
   :widths: 20 25 55
   :header-rows: 1

   * - Interface
     - Best For
     - Description
   * - **C++ Direct**
     - High-performance applications
     - Full access to all TACS functionality with maximum performance
   * - **pyTACS**
     - General Python users
     - High-level Python interface for structural analysis and optimization
   * - **MPhys**
     - Multidisciplinary optimization
     - Integration with OpenMDAO for MDO workflows
   * - **caps2tacs**
     - CAD integration
     - Direct integration with CAPS for geometry-driven analysis
   * - **MACH**
     - Library for wrapping pyTACS interface into MDOLab's MACH framework
     - Useful for structural and aerostructural optimization with geometric design variables

Choosing the Right Interface
----------------------------

**For New Users:**
Start with **pyTACS** - it provides the most intuitive interface for structural analysis with comprehensive documentation and examples.

**For Optimization:**
Use **MPhys** if you need to integrate TACS with other disciplines in an OpenMDAO-based optimization framework.

**For CAD Integration:**
Use **caps2tacs** if you want to perform analysis directly from CAD geometry without manual mesh generation.

**For MDOLab MACH Integration:**
Use **MACH** if you need to interface TACS with MDOLab MACH library codes, particularly for aerostructural coupling or other multidisciplinary analysis workflows.

**For Performance-Critical Applications:**
Use the **C++ Direct** interface for maximum performance and full control over the analysis process.

Interface Documentation
-----------------------

.. toctree::
  :maxdepth: 1

  core/TACS
  pytacs/pytacs
  mphys/mphys
  caps2tacs/caps2tacs
  mach/mach
