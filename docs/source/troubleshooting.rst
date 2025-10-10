Troubleshooting Guide
====================

This guide helps you resolve common issues when using TACS. If you don't find your issue here, check the `GitHub Issues <https://github.com/smdogroup/tacs/issues>`_ for community support.

Installation Issues
-------------------

**ImportError: No module named 'tacs'**

*Cause*: TACS is not installed or not in the Python path.

*Solutions*:

- Ensure conda environment is activated: ``conda activate tacs``
- Reinstall TACS: ``mamba install -c conda-forge -c smdogroup tacs``
- For source installation, ensure Python interface is built: ``pip install -e .[all]`` or ``make interface``

**ImportError: libtacs.so not found**

*Cause*: TACS C++ library is not in the library path.

*Solutions*:

.. code-block:: bash

   # Add TACS lib directory to LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/path/to/tacs/lib:$LD_LIBRARY_PATH
   
   # Or add to your .bashrc for permanent fix
   echo 'export LD_LIBRARY_PATH=/path/to/tacs/lib:$LD_LIBRARY_PATH' >> ~/.bashrc

**MPI compilation errors**

*Cause*: MPI compiler not found or misconfigured.

*Solutions*:

- Ensure MPI is installed: ``which mpicc``
- Check MPI compiler in Makefile.in: ``CXX = mpicxx``
- For Intel MPI: ``CXX = icpc -lmpi``

Analysis Issues
---------------

**Analysis fails to converge**

*Cause*: Poor conditioning, incorrect boundary conditions, or excessive loads.

*Solutions*:

- Check boundary conditions are properly applied
- Verify loads are reasonable for the structure
- Ensure material properties are correct
- Try reducing load magnitude
- Check for rigid body modes

**Element callback errors**

*Cause*: Mismatch between BDF element types and callback function.

*Solutions*:

- Verify element types in BDF file match callback function
- Add error handling for unexpected element types:

  .. code-block:: python

     for descript in elemDescripts:
         if descript == 'CQUAD4':
             elem = elements.Quad4Shell(transform, con)
         else:
             raise ValueError(f"Unexpected element type: {descript}")

Postprocessing Issues
--------------------

**f5tovtk/f5totec not found**

*Cause*: Postprocessing utilities not installed or not in PATH.

*Solutions*:

.. code-block:: bash

   # Add to PATH
   export PATH=/path/to/tacs/extern/f5tovtk:$PATH
   export PATH=/path/to/tacs/extern/f5totec:$PATH
   
   # Or use full path
   /path/to/tacs/extern/f5tovtk/f5tovtk solution.f5

**Conversion fails with large files**

*Cause*: Insufficient disk space or memory.

*Solutions*:

- Use essential output flags only
- Convert files individually

Platform-Specific Issues
------------------------

**Windows/WSL Issues**

*Cause*: Path or permission issues in WSL.

*Solutions*:

- Use forward slashes in paths
- Ensure proper file permissions
- Use WSL2 for better performance
- Consider Docker alternative

**macOS Issues**

*Cause*: Library path or compiler issues.

*Solutions*:

- Use conda installation for easier setup
- Ensure Xcode command line tools installed
- Check library paths: ``otool -L libtacs.so``

Debugging Tips
--------------

**Check TACS Version**

.. code-block:: python

   import tacs
   print("TACS version:", tacs.__version__)

**Verify Installation**

.. code-block:: python

   import tacs
   from tacs import functions, constitutive, elements, pyTACS
   print("All modules imported successfully")

Getting Help
------------

**Community Resources:**

- `GitHub Issues <https://github.com/smdogroup/tacs/issues>`_: Report bugs and request features

**When Reporting Issues:**

Include the following information:
- TACS version: ``tacs.__version__``
- Python version: ``python --version``
- Operating system and version
- Complete error message and traceback
- Minimal code example that reproduces the issue
- Expected vs. actual behavior

**Useful Commands:**

.. code-block:: bash

   # Check TACS installation
   python -c "import tacs; print(tacs.__version__)"
   
   # Check MPI
   mpirun --version
   
   # Check libraries
   ldd libtacs.so  # Linux
   otool -L libtacs.so  # macOS
   
   # Check environment
   echo $LD_LIBRARY_PATH
   echo $PATH
