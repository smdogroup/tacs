Interfaces
**********
TACS is implemented in C++, so the interface through C++ contains all publicly accessible class member functions.
A Python-level interface is created by wrapping the most important classes and functions into Python using Cython.
There are currently two supported approaches for interfacing with this functionality in Python: direct and pyTACS.
These two methods are documented below.

.. toctree::
  :maxdepth: 1

  core/TACS
  pytacs/pytacs