Interfaces
**********
TACS is implemented in C++, so the interface through C++ contains all publicly accessible class member functions.
A Python-level interface is created by wrapping the most important classes and functions into Python using Cython.
There are currently three supported approaches for interfacing with this functionality in Python: direct, pyTACS, and MPhys.
These methods are documented below.

.. toctree::
  :maxdepth: 1

  core/TACS
  pytacs/pytacs
  mphys/mphys
  caps2tacs/caps2tacs
