"""
TACS is a parallel finite-element package for analysis and
gradient-based design optimization.
"""

import os

__all__ = []


def get_cython_include():
    """
    Get the include directory for the Cython .pxd files in TACS
    """
    return [os.path.abspath(os.path.dirname(__file__))]


def get_include():
    """
    Get the include directory for the Cython .pxd files in TACS
    """
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_inc_dirs = [
        "src",
        "src/bpmat",
        "src/elements",
        "src/elements/dynamics",
        "src/elements/basis",
        "src/constitutive",
        "src/functions",
        "src/io",
        "extern/AMD/Include",
        "extern/UFconfig",
        "extern/metis/include",
    ]

    inc_dirs = []
    for path in rel_inc_dirs:
        inc_dirs.append(os.path.join(root_path, path))

    return inc_dirs


def get_libraries():
    """
    Get the library directories
    """
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_lib_dirs = ["lib"]
    libs = ["tacs"]
    lib_dirs = []
    for path in rel_lib_dirs:
        lib_dirs.append(os.path.join(root_path, path))

    return lib_dirs, libs


# Import pytacs modules
from . import pytacs
from .pytacs import pyTACS
from . import problems

__all__.extend(["pytacs", "pyTACS", "problems"])
