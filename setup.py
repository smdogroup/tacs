import os
from subprocess import check_output
import sys

# Numpy/mpi4py must be installed prior to installing TACS
import numpy
import mpi4py

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize

# Convert from local to absolute directories
def get_global_dir(files):
    tacs_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tacs_root, f))
    return new

def get_mpi_flags():
    # Split the output from the mpicxx command
    args = check_output(['mpicxx', '-show']).decode('utf-8').split()

    # Determine whether the output is an include/link/lib command
    inc_dirs, lib_dirs, libs = [], [], []
    for flag in args:
        if flag[:2] == '-I':
            inc_dirs.append(flag[2:])
        elif flag[:2] == '-L':
            lib_dirs.append(flag[2:])
        elif flag[:2] == '-l':
            libs.append(flag[2:])

    return inc_dirs, lib_dirs, libs

inc_dirs, lib_dirs, libs = get_mpi_flags()

# Add tacs-dev/lib as a runtime directory
runtime_lib_dirs = get_global_dir(['lib'])

# Relative paths for the include/library directories
rel_inc_dirs = ['src', 'src/bpmat', 'src/elements', 
                'src/constitutive', 'src/functions', 'src/io']
rel_lib_dirs = ['lib']
libs.extend(['tacs'])

# Convert from relative to absolute directories
inc_dirs.extend(get_global_dir(rel_inc_dirs))
lib_dirs.extend(get_global_dir(rel_lib_dirs))

# This should be made more general so that you can specify alternate
# locations for the installation of AMD/METIS
default_ext_inc = ['extern/AMD/Include', 
                   'extern/UFconfig', 
                   'extern/metis/include']
inc_dirs.extend(get_global_dir(default_ext_inc))

# Add the numpy/mpi4py directories
inc_dirs.extend([numpy.get_include(), mpi4py.get_include()])

exts = []
for mod in ['TACS', 'elements', 'constitutive', 'functions']:
    exts.append(Ext('tacs.%s'%(mod), sources=['tacs/%s.pyx'%(mod)],
                    include_dirs=inc_dirs, libraries=libs, 
                    library_dirs=lib_dirs, runtime_library_dirs=runtime_lib_dirs))
for e in exts:
    e.cython_directives = {"embedsignature": True,
                           "binding":True}
    
setup(name='tacs',
      version=0.1,
      description='Parallel finite-element analysis package',
      author='Graeme J. Kennedy',
      author_email='graeme.kennedy@ae.gatech.edu',
      ext_modules=cythonize(exts, include_path=inc_dirs))