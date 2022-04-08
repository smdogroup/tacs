import os
from subprocess import check_output
import sys

# Numpy/mpi4py must be installed prior to installing TACS
import numpy
import mpi4py

# Import distutils
from setuptools import setup, find_packages
from distutils.core import Extension as Ext
from Cython.Build import cythonize
from Cython.Compiler import Options

Options.docstrings = True

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
rel_inc_dirs = ['src', 'src/bpmat', 'src/elements', 'src/elements/dynamics',
                'src/elements/shell', 'src/elements/basis',
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
    e.cython_directives = {'embedsignature': True,
                           'binding': True}

tacs_root = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(tacs_root, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='tacs',
      version=0.1,
      description='Parallel finite-element analysis package',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Graeme J. Kennedy',
      author_email='graeme.kennedy@ae.gatech.edu',
      install_requires=[
          # Make sure the user's numpy version are at least more current than the build versions
          f'numpy>={numpy.__version__}',
          'mpi4py>=3.0.2',
          'scipy>=1.2.1',
          'pynastran>=1.3.3',
          # This package isn't actually a requirement of tacs, but of pynastran
          # Right now 2.0.0 breaks pynastran and the pip package hasn't been updated to fix this
          'nptyping<2.0.0',
      ],
      extras_require={
        'testing': ['testflo'],
        'docs': ['sphinx', 'breathe', 'sphinxcontrib-programoutput'],
        'mphys': ['mphys', 'openmdao'],
      },
      packages=find_packages(include=['tacs*']),
      ext_modules=cythonize(exts, include_path=inc_dirs))
