MPhys
*****
`MPhys <https://pypi.org/project/mphys/>`_ is a package that standardizes high-fidelity multiphysics problems in OpenMDAO.
MPhys provides a convenient interface for connecting TACS problems to optimizers and other physics disciplinary solvers for coupled analysis.
MPhys is an optional TACS dependency, so in order to use the module within TACS it will need to be installed manually along with its dependencies.
Assuming the user has already installed TACS, the only additional dependency the user will need to install would be `petsc4py <https://pypi.org/project/petsc4py/>`_.
These packages can be installed easily in a conda environment using the commands below:

::

    pip install mphys
    conda install -c conda-forge petsc4py


The TACS MPhys interface consists of one main classes: a builder class called :class:`~tacs.mphys.builder.TacsBuilder`.
For more information on the general MPhys interface users should see the MPhys `docs <https://openmdao.github.io/mphys/>`_.
The details of the :class:`~tacs.mphys.builder.TacsBuilder` interface will be discussed in the sections below.

.. toctree::
  :maxdepth: 1

  builder