.. TACS documentation master file, created by
   sphinx-quickstart on Wed Oct 10 15:04:01 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TACS Overview
=============

The Toolkit for the Analysis of Composite Structures (TACS) is a parallel finite-element code for analysis and gradient-based design of advanced structures. Originally, TACS was primarily designed as a tool for the analysis of shell structures, such as wing-boxes. More recently it has been extended to perform topology optimization of large three-dimensional structures using gradient-based methods.

TACS has been under continuous development since 2010 by the `Structural and Multidisciplinary Design Optimization group at Georgia Tech <http://gkennedy.gatech.edu>`_ and by the `Multidisciplinary Design Optimization Lab at the University of Michigan <http://mdolab.engin.umich.edu/>`_.

Getting Started
===============
.. toctree::
   :maxdepth: 2

   install
   interfaces

Examples
========
.. toctree::
   :maxdepth: 1

   examples/Example-Plate
   examples/Example-Transient_Battery
   examples/Example-Beam_Optimization
   examples/Example-Composite_Optimization

References
==========
.. toctree::
   :maxdepth: 2

   theory/theory
   core/core

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
