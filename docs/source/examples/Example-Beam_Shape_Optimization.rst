2D Beam Shape Optimization with MACH
************************************
.. note:: The script for this example can be found under the ``examples/beam/shape_opt.py`` file.

This example demonstrates TACS structural shape optimization using the
:ref:`mach/mach:MACH` interface.
It considers the same cantilevered beam with a tip shear load as the
:ref:`examples/Example-Beam_Optimization:Beam optimization with MPhys` example, but differs
in two key ways:

#. The beam is modeled in 2D using shell elements rather than 1D beam elements,
#. The geometry is optimized by physically warping the finite-element mesh via a
   free-form deformation (FFD) volume rather than by adjusting 1D cross-sectional
   properties.

The beam is discretized using 1001 shell elements along its span and depth.

The optimization problem is:

  **Minimize** the mass of the beam with respect to the depth of the cross-section along the span,
  subject to a maximum stress constraint dictated by the material's yield stress.

In order to change the shape of the FEM, we use a FFD volume
parameterization scheme provided by the `pyGeo <https://github.com/mdolab/pygeo>`_ library.
TACS :ref:`mach/mach:MACH` module is used to link the static problem with pyGeo's FFD volume parameterization.

The figure below shows the initial (unoptimized) beam with the von Mises failure contour and
the FFD box with its control points encapsulating the structure:

.. image:: images/Beam_FFD.png
   :width: 700
   :alt: Initial beam failure contour with FFD volume box and control points

Each spanwise control point section will be parameterized in such a way that the depth of the beam can be optimized at each station.

As was shown in :ref:`examples/Example-Beam_Optimization:Beam optimization with MPhys` example an analytic solution for the optimal spanwise depth profile can be derived and is given below:

.. math::
    d(x) = \sqrt{\frac{6 V (L - x)}{t \, \sigma_y}}

The optimization will be driven by SLSQP via pyoptsparse, with gradients supplied by
TACS' adjoint solver through the :class:`~tacs.mach.struct_problem.StructProblem` wrapper.

First, import required libraries:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:imports-start]
   :end-before: # [docs:imports-end]

Next, define the problem parameters and file paths:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:parameters-start]
   :end-before: # [docs:parameters-end]

Now, define the element callback function used to setup TACS element objects and design variables:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:element-callback-start]
   :end-before: # [docs:element-callback-end]

Create and initialize the pyTACS assembler:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:pytacs-init-start]
   :end-before: # [docs:pytacs-init-end]

Set up the FFD and geometric design variables using `pyGeo <https://github.com/mdolab/pygeo>`_'s DVGeometry:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:dvgeo-setup-start]
   :end-before: # [docs:dvgeo-setup-end]

Create the static problem and add functions of interest:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:static-problem-start]
   :end-before: # [docs:static-problem-end]

Wrap the static problem with the :class:`~tacs.mach.struct_problem.StructProblem` using the MACH interface.
Passing ``DVGeo`` here registers the structural node coordinates with the FFD volume;
nodes are updated automatically before each solve when design variables change:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:struct-problem-start]
   :end-before: # [docs:struct-problem-end]

Define the objective and constraint evaluation function:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:struct-obj-start]
   :end-before: # [docs:struct-obj-end]

Define the sensitivity evaluation function.
:meth:`~tacs.mach.struct_problem.StructProblem.evalFunctionsSens` folds in the DVGeo
chain-rule term automatically, producing sensitivities keyed by the geometric DV name
(``"depth"``).  The structural DV sensitivity (keyed ``"struct"``) is
popped out because it is not used in the pyoptsparse optimization problem:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:struct-sens-start]
   :end-before: # [docs:struct-sens-end]

Set up the optimization problem using pyoptsparse.
:meth:`~tacs.mach.struct_problem.StructProblem.addVariablesPyOpt` registers the TACS
structural design variables (none in this case, since ``tNum=-1``), and
``DVGeo.addVariablesPyOpt`` registers the FFD ``"depth"`` DVs.
The stress constraint is added as a nonlinear inequality using the KS failure aggregation:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:opt-setup-start]
   :end-before: # [docs:opt-setup-end]

Finally, run the optimization:

.. literalinclude:: ../../../examples/beam/shape_opt.py
   :language: python
   :start-after: # [docs:run-opt-start]
   :end-before: # [docs:run-opt-end]

Results
-------

The optimization minimizes the beam mass while keeping the KS failure index below 1.0.
The converged depth profile decreases from root to tip, matching the analytical solution
:math:`d(x) = \sqrt{6V(L-x)/(t\,\sigma_y)}`.

.. image:: images/Beam_Shape_Opt.png
   :width: 700
   :alt: Optimized beam failure contours with the idealized analytic depth profile overlaid
