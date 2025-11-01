2D Beam Shape Optimization with MACH
************************************
.. note:: The script for this example can be found under the `examples/beam/` directory.

This example demonstrates TACS structural optimization capabilities using the MACH interface.
The beam model is a rectangular beam, cantilevered, with a shear load applied at the tip.
The beam is discretized using shell elements along its span and depth.

The optimization problem is as follows:
Minimize the mass of the beam with respect to the depth of the cross-section along the span,
subject to a maximum stress constraint dictated by the material's yield stress.

In order to change the shape of the FEM, we use a free-form deformation (FFD) volume
parameterization scheme provided by the pyGeo library.

An approximate analytical solution can be derived from beam theory,
by realizing that the stress at any spanwise cross-section in the beam
can be found independently using:
    sigma(x,y) = y*M(x)/I
An analytical solution for this problem can be shown to be:
    d(x) = sqrt(6*V*(L-x)/(t*sigma_y))

The optimization is setup using TACS' MACH module, which provides a wrapper
for pyOptsparse optimization.

.. image:: images/beam_shape_initial.png
  :width: 600
  :alt: Initial beam geometry and mesh

First, import required libraries:

.. code-block:: python

  import numpy as np
  import os

  from pygeo import DVGeometry
  from pyoptsparse import Optimization, SNOPT

  from tacs.mach import StructProblem
  from tacs import pyTACS
  from tacs import elements, constitutive, functions

Next, define the problem parameters and file paths:

.. code-block:: python

  bdf_file = os.path.join(os.path.dirname(__file__), 'Slender_Beam.bdf')
  ffd_file = os.path.join(os.path.dirname(__file__), 'ffd_8_linear.fmt')

  # Beam thickness
  t = 0.01            # m
  # Length of beam
  L = 1.0

  # Material properties
  rho = 2780.0 # kg /m^3
  E = 70.0e9
  nu = 0.0
  ys = 420.0e6

  # Shear force applied at tip
  V = 2.5E4

Now, define the element callback function used to setup TACS element objects and design variables:

.. code-block:: python

  def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
      # Setup (isotropic) property and constitutive objects
      prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
      con = constitutive.IsoShellConstitutive(prop, t=t, tNum=-1)
      # TACS shells are sometimes a little overly-rigid in shear
      # We can reduce this effect by decreasing the drilling regularization
      con.setDrillingRegularization(0.1)
      refAxis = np.array([1.0, 0.0, 0.0])
      transform = elements.ShellRefAxisTransform(refAxis)
      elem = elements.Quad4Shell(transform, con)
      return elem

Create and initialize the pyTACS assembler:

.. code-block:: python

  FEAAssembler = pyTACS(bdf_file)
  FEAAssembler.initialize(element_callback)

Set up the geometric design variables using pyGeo's DVGeometry:

.. code-block:: python

  DVGeo = DVGeometry(fileName=ffd_file)
  # Create reference axis
  nRefAxPts = DVGeo.addRefAxis(name="centerline", alignIndex='i', yFraction=0.5)

  # Set up global design variables
  def depth(val, geo):
      for i in range(nRefAxPts):
          geo.scale_y["centerline"].coef[i] = val[i]

  DVGeo.addGlobalDV(dvName="depth", value=np.ones(nRefAxPts), func=depth, 
                    lower=1e-3, upper=10.0, scale=20.0)

.. image:: images/beam_ffd_setup.png
  :width: 600
  :alt: FFD volume and control points for geometric parameterization

Create the static problem and add functions of interest:

.. code-block:: python

  staticProb = FEAAssembler.createStaticProblem("tip_shear")
  # Add TACS Functions
  staticProb.addFunction('mass', functions.StructuralMass)
  staticProb.addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.0, ksWeight=100.0)
  # Add forces to static problem
  staticProb.addLoadToNodes(1112, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)

Create the StructProblem using the MACH interface:

.. code-block:: python

  structProb = StructProblem(staticProb, FEAAssembler, DVGeo=DVGeo)

Define the objective and constraint evaluation function:

.. code-block:: python

  def structObj(x):
      """Evaluate the objective and constraints"""
      funcs = {}
      structProb.setDesignVars(x)
      DVGeo.setDesignVars(x)
      structProb.solve()
      structProb.evalFunctions(funcs)
      structProb.writeSolution()
      if structProb.comm.rank == 0:
          print(x)
          print(funcs)

      return funcs, False

Define the sensitivity evaluation function:

.. code-block:: python

  def structSens(x, funcs):
      """Evaluate the objective and constraint sensitivities"""
      funcsSens = {}
      structProb.evalFunctionsSens(funcsSens)
      for func in funcsSens:
          funcsSens[func].pop("struct")
      return funcsSens, False

Set up the optimization problem using pyOptsparse:

.. code-block:: python

  # Now we create the structural optimization problem:
  optProb = Optimization("Mass min", structObj)
  optProb.addObj("tip_shear_mass")
  structProb.addVariablesPyOpt(optProb)
  DVGeo.addVariablesPyOpt(optProb)
  optProb.addCon("tip_shear_ks_vmfailure", upper=1.0)

  optProb.printSparsity()

  optOptions = {
      "Major feasibility tolerance": 1e-4,
      "Major optimality tolerance": 1e-4,
      "Major iterations limit": 200,
      "Minor iterations limit": 150000,
      "Iterations limit": 1000000,
      "Major step limit": 0.1,
      "Function precision": 1.0e-8,
      "Problem Type": "Minimize",
      "New superbasics limit": 500,
      "Penalty parameter": 1e3,
  }

  opt = SNOPT(options=optOptions)

Finally, run the optimization:

.. code-block:: python

  # Finally run the actual optimization
  sol = opt(optProb, sens=structSens, storeSens=False)

.. image:: images/beam_shape_optimized.png
  :width: 600
  :alt: Optimized beam geometry showing depth variation along span

Results
-------

The optimization successfully minimizes the beam mass while satisfying the stress constraint.
The optimal depth profile follows the analytical solution, with the beam depth decreasing
from the root to the tip to maintain constant stress along the span.

Key features demonstrated in this example:

- **Geometric Design Variables**: Using pyGeo's DVGeometry for shape parameterization
- **MACH Interface**: Integration with pyOptsparse for optimization
- **Constraint Handling**: Stress constraints using KS failure criteria
- **Sensitivity Analysis**: Adjoint-based gradient computation for optimization
- **FFD Parameterization**: Free-form deformation for smooth shape changes

The MACH interface provides a clean separation between the structural analysis (TACS)
and the optimization framework (pyOptsparse), making it easy to set up complex
multidisciplinary optimization problems.
