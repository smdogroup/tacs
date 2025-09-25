MACH Interface
==============

The MACH interface provides a specialized interface for coupling TACS with the MDOLab's MACH library codes, particularly for aerostructural analysis and other multidisciplinary optimization workflows. 
This interface is built around the :class:`StructProblem` class, which serves as a bridge between TACS structural analysis capabilities and external solvers.

Overview
--------

The MACH interface is designed for scenarios where TACS needs to be integrated with external computational libraries, such as:

- The aerodynamic solver ADflow for aerostructural coupling
- The geometric parametrization library pygeo for geometric design variables
- The optimization library pyoptsparse for optimization

Key Features
------------

- **External Force Integration**: Seamless handling of forces from external solvers
- **Adjoint Capabilities**: Full support for adjoint-based sensitivity analysis
- **Design Variable Management**: Comprehensive design variable handling with geometric and material parameters
- **Constraint Support**: Built-in support for structural constraints

StructProblem Class
-------------------

The :class:`StructProblem` class is the core component of the MACH interface. It extends the base structural problem functionality to provide enhanced capabilities for external library integration.

.. class:: StructProblem(staticProblem, FEAAssembler, DVGeo=None, loadFile=None)

   Base class for TACS problem types interfacing with external libraries.

   **Parameters:**
   
   * **staticProblem** (*tacs.problems.StaticProblem*): StaticProblem object used for modeling and solving static cases
   * **FEAAssembler** (*tacs.pytacs.pyTACS*): pyTACS assembler object
   * **DVGeo** (*pygeo.DVGeometry, optional*): Object responsible for manipulating geometric design variables
   * **loadFile** (*str, optional*): Filename of external load file for coupling

   **Key Methods:**

   .. method:: solve(damp=1.0, useAitkenAcceleration=False, dampLB=0.2, loadScale=1.0)

      Solve the static problem with optional damping and acceleration for coupling.

      **Parameters:**
      
      * **damp** (*float*): Damping factor for solution update
      * **useAitkenAcceleration** (*bool*): Enable Aitken acceleration for aerostructural problems
      * **dampLB** (*float*): Lower bound for damping parameter
      * **loadScale** (*float*): Scaling factor for external loads

   .. method:: setExternalForce(Fext)

      Set external coupling forces from external solvers.

      **Parameters:**
      
      * **Fext** (*numpy.ndarray or tacs.TACS.Vec*): External force vector

   .. method:: getExternalForce(Fext=None)

      Retrieve current external coupling forces.

      **Returns:**
      
      * **Fext** (*numpy.ndarray*): Current external force vector

   .. method:: readExternalForceFile(fileName)

      Read external loads from file for coupling workflows.

      **Parameters:**
      
      * **fileName** (*str*): Filename of load file (.dat or .bdf extension)

   .. method:: writeExternalForceFile(outputDir=None, baseName=None, number=None)

      Write external loads to file for coupling workflows.

   .. method:: addConstraint(constr)

      Add TACS constraint object to the problem.

      **Parameters:**
      
      * **constr** (*tacs.constraints.TACSConstraint*): Constraint object to add

   .. method:: setDVGeo(DVGeo, pointSetKwargs=None)

      Set the DVGeometry object for geometric design variable manipulation.

      **Parameters:**
      
      * **DVGeo** (*pygeo.DVGeometry*): Geometry manipulation object
      * **pointSetKwargs** (*dict, optional*): Keyword arguments for point set creation

Usage Examples
--------------

Basic Structural Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from tacs.mach import StructProblem
   from tacs import pyTACS
   from tacs import elements, constitutive, functions
   
   # Create pyTACS assembler
   FEAAssembler = pyTACS('model.bdf')
   
   # Element callback function
   def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
       prop = constitutive.MaterialProperties(rho=2700.0, E=70e9, nu=0.3)
       con = constitutive.IsoShellConstitutive(prop, t=0.01, tNum=dvNum)
       transform = elements.ShellRefAxisTransform([1.0, 0.0, 0.0])
       elem = elements.Quad4Shell(transform, con)
       return elem
   
   # Initialize assembler
   FEAAssembler.initialize(elemCallBack)
   
   # Create static problem
   staticProblem = FEAAssembler.createStaticProblem('analysis')
   staticProblem.addFunction('mass', functions.StructuralMass)
   
   # Create StructProblem for MACH interface
   structProblem = StructProblem(staticProblem, FEAAssembler)
   
   # Solve problem
   structProblem.solve()
   
   # Evaluate functions
   funcs = {}
   structProblem.evalFunctions(funcs)
   print(f"Structural mass: {funcs['analysis_mass']}")

Geometric Design Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pygeo import DVGeometry
   
   # Create DVGeometry object
   DVGeo = DVGeometry('ffd.fmt')
   
   # Create StructProblem with geometric design variables
   structProblem = StructProblem(staticProblem, FEAAssembler, DVGeo=DVGeo)
   
   # Set design variables
   designVars = {'ffd': np.array([0.1, -0.05, 0.02])}
   structProblem.setDesignVars(designVars)
   
   # Solve and evaluate sensitivities
   structProblem.solve()
   funcs = {}
   funcsSens = {}
   structProblem.evalFunctions(funcs)
   structProblem.evalFunctionsSens(funcsSens)

Constraint Handling
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from tacs.constraints import AdjacencyConstraint, DVConstraint
   
   # Create StructProblem
   structProblem = StructProblem(staticProblem, FEAAssembler)
   
   # Add adjacency constraint
   adjCon = FEAAssembler.createAdjacencyConstraint("AdjCon")
   adjCon.addConstraint(
       conName="thicknessAdj",
       compIDs=compIDs,
       lower=-0.001,
       upper=0.001,
       dvIndex=0
   )
   structProblem.addConstraint(adjCon)
   
   # Add design variable constraint
   dvCon = FEAAssembler.createDVConstraint("DVCon")
   dvCon.addConstraint(
       conName="thicknessRatio",
       upper=0.0,
       dvIndices=[0, 1],
       dvWeights=[1.0, -2.0]
   )
   structProblem.addConstraint(dvCon)
   
   # Evaluate constraints
   funcs = {}
   structProblem.evalConstraints(funcs)

Related Documentation
---------------------

- :doc:`../pytacs/pytacs` - pyTACS interface documentation
- :doc:`../core/TACS` - Core TACS C++ interface documentation
