TacsBuilder class
=================
In large multiphysics problems, creation and connection of the OpenMDAO can be complicated and time-consuming. The design of MPhys is based on builder classes in order to reduce the burden on the user. Most of the assembly of the OpenMDAO model with MPhys is handled by a set of builder helper objects.

Developers wishing to integrate their code to MPhys should subclass the builder and implement the methods relevant to their code. Not all builders need to implement all the methods. For example, a transfer scheme builder may not need a precoupling post coupling subsystem in the scenario.

Options
-------

assembler_setup function
^^^^^^^^^^^^^^^^^^^^^^^^
:class:`~tacs.pytacs.pyTACS` can be initialized py passing a user-defined ``elemCallBack`` function to
:meth:`pyTACS.initialize <tacs.pytacs.pyTACS.initialize>`, that will be used to setup the correct
TACS elements at runtime.

The ``assembler_setup`` function should have the following structure:

  .. function:: assembler_setup(fea_assembler)

    User-defined function used by :class:`~tacs.pytacs.pyTACS` to setup a :ref:`Element<core/elements:Element classes>`
    for each element type in a given component (property) group.

    :param fea_assembler: pyTACS assembler created by Builder.
    :type fea_assembler: tacs.pytacs.pyTACS

  .. code-block:: python

    def assembler_setup(fea_assembler):
        fea_assembler.assignMassDV("engine_mass", 10401)
        fea_assembler.assignMassDV("fuel_mass", 10404)


elem_callback function
^^^^^^^^^^^^^^^^^^^^^^
See :func:`~tacs.pytacs.elemCallBack`.

  .. code-block:: python

    def elem_callback(dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs):
        # Material properties
        rho = 2500.0        # density kg/m^3
        E = 70e9            # Young's modulus (Pa)
        nu = 0.3            # Poisson's ratio
        ys = 464.0e6        # yield stress

        # Plate geometry
        tplate = 0.005    # 5 mm

        # Set up material properties
        prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
        # Set up constitutive model
        con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dv_num)
        # Set the transform used to define shell stresses, None defaults to NaturalShellTransform
        transform = None
        # Set up tacs element for every entry in elem_descripts
        # According to the bdf file, elem_descripts should always be ["CQUAD4"]
        elem_list = []
        for descript in elem_descripts:
            if descript == 'CQUAD4':
                elem = elements.Quad4Shell(transform, con)
            else: # Add a catch for any unexpected element types
                raise ValueError(f"Unexpected element of type {descript}.")
        return elem_list


problem_setup function
^^^^^^^^^^^^^^^^^^^^^^
This function is called each time a new MPhys Scenario is created.
This function sets up problem adding fixed loads, modifying options, and adding eval functions.

The function should have the following structure:

  .. function:: problem_setup(scenario_name, fea_assembler, problem)

    User-defined function used by :class:`~tacs.pytacs.pyTACS` to setup a :ref:`Element<core/elements:Element classes>`
    for each element type in a given component (property) group.

    :param scenario_name: The name of the mphys Scenario calling the function.
    :type scenario_name: str
    :param fea_assembler: pyTACS assembler created by Builder.
    :type fea_assembler: tacs.pytacs.pyTACS
    :param problem: The component description label read in from optional
             formatted comments in BDF file.
    :type problem: tacs.problems.StaticProblem

  .. code-block:: python

     def problem_setup(scenario_name, fea_assembler, problem):
        """
        Helper function to add fixed forces and eval functions
        to structural problems used in tacs builder
        """
        # Set scenario to its own output directory
        problem.setOption('outputDir', scenario_name)

        # Only include mass from elements that belong to pytacs components (i.e. skip concentrated masses)
        comp_ids = fea_assembler.selectcomp_ids(nGroup=-1)
        problem.addFunction('struct_mass', functions.StructuralMass, comp_ids=comp_ids)
        problem.addFunction('ks_vmfailure', functions.KSFailure,
                            safetyFactor=1.5, ksWeight=100.0)

        # load factor
        if scenario_name == "maneuver_2_5g":
          n = 2.5
        elif scenario_name == "maneuver_m1g":
          n = -1.0
        else:
          n = 1.0
        # Add gravity load
        g = n * np.array([0.0, 0.0, -9.81])  # m/s^2
        problem.addInertialLoad(g)

API Reference
-------------

.. autoclass:: tacs.mphys.builder.TacsBuilder
   :members:
   :undoc-members:
