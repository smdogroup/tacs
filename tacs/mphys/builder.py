import copy
import warnings

from mphys.builder import Builder
import numpy as np

from tacs.pytacs import pyTACS
from tacs.mphys.mesh import TacsMeshGroup
from tacs.mphys.precoupling import TacsPrecouplingGroup
from tacs.mphys.coupling import TacsCouplingGroup
from tacs.mphys.postcoupling import TacsPostcouplingGroup


class TacsBuilder(Builder):
    def __init__(
        self,
        mesh_file,
        assembler_setup=None,
        element_callback=None,
        problem_setup=None,
        constraint_setup=None,
        buckling_setup=None,
        pytacs_options=None,
        check_partials=False,
        conduction=False,
        coupled=True,
        write_solution=True,
        separate_mass_dvs=False,
        res_ref=None,
    ):
        """
        Create the Builder responsible for creating MPhys Scenario's based on TACS Analyses.

        Parameters
        ----------
        mesh_file : str or pyNastran.bdf.bdf.BDF
            The BDF file or a pyNastran BDF object to load.
        assembler_setup : collections.abc.Callable, optional
            User-defined callback function for modifying pyTACS assembler prior to initialization.
            This is used for adding point mass DVs to the model. Defaults to None.
            follows::

                def assembler_setup(fea_assembler):

            fea_assembler is the uninitialized pyTACS assembler instance created by the builder.
        element_callback : collections.abc.Callable, optional
            User-defined callback function for setting up TACS elements and element DVs. Defaults to None.
            See :ref:`pytacs/pytacs_module:Initializing with elemCallBack` for more info.
            The calling sequence for elem_callback **must** be as follows::

                def elem_callback(dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs):

            The dv_num is the current counter which must be used by the
            user when creating a constitutive object with design
            variables.

            comp_id is the ID number used by tacs to reference this property group.
            Use kwargs['propID'] to get the corresponding Nastran property ID that
            is read in from the BDF.

            comp_descript is the component description label read in from optional
            formatted comments in BDF file

            elem_descripts are the name of the elements belonging to this group
            (e.g. CQUAD4, CTRIA3, CTETRA, etc). This value will be a list since
            one component may contain multiple compatible element types.
            Example: ['CQUAD4', CTRIA3']

            global_dvs is a dictionary containing information about any
            global DVs that have been added.
            elem_callback must return a list containing as many TACS element
            objects as there are element types in elemDescripts (one for each).
        problem_setup : collections.abc.Callable, optional
            This function is called each time a new MPhys Scenario is created. This function sets up the problem by adding
            fixed loads, modifying options, and adding eval functions. The function should have the following structure::

              def problem_setup(scenario_name, fea_assembler, problem):

            scenario_name is a str denoting which MPhys Scenario the problem is currently being set up for.
            fea_assembler is the uninitialized pyTACS assembler instance created by the builder.
            problem is the tacs.BaseProblem class being setup for this scenario.
        constraint_setup : collections.abc.Callable, optional
            This function is called each time a new MPhys Scenario is created. This function sets up a series of
            constraints to be run after at the end of an MPhys analysis.
            The function should have the following structure::

                def constraint_setup(scenario_name, fea_assembler, constraints):

        buckling_setup : collections.abc.Callable, optional
            This function is called each time a new MPhys Scenario is created. This function sets up a buckling problem
            to be run after at the end of an MPhys analysis. The function should have the following structure::

                def buckling_setup(scenario_name, fea_assembler)

        pytacs_options : dict, optional
            Options dictionary passed to pyTACS assembler.
        check_partials : bool, optional
            This flag allows TACS components partial derivative routines to be evaluated in forward mode without raising
            an error. This lets OpenMDAO's check_partials routine run without errors, allowing users to check TACS'
            reverse derivatives. The forward derivative checks will still not be meaningful since TACS only supports
            reverse mode. Defaults to False.
        conduction : bool, optional
            Flag to determine weather TACS component represents a thermal (True) or structural (False) analysis.
            Defaults to False.
        coupled : bool, optional
            Flag to determine of if multidisciplinary coupling variables should be turned on
            (used in aerostructural/thermostructural analyses). Defaults to True.
        write_solution : bool, optional
            Flag to determine whether to write out TACS solutions to f5 file each design iteration. Defaults to True.
        separate_mass_dvs : bool, optional
            Flag to determine if TACS' mass dvs should be lumped into the struct_dv input vector (False) or
            split into separate OpenMDAO inputs based on their assigned names (True). Defaults to False.
        res_ref : float, optional
            Reference residual norm to be used by OpenMDAO's residual scaling. Can be useful for ensuring residuals
            from different coupled disciplines are of a similar order of magnitude. Defaults to None.

        Examples
        --------
        assembler_setup:
            >>>  def assembler_setup(fea_assembler):
            >>>      # Assign dvs for point mass elements 10401 and 10402
            >>>      # to vary mass values during optimization
            >>>      fea_assembler.assignMassDV("engine_mass", 10401)
            >>>      fea_assembler.assignMassDV("fuel_mass", 10402)
        element_callback:
            >>>  def elem_callback(dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs):
            >>>      # Material properties
            >>>      rho = 2500.0        # density kg/m^3
            >>>      E = 70e9            # Young's modulus (Pa)
            >>>      nu = 0.3            # Poisson's ratio
            >>>      ys = 464.0e6        # yield stress
            >>>
            >>>      # Plate geometry
            >>>      tplate = 0.005    # 5 mm
            >>>
            >>>      # Set up material properties
            >>>      prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            >>>      # Set up constitutive model
            >>>      con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dv_num)
            >>>      # Set the transform used to define shell stresses, None defaults to NaturalShellTransform
            >>>      transform = None
            >>>      # Set up tacs element for every entry in elem_descripts
            >>>      # According to the bdf file, elem_descripts should always be ["CQUAD4"]
            >>>      elem_list = []
            >>>      for descript in elem_descripts:
            >>>          if descript == 'CQUAD4':
            >>>              elem = elements.Quad4Shell(transform, con)
            >>>          else: # Add a catch for any unexpected element types
            >>>              raise ValueError(f"Unexpected element of type {descript}.")
            >>>      return elem_list
        problem_setup:
            >>> def problem_setup(scenario_name, fea_assembler, problem):
            >>>     # Set scenario to its own output directory
            >>>     problem.setOption('outputDir', scenario_name)
            >>>
            >>>     # Only include mass from elements that belong to pytacs components (i.e. skip concentrated masses)
            >>>     comp_ids = fea_assembler.selectCompIDs(nGroup=-1)
            >>>     problem.addFunction('struct_mass', functions.StructuralMass, comp_ids=comp_ids)
            >>>     problem.addFunction('ks_vmfailure', functions.KSFailure,
            >>>                         safetyFactor=1.5, ksWeight=100.0)
            >>>
            >>>     # load factor
            >>>     if scenario_name == "maneuver_2_5g":
            >>>       n = 2.5
            >>>     elif scenario_name == "maneuver_m1g":
            >>>       n = -1.0
            >>>     else:
            >>>       n = 1.0
            >>>     # Add gravity load
            >>>     g = n * np.array([0.0, 0.0, -9.81])  # m/s^2
            >>>     problem.addInertialLoad(g)

        constraint_setup:
            >>> def constraint_setup(scenario_name, fea_assembler, constraints):
            >>>     # Add constraint on enclosed volume of structure
            >>>     constr = fea_assembler.createVolumeConstraint("constraints")
            >>>     constr.addConstraint("volume")
            >>>     constraints.append(constr)

        buckling_setup:
            >>> def buckling_setup(scenario_name, fea_assembler):
            >>>     # Add buckling analysis only to 2.5g maneuver scenario
            >>>     if scenario_name == "maneuver_2_5g":
            >>>         problem = fea_assembler.createBucklingProblem(
            >>>             "buckling", sigma=1.0, numEigs=2
            >>>         )
            >>>         problem.setOption("L2Convergence", 1e-20)
            >>>         problem.setOption("L2ConvergenceRel", 1e-20)
            >>>         return problem

        """
        if isinstance(mesh_file, dict):
            warnings.warn(
                "The signature for TacsBuilder has changed. Arguments such as 'mesh_file' must be passed directly as "
                "arguments to TacsBuilder. Please see the TacsBuilder docstring for more info. This will become an "
                "error in tacs 3.7.0",
                DeprecationWarning,
            )
            options = mesh_file
            # Make deep copy of dict so we can modify it
            pytacs_options = copy.deepcopy(options)
            mesh_file = pytacs_options.pop("mesh_file")

            # Load optional user-defined callback function for setting up tacs elements
            if "assembler_setup" in pytacs_options:
                assembler_setup = pytacs_options.pop("assembler_setup")

            # Load optional user-defined callback function for setting up tacs elements
            if "element_callback" in pytacs_options:
                element_callback = pytacs_options.pop("element_callback")

            # Load optional user-defined callback function for setting up tacs elements
            if "problem_setup" in pytacs_options:
                problem_setup = pytacs_options.pop("problem_setup")

            # Load optional user-defined callback function for setting up constraints
            if "constraint_setup" in pytacs_options:
                constraint_setup = pytacs_options.pop("constraint_setup")

            # Load optional user-defined callback function for setting up buckling problem
            if "buckling_setup" in pytacs_options:
                buckling_setup = pytacs_options.pop("buckling_setup")

        self.mesh_file = mesh_file
        self.assembler_setup = assembler_setup
        self.element_callback = element_callback
        self.problem_setup = problem_setup
        self.constraint_setup = constraint_setup
        self.buckling_setup = buckling_setup
        self.pytacs_options = pytacs_options
        self.check_partials = check_partials
        self.conduction = conduction
        self.coupled = coupled
        self.write_solution = write_solution
        self.separate_mass_dvs = separate_mass_dvs
        self.res_ref = res_ref

    def initialize(self, comm):
        """
        Initialize pyTACS.
        This method will be called when the MPI comm is available

        Parameters
        ----------
        comm : :class:`~mpi4py.MPI.Comm`
            The communicator object created for this xfer object instance.
        """
        # Create pytacs instance
        self.fea_assembler = pyTACS(
            self.mesh_file, options=self.pytacs_options, comm=comm
        )
        self.comm = comm

        # Do any pre-initialize setup requested by user
        if self.assembler_setup is not None:
            self.assembler_setup(self.fea_assembler)

        # Set up elements and TACS assembler
        self.fea_assembler.initialize(self.element_callback)

    def get_coupling_group_subsystem(self, scenario_name=None):
        """
        The subsystem that this builder will add to the CouplingGroup

        Parameters
        ----------
        scenario_name : str, optional
            The name of the scenario calling the builder.

        Returns
        -------
        subsystem : openmdao.api.Group
            The openmdao subsystem that handles all the computations for
            this solver. Transfer schemes can return multiple subsystems
        """
        return TacsCouplingGroup(
            fea_assembler=self.fea_assembler,
            conduction=self.conduction,
            check_partials=self.check_partials,
            coupled=self.coupled,
            scenario_name=scenario_name,
            problem_setup=self.problem_setup,
            res_ref=self.res_ref,
        )

    def get_mesh_coordinate_subsystem(self, scenario_name=None):
        """
        The subsystem that contains the subsystem that will return the mesh
        coordinates

        Parameters
        ----------
        scenario_name : str, optional
            The name of the scenario calling the builder.

        Returns
        -------
        mesh : :class:`~openmdao.api.Component` or :class:`~openmdao.api.Group`
            The openmdao subsystem that has an output of coordinates.
        """
        return TacsMeshGroup(fea_assembler=self.fea_assembler)

    def get_pre_coupling_subsystem(self, scenario_name=None):
        """
        Method that returns the openmdao subsystem to be added to each scenario before the coupling group

        Parameters
        ----------
        scenario_name : str, optional
            The name of the scenario calling the builder.

        Returns
        -------
        subsystem : openmdao.api.Group
        """
        initial_dvs = self.get_initial_dvs()
        return TacsPrecouplingGroup(
            fea_assembler=self.fea_assembler,
            initial_dv_vals=initial_dvs,
            separate_mass_dvs=self.separate_mass_dvs,
        )

    def get_post_coupling_subsystem(self, scenario_name=None):
        """
        Method that returns the openmdao subsystem to be added to each scenario after the coupling group

        Parameters
        ----------
        scenario_name : str, optional
            The name of the scenario calling the builder.

        Returns
        -------
        subsystem : openmdao.api.Group
        """
        return TacsPostcouplingGroup(
            fea_assembler=self.fea_assembler,
            check_partials=self.check_partials,
            conduction=self.conduction,
            write_solution=self.write_solution,
            scenario_name=scenario_name,
            problem_setup=self.problem_setup,
            constraint_setup=self.constraint_setup,
            buckling_setup=self.buckling_setup,
        )

    def get_ndof(self):
        """
        The number of degrees of freedom used at each output location.

        Returns
        -------
        ndof : int
            number of degrees of freedom of each node in the domain
        """
        return self.fea_assembler.getVarsPerNode()

    def get_number_of_nodes(self):
        """
        Get the number of nodes on this processor,
        not including lagrange multiplier nodes
        """
        nnodes = self.fea_assembler.getNumOwnedNodes()
        nmult = self.fea_assembler.getNumOwnedMultiplierNodes()
        return nnodes - nmult

    def get_initial_dvs(self):
        """
        Get an array holding all dvs values that have been added to TACS
        """
        local_dvs = self.fea_assembler.getOrigDesignVars()
        all_local_dvs = self.comm.allgather(local_dvs)
        global_dvs = np.concatenate(all_local_dvs)
        return global_dvs.astype(float)

    def get_dv_bounds(self):
        """Get arrays containing the lower and upper bounds for the design variables,
        in the form needed by OpenMDAO's `add_design_variable` method.

        Returns
        -------
        list of ndarray
            lower and upper bounds for the design variables
        """
        local_lb, local_ub = self.fea_assembler.getDesignVarRange()
        all_lb = self.comm.allgather(local_lb)
        global_lbs = np.concatenate(all_lb)
        all_ub = self.comm.allgather(local_ub)
        global_ubs = np.concatenate(all_ub)
        return global_lbs.astype(float), global_ubs.astype(float)

    def get_dv_scalers(self):
        """Get an array containing the scaling factors for the design
        variables, in the form needed by OpenMDAO's `add_design_variable`
        method.

        Returns
        -------
        array
            Scaling values
        """
        return np.array(self.fea_assembler.scaleList)

    def get_ndv(self):
        """
        Get total number of structural design variables across all procs
        """
        return self.fea_assembler.getTotalNumDesignVars()

    def get_solver(self):
        # this method is only used by the RLT transfer scheme
        return self.fea_assembler.assembler

    def get_fea_assembler(self):
        """
        Returns underlying pytacs object.
        """
        return self.fea_assembler

    def get_tagged_indices(self, tags):
        """
        Method that returns grid IDs for a list of body/boundary tags.

        Parameters
        ----------
        tags : list[str]

        Returns
        -------
        grid_ids : list[int]
            list of grid IDs that correspond to given body/boundary tags
        """
        # Select all compIDs
        if tags == -1 or tags == [-1]:
            nnodes = self.fea_assembler.getNumOwnedNodes()
            # Select all node IDs
            masked_local_nodes = np.arange(nnodes)

        # Get the compIDs associated with tags
        else:
            tagged_comps = self.fea_assembler.selectCompIDs(include=tags)
            # Select local node IDs for tags
            masked_local_nodes = self.fea_assembler.getLocalNodeIDsForComps(
                tagged_comps
            )

        # Select local node IDs and multiplier node IDs
        local_mnodes = self.fea_assembler.getLocalMultiplierNodeIDs()

        # Loop through the multiplier nodes and remove them
        masked_local_nodes = list(masked_local_nodes)
        for mult_node in local_mnodes:
            if mult_node in masked_local_nodes:
                masked_local_nodes.remove(mult_node)
        masked_local_nodes = np.array(masked_local_nodes)

        # Loop through the multiplier nodes and offset for the multiplier nodes we removed
        for mult_node in reversed(local_mnodes):
            masked_local_nodes[masked_local_nodes > mult_node] -= 1

        return list(masked_local_nodes)
