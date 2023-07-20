import numpy as np
import copy

from mphys.builder import Builder

from .. import pyTACS
from .mesh import TacsMeshGroup
from .precoupling import TacsPrecouplingGroup
from .coupling import TacsCouplingGroup
from .postcoupling import TacsPostcouplingGroup


class TacsBuilder(Builder):
    def __init__(
        self,
        options,
        check_partials=False,
        conduction=False,
        coupled=True,
        write_solution=True,
        separate_mass_dvs=False,
    ):
        self.options = copy.deepcopy(options)
        self.check_partials = check_partials
        # Flag to switch to tacs conduction solver (False->structural)
        self.conduction = conduction
        # Flag to turn on f5 file writer
        self.write_solution = write_solution
        # Flag to turn on coupling variables
        self.coupled = coupled
        # Flag to separate point mass dvs from struct dvs in openmdao input array
        self.separate_mass_dvs = separate_mass_dvs

    def initialize(self, comm):
        pytacs_options = copy.deepcopy(self.options)
        bdf_file = pytacs_options.pop("mesh_file")

        # Load optional user-defined callback function for setting up tacs elements
        if "assembler_setup" in pytacs_options:
            assembler_setup = pytacs_options.pop("assembler_setup")
        else:
            assembler_setup = None

        # Load optional user-defined callback function for setting up tacs elements
        if "element_callback" in pytacs_options:
            element_callback = pytacs_options.pop("element_callback")
        else:
            element_callback = None

        # Load optional user-defined callback function for setting up tacs elements
        if "problem_setup" in pytacs_options:
            self.problem_setup = pytacs_options.pop("problem_setup")
        else:
            self.problem_setup = None

        # Load optional user-defined callback function for setting up constraints
        if "constraint_setup" in pytacs_options:
            self.constraint_setup = pytacs_options.pop("constraint_setup")
        else:
            self.constraint_setup = None

        # Load optional user-defined callback function for setting up buckling problem
        if "buckling_setup" in pytacs_options:
            self.buckling_setup = pytacs_options.pop("buckling_setup")
        else:
            self.buckling_setup = None

        # Create pytacs instance
        self.fea_assembler = pyTACS(bdf_file, options=pytacs_options, comm=comm)
        self.comm = comm

        # Do any pre-initialize setup requested by user
        if assembler_setup is not None:
            assembler_setup(self.fea_assembler)

        # Set up elements and TACS assembler
        self.fea_assembler.initialize(element_callback)

    def get_coupling_group_subsystem(self, scenario_name=None):
        return TacsCouplingGroup(
            fea_assembler=self.fea_assembler,
            conduction=self.conduction,
            check_partials=self.check_partials,
            coupled=self.coupled,
            scenario_name=scenario_name,
            problem_setup=self.problem_setup,
        )

    def get_mesh_coordinate_subsystem(self, scenario_name=None):
        return TacsMeshGroup(fea_assembler=self.fea_assembler)

    def get_pre_coupling_subsystem(self, scenario_name=None):
        initial_dvs = self.get_initial_dvs()
        return TacsPrecouplingGroup(
            fea_assembler=self.fea_assembler,
            initial_dv_vals=initial_dvs,
            separate_mass_dvs=self.separate_mass_dvs,
        )

    def get_post_coupling_subsystem(self, scenario_name=None):
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
            masked_local_nodes = self.fea_assembler.getLocalNodeIDsForComps(tagged_comps)

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
