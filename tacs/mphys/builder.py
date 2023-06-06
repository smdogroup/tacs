import numpy as np
import copy

from openmdao.utils.mpi import MPI
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
        if MPI is not None and self.comm.size > 1:
            # Get DVs locally owned by this processor
            local_dvs = self.fea_assembler.getOrigDesignVars()
            local_dvs = local_dvs.astype(float)
            # Size of design variable on this processor
            local_ndvs = self.fea_assembler.getNumDesignVars()
            # Size of design variable vector on each processor
            dv_sizes = self.comm.allgather(local_ndvs)
            # Offsets for global design variable vector
            offsets = np.zeros(self.comm.size, dtype=int)
            offsets[1:] = np.cumsum(dv_sizes)[:-1]
            # Gather the portions of the design variable array distributed across each processor
            tot_ndvs = sum(dv_sizes)
            global_dvs = np.zeros(tot_ndvs, dtype=local_dvs.dtype)
            self.comm.Allgatherv(local_dvs, [global_dvs, dv_sizes, offsets, MPI.DOUBLE])
            # return the global dv array
            return global_dvs
        else:
            return self.fea_assembler.getOrigDesignVars()

    def get_dv_bounds(self):
        """Get arrays containing the lower and upper bounds for the design variables,
        in the form needed by OpenMDAO's `add_design_variable` method.

        Returns
        -------
        list of ndarray
            lower and upper bounds for the design variables
        """
        if MPI is not None and self.comm.size > 1:
            # Get bounds owned by this processor
            local_dv_bounds = self.fea_assembler.getDesignVarRange()
            local_dv_bounds = list(local_dv_bounds)
            local_dv_bounds[0] = local_dv_bounds[0].astype(float)
            local_dv_bounds[1] = local_dv_bounds[1].astype(float)

            # Size of design variable on this processor
            local_ndvs = self.fea_assembler.getNumDesignVars()
            # Size of design variable vector on each processor
            dv_sizes = self.comm.allgather(local_ndvs)
            # Offsets for global design variable vector
            offsets = np.zeros(self.comm.size, dtype=int)
            offsets[1:] = np.cumsum(dv_sizes)[:-1]
            # Gather the portions of the design variable array distributed across each processor
            tot_ndvs = sum(dv_sizes)
            global_dv_bounds = []
            for ii in [0, 1]:
                global_dv_bounds.append(
                    np.zeros(tot_ndvs, dtype=local_dv_bounds[ii].dtype)
                )
                self.comm.Allgatherv(
                    local_dv_bounds[ii],
                    [global_dv_bounds[ii], dv_sizes, offsets, MPI.DOUBLE],
                )
            # return the global dv array
            return global_dv_bounds
        else:
            return self.fea_assembler.getDesignVarRange()

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
