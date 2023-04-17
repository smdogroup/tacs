# A demonstration of a simple optimization problem in TACS: minimize the
# mass of the CRM model subject to a global stress aggregate constraint
# enforcing that the maximum stress at any quadrature point is less than
# a specified upper bound.

# Import necessary libraries
import numpy as np
from mpi4py import MPI
from paropt import ParOpt

from tacs import pyTACS, elements, constitutive, functions


class uCRM_VonMisesMassMin(ParOpt.Problem):
    """
    Mass minimization with a von Mises stress constraint
    """

    def __init__(self, comm, bdf_name):
        self.comm = comm

        # Set constitutive properties
        rho = 2500.0  # density, kg/m^3
        E = 70e9  # elastic modulus, Pa
        nu = 0.3  # poisson's ratio
        ys = 350e6  # yield stress, Pa
        min_thickness = 0.002
        max_thickness = 0.20
        thickness = 0.02

        # Callback function used to setup TACS element objects and DVs
        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set one thickness dv for every component
            con = constitutive.IsoShellConstitutive(
                prop, t=thickness, tNum=dvNum, tlb=min_thickness, tub=max_thickness
            )

            # This is an all-quad mesh, return a quad shell
            transform = None
            elem = elements.Quad4Shell(transform, con)
            return elem

        # Create pyTACS assembler object
        self.fea_assembler = pyTACS(bdf_name, comm=comm)
        # Setup elements
        self.fea_assembler.initialize(elemCallBack)

        # Create a static problem with a uniform z load applied to all nodes
        self.static_problem = self.fea_assembler.createStaticProblem("uniform_z_load")
        F = self.fea_assembler.createVec()
        F[2::6] += 100.0
        self.static_problem.addLoadToRHS(F)
        # Add failure function (con)
        self.static_problem.addFunction(
            "ks_failure", functions.KSFailure, ksWeight=50.0
        )
        # add mass (obj)
        self.static_problem.addFunction("mass", functions.StructuralMass)

        # Scale the mass objective so that it is O(10)
        self.mass_scale = 1e-3

        # Scale the thickness variables so that they are measured in
        # mm rather than meters
        self.thickness_scale = 1000.0

        # The number of thickness variables in the problem
        self.nvars = self.fea_assembler.getNumDesignVars()

        # The number of constraints (1 global stress constraint that
        # will use the KS function)
        self.ncon = 1

        # Initialize the base class - this will run the same problem
        # on all processors
        super(uCRM_VonMisesMassMin, self).__init__(MPI.COMM_SELF, self.nvars, self.ncon)

        # Set the inequality options for this problem in ParOpt:
        # The dense constraints are inequalities c(x) >= 0 and
        # use both the upper/lower bounds
        self.setInequalityOptions(sparse_ineq=False, use_lower=True, use_upper=True)

        self.iter_count = 0

        return

    def getVarsAndBounds(self, x, lb, ub):
        """Set the values of the bounds"""
        xvals = self.static_problem.getDesignVars()
        x[:] = self.thickness_scale * xvals

        xlb, xub = self.static_problem.getDesignVarRange()
        lb[:] = self.thickness_scale * xlb
        ub[:] = self.thickness_scale * xub

        return

    def evalObjCon(self, x):
        """Evaluate the objective and constraint"""
        # Evaluate the objective and constraints
        fail = 0
        con = np.zeros(1)
        func_dict = {}

        # Set the new design variable values
        self.static_problem.setDesignVars(x[:] / self.thickness_scale)

        # Solve
        self.static_problem.solve()

        # Evaluate the function
        self.static_problem.evalFunctions(func_dict)

        if self.comm.rank == 0:
            print(func_dict)

        # Set the mass as the objective
        fobj = self.mass_scale * func_dict["uniform_z_load_mass"]

        # Set the KS function (the approximate maximum ratio of the
        # von Mises stress to the design stress) so that
        # it is less than or equal to 1.0
        con[0] = (
            1.0 - func_dict["uniform_z_load_ks_failure"]
        )  # ~= 1.0 - max (sigma/design) >= 0

        return fail, fobj, con

    def evalObjConGradient(self, x, g, A):
        """Evaluate the objective and constraint gradient"""
        fail = 0
        sens_dict = {}

        # Set the new design variable values
        self.static_problem.setDesignVars(x[:] / self.thickness_scale)

        # Evaluate the function sensitivities
        self.static_problem.evalFunctionsSens(sens_dict)

        # objective gradient
        g[:] = (
            self.mass_scale
            * sens_dict["uniform_z_load_mass"]["struct"]
            / self.thickness_scale
        )

        # Set the constraint gradient
        A[0][:] = (
            -sens_dict["uniform_z_load_ks_failure"]["struct"] / self.thickness_scale
        )

        # Write out the solution file every 10 iterations
        if self.iter_count % 10 == 0:
            self.static_problem.writeSolution()
        self.iter_count += 1

        return fail

    def writeSolutionBDF(self, file_name):
        # Write loads and optimized panel properties to BDF
        self.fea_assembler.writeBDF(file_name, self.static_problem)


# Load structural mesh from BDF file
tacs_comm = MPI.COMM_WORLD
bdf_name = "CRM_box_2nd.bdf"

crm_opt = uCRM_VonMisesMassMin(tacs_comm, bdf_name)

# Set up the optimization problem
options = {
    "algorithm": "tr",
    "tr_init_size": 0.05,
    "tr_min_size": 1e-6,
    "tr_max_size": 10.0,
    "tr_eta": 0.25,
    "tr_infeas_tol": 1e-6,
    "tr_l1_tol": 1e-3,
    "tr_linfty_tol": 0.0,
    "tr_adaptive_gamma_update": True,
    "tr_max_iterations": 1000,
    "max_major_iters": 100,
    "penalty_gamma": 1e3,
    "qn_subspace_size": 10,
    "qn_type": "bfgs",
    "abs_res_tol": 1e-8,
    "starting_point_strategy": "affine_step",
    "barrier_strategy": "mehrotra_predictor_corrector",
    "use_line_search": False,
}
opt = ParOpt.Optimizer(crm_opt, options)

# Run the optimization
opt.optimize()

# Get the optimized point
x, z, zw, zl, zu = opt.getOptimizedPoint()

# Write final design to BDF
crm_opt.writeSolutionBDF("opt_sol.bdf")
