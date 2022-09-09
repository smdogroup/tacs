# A demonstration of a simple optimization problem in TACS: minimize the
# mass of the CRM model subject to a global stress aggregate constraint
# enforcing that the maximum stress at any quadrature point is less than
# a specified upper bound.
from __future__ import print_function

# Import necessary libraries
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from paropt import ParOpt


class uCRM_VonMisesMassMin(ParOpt.pyParOptProblem):
    """
    Mass minimization with a von Mises stress constraint
    """

    def __init__(self, comm, bdf_name):
        self.comm = comm
        struct_mesh = TACS.MeshLoader(self.comm)
        struct_mesh.scanBDFFile(bdf_name)

        # Set constitutive properties
        rho = 2500.0  # density, kg/m^3
        E = 70e9  # elastic modulus, Pa
        nu = 0.3  # poisson's ratio
        kcorr = 5.0 / 6.0  # shear correction factor
        ys = 350e6  # yield stress, Pa
        min_thickness = 0.002
        max_thickness = 0.20
        thickness = 0.02

        # Loop over components, creating stiffness and element object for each
        num_components = struct_mesh.getNumComponents()
        for i in range(num_components):
            descriptor = struct_mesh.getElementDescript(i)

            # Set the design variable index
            design_variable_index = i
            stiff = constitutive.isoFSDT(
                rho,
                E,
                nu,
                kcorr,
                ys,
                thickness,
                design_variable_index,
                min_thickness,
                max_thickness,
            )
            element = None

            # Create the element object
            if descriptor in ["CQUAD", "CQUADR", "CQUAD4"]:
                element = elements.MITCShell(2, stiff, component_num=i)
            struct_mesh.setElement(i, element)

        # Create tacs assembler object from mesh loader
        self.assembler = struct_mesh.createTACS(6)

        # Create the KS Function
        ksWeight = 50.0
        self.funcs = [
            functions.StructuralMass(self.assembler),
            functions.KSFailure(self.assembler, ksWeight=ksWeight),
        ]

        # Create the forces
        self.forces = self.assembler.createVec()
        force_array = self.forces.getArray()
        force_array[2::6] += 100.0  # uniform load in z direction
        self.assembler.applyBCs(self.forces)

        # Set up the solver
        self.ans = self.assembler.createVec()
        self.res = self.assembler.createVec()
        self.adjoint = self.assembler.createVec()
        self.dfdu = self.assembler.createVec()
        self.mat = self.assembler.createFEMat()
        self.pc = TACS.Pc(self.mat)
        subspace = 100
        restarts = 2
        self.gmres = TACS.KSM(self.mat, self.pc, subspace, restarts)

        # Scale the mass objective so that it is O(10)
        self.mass_scale = 1e-3

        # Scale the thickness variables so that they are measured in
        # mm rather than meters
        self.thickness_scale = 1000.0

        # The number of thickness variables in the problem
        self.nvars = num_components

        # The number of constraints (1 global stress constraint that
        # will use the KS function)
        self.ncon = 1

        # Initialize the base class - this will run the same problem
        # on all processors
        super(uCRM_VonMisesMassMin, self).__init__(MPI.COMM_SELF, self.nvars, self.ncon)

        # Set the inequality options for this problem in ParOpt:
        # The dense constraints are inequalities c(x) >= 0 and
        # use both the upper/lower bounds
        self.setInequalityOptions(dense_ineq=True, use_lower=True, use_upper=True)

        # For visualization
        flag = (
            TACS.ToFH5.NODES
            | TACS.ToFH5.DISPLACEMENTS
            | TACS.ToFH5.STRAINS
            | TACS.ToFH5.EXTRAS
        )
        self.f5 = TACS.ToFH5(self.assembler, TACS.PY_SHELL, flag)
        self.iter_count = 0

        return

    def getVarsAndBounds(self, x, lb, ub):
        """Set the values of the bounds"""
        xvals = np.zeros(self.nvars, TACS.dtype)
        self.assembler.getDesignVars(xvals)
        x[:] = self.thickness_scale * xvals

        xlb = np.zeros(self.nvars, TACS.dtype)
        xub = np.zeros(self.nvars, TACS.dtype)
        self.assembler.getDesignVarRange(xlb, xub)
        lb[:] = self.thickness_scale * xlb
        ub[:] = self.thickness_scale * xub

        return

    def evalObjCon(self, x):
        """Evaluate the objective and constraint"""
        # Evaluate the objective and constraints
        fail = 0
        con = np.zeros(1)

        # Set the new design variable values
        self.assembler.setDesignVars(x[:] / self.thickness_scale)

        # Assemble the Jacobian and factor the matrix
        alpha = 1.0
        beta = 0.0
        gamma = 0.0
        self.assembler.zeroVariables()
        self.assembler.assembleJacobian(alpha, beta, gamma, self.res, self.mat)
        self.pc.factor()

        # Solve the linear system and set the varaibles into TACS
        self.gmres.solve(self.forces, self.ans)
        self.assembler.setVariables(self.ans)

        # Evaluate the function
        fvals = self.assembler.evalFunctions(self.funcs)

        # Set the mass as the objective
        fobj = self.mass_scale * fvals[0]

        # Set the KS function (the approximate maximum ratio of the
        # von Mises stress to the design stress) so that
        # it is less than or equal to 1.0
        con[0] = 1.0 - fvals[1]  # ~= 1.0 - max (sigma/design) >= 0

        return fail, fobj, con

    def evalObjConGradient(self, x, g, A):
        """Evaluate the objective and constraint gradient"""
        fail = 0

        # Evaluate the derivative of the mass and place it in the
        # objective gradient
        gx = np.zeros(self.nvars, TACS.dtype)
        self.assembler.evalDVSens(self.funcs[0], gx)
        g[:] = self.mass_scale * gx / self.thickness_scale

        # Compute the total derivative w.r.t. material design variables
        dfdx = np.zeros(self.nvars, TACS.dtype)
        product = np.zeros(self.nvars, TACS.dtype)

        # Compute the derivative of the function w.r.t. the state
        # variables
        self.assembler.evalDVSens(self.funcs[1], dfdx)
        self.assembler.evalSVSens(self.funcs[1], self.dfdu)
        self.gmres.solve(self.dfdu, self.adjoint)

        # Compute the product of the adjoint with the derivative of the
        # residuals
        self.assembler.evalAdjointResProduct(self.adjoint, product)

        # Set the constraint gradient
        A[0][:] = -(dfdx - product) / self.thickness_scale

        # Write out the solution file every 10 iterations
        if self.iter_count % 10 == 0:
            self.f5.writeToFile("ucrm_iter%d.f5" % (self.iter_count))
        self.iter_count += 1

        return fail


# Load structural mesh from BDF file
tacs_comm = MPI.COMM_WORLD
bdf_name = "CRM_box_2nd.bdf"

crm_opt = uCRM_VonMisesMassMin(tacs_comm, bdf_name)

# Set up the optimization problem
max_lbfgs = 5
opt = ParOpt.pyParOpt(crm_opt, max_lbfgs, ParOpt.BFGS)
opt.setOutputFile("crm_opt.out")

# Set optimization parameters
opt.checkGradients(1e-6)

# Set optimization parameters
opt.setArmijoParam(1e-5)
opt.optimize()

# Get the optimized point
x, z, zw, zl, zu = opt.getOptimizedPoint()
