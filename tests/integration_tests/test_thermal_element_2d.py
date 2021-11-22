import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

'''
The nominal case is a heat conduction problem of a
1m radius plate with a Dirichilet boundary condition applied at the edges,
such that:
    T(theta) = T0 + dT * sin(2*theta)
    T0 = 70 C
    dT = 30 C
The problem is then to solve for the temperature within the boundaries of the plate.
The problem basically boils down to Laplaces problem:
    grad**2 T = 0

test KSTemperature, StructuralMass, and AverageTemperature functions and sensitivities
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/circ-plate-dirichlet-bcs.bdf")

FUNC_REFS = {'steady_state_avg_temp': 69.88016093991516, 'steady_state_ks_temp': 98.74014374789108,
             'steady_state_mass': 39.20272476980967,

             'transient_avg_temp': 396.66762638787424, 'transient_ks_temp': 97.83730882226564,
             'transient_mass': 392.027247698097}


# Radius of plate
R = 1.0
# Area of plate
area = np.pi * R ** 2

# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_pytacs(self, comm, dtype):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Instantiate FEA Assembler
        struct_options = {'outputElement': TACS.PLANE_STRESS_ELEMENT,
                          # Finer tol needed to pass complex sens test
                          'L2Convergence': 1e-16}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elem_call_back(dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs):
            # Material properties
            rho = 2500.0  # density kg/m^3
            kappa = 230.0  # Thermal conductivity W/(m⋅K)

            # Plate geometry
            tplate = 0.005  # 5 mm

            # Setup property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, kappa=kappa)
            # Set one thickness dv for every component
            con = constitutive.PlaneStressConstitutive(prop, t=tplate, tNum=dv_num)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elem_list = []
            model = elements.HeatConduction2D(con)
            for elem_descript in elem_descripts:
                if elem_descript in ['CQUAD4', 'CQUADR']:
                    basis = elements.LinearQuadBasis()
                elif elem_descript in ['CTRIA3', 'CTRIAR']:
                    basis = elements.LinearTriangleBasis()
                else:
                    print("Uh oh, '%s' not recognized" % (elem_descript))
                elem = elements.Element2D(model, basis)
                elem_list.append(elem)

            # Add scale for thickness dv
            scale = [100.0]
            return elem_list, scale

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

        # Define perturbation array that moves all nodes on plate
        xpts = fea_assembler.getOrigNodes()
        xpts_pert_vec[:] = xpts

        return

    def setup_funcs(self, fea_assembler, problems):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        # Add Functions
        for problem in problems:
            problem.addFunction('mass', functions.StructuralMass)
            problem.addFunction('ks_temp', functions.KSTemperature,
                                   ksWeight=100.0)
            problem.addFunction('avg_temp', functions.AverageTemperature, volume=area)
        func_list = ['mass', 'ks_temp', 'avg_temp']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        tacs_probs = []

        # Create static problem, loads are already applied through BCs
        sp = fea_assembler.createStaticProblem(name='steady_state')
        tacs_probs.append(sp)

        # Create transient problem, loads are already applied through BCs
        tp = fea_assembler.createTransientProblem(name='transient', tInit=0.0, tFinal=10.0, numSteps=25)
        tacs_probs.append(tp)

        return tacs_probs
