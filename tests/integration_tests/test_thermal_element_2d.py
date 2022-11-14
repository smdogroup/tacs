import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
The nominal case is a heat conduction problem of a
1m radius plate with a Dirichilet boundary condition applied at the edges,
such that:
    T(theta) = T0 + dT * sin(2*theta)
    T0 = 70 C
    dT = 30 C
The problem is then to solve for the temperature within the boundaries of the plate.
The problem basically boils down to Laplaces problem:
    grad**2 T = 0

test KSTemperature, StructuralMass, CenterOfMass, MomentOfInertia, and AverageTemperature functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/circ-plate-dirichlet-bcs.bdf")

# Radius of plate
R = 1.0
# Area of plate
area = np.pi * R**2

# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "steady_state_avg_temp": 69.8801609399151,
        "steady_state_ks_temp": 98.74014374789103,
        "steady_state_mass": 39.20272476980967,
        "steady_state_x_cg": 2.815920682086164e-10,
        "steady_state_y_cg": 2.826318842831093e-10,
        "steady_state_Ixx": 9.783919839192055,
        "steady_state_Ixy": 2.640029789051368e-08,
        "steady_state_Iyy": 9.783919795061697,
        "transient_avg_temp": 79333.52527756922,
        "transient_ks_temp": 97.89029199587861,
        "transient_mass": 78405.4495396193,
        "transient_x_cg": 2.8159209496185734e-10,
        "transient_y_cg": 2.826321489895307e-10,
        "transient_Ixx": 19567.83967838426,
        "transient_Ixy": 5.280059696105566e-05,
        "transient_Iyy": 19567.83959012328,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Instantiate FEA Assembler
        struct_options = {
            "outputElement": TACS.PLANE_STRESS_ELEMENT,
            # Finer tol needed to pass complex sens test
            "L2Convergence": 1e-16,
        }

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
        ):
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
                if elem_descript in ["CQUAD4", "CQUADR"]:
                    basis = elements.LinearQuadBasis()
                elif elem_descript in ["CTRIA3", "CTRIAR"]:
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

        tacs_probs = []

        # Create static problem, loads are already applied through BCs
        sp = fea_assembler.createStaticProblem(name="steady_state")
        tacs_probs.append(sp)

        # Create transient problem, loads are already applied through BCs
        tp = fea_assembler.createTransientProblem(
            name="transient", tInit=0.0, tFinal=2000.0, numSteps=25
        )
        tacs_probs.append(tp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_temp", functions.KSTemperature, ksWeight=100.0)
            problem.addFunction("avg_temp", functions.AverageTemperature, volume=area)
            problem.addFunction(
                "x_cg", functions.CenterOfMass, direction=[1.0, 0.0, 0.0]
            )
            problem.addFunction(
                "y_cg", functions.CenterOfMass, direction=[0.0, 1.0, 0.0]
            )
            problem.addFunction(
                "Ixx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0],
                direction2=[1.0, 0.0],
                aboutCG=True,
            )
            problem.addFunction(
                "Ixy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0],
                direction2=[1.0, 0.0],
                aboutCG=True,
            )
            problem.addFunction(
                "Iyy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0],
                direction2=[0.0, 1.0],
                aboutCG=True,
            )

        return tacs_probs, fea_assembler
