import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
Test heat conduction of a problem with non-uniform boundary conditions with initial conditions everywhere else.
This is a unit-square domain with dT=100 on the left edge, dT=200 on the right edge, and initial condition of dT=150 in between.

This example is a replication of the example here:
https://kitchingroup.cheme.cmu.edu/blog/2013/03/07/Transient-heat-conduction-partial-differential-equations/
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/unit_plate.bdf")

# Area of plate
area = 1.0

# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "steady_state_avg_temp": 150.00000000000057,
        "steady_state_ks_temp": 198.51811570667203,
        "steady_state_mass": 0.9999999999999951,
        "transient_avg_temp": 750.0000000001098,
        "transient_ks_temp": 198.36960589918542,
        "transient_mass": 4.999999999999405,
    }

    def setup_tacs_problems(self, comm):

        # Instantiate FEA Assembler
        struct_options = {  # Finer tol needed to pass complex sens test
            "L2Convergence": 1e-16
        }

        fea_assembler = pytacs.pyTACS(bdf_file, comm)

        # Specify the plate thickness
        tplate = 1.0

        # Material properties
        rho = 1.0  # density kg/m^3
        kappa = 0.02  # Thermal conductivity W/(m⋅K)
        cp = 1.0  # Specific heat J/(kg⋅K)

        # The callback function to define the element properties
        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
        ):

            # Setup property and constitutive objects
            prop = constitutive.MaterialProperties(
                rho=rho, kappa=kappa, specific_heat=cp
            )

            # Set one thickness value for every component
            con = constitutive.PlaneStressConstitutive(prop, t=tplate, tNum=-1)

            model = elements.HeatConduction2D(con)
            basis = elements.LinearQuadBasis()
            elem = elements.Element2D(model, basis)

            return elem

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        tacs_probs = []

        # Create static problem, loads are already applied through BCs
        sp = fea_assembler.createStaticProblem(name="steady_state")
        tacs_probs.append(sp)

        # Create transient problem, loads are already applied through BCs
        tp = fea_assembler.createTransientProblem(
            name="transient", tInit=0.0, tFinal=5.0, numSteps=100
        )
        # Set the initial conditions
        tp.setInitConditions(vars=150.0)
        tacs_probs.append(tp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_temp", functions.KSTemperature, ksWeight=ksweight)
            problem.addFunction("avg_temp", functions.AverageTemperature, volume=area)

        return tacs_probs, fea_assembler
