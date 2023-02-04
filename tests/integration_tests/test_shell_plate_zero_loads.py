import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

""""
The nominal case is a 0.5m x 0.5m flat plate under no load. The solution to this is trivial, 
but its still worth testing. The plate comprises 100 CQUAD4 elements and test KSFailure, 
KSDisplacement, StructuralMass, and Compliance functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")


# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "zero_load_static_compliance": 0.0,
        "zero_load_static_ks_disp": -0.13862943611198988,
        "zero_load_static_ks_vmfailure": -0.13169796430639044,
        "zero_load_static_mass": 3.125,
        "zero_load_transient_compliance": 0.0,
        "zero_load_transient_ks_disp": -0.13862943611199358,
        "zero_load_transient_ks_vmfailure": -0.13169796430639413,
        "zero_load_transient_mass": 3.12499999999998,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 2e-1
            self.atol = 1e-4
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Material properties
            rho = 2500.0  # density kg/m^3
            E = 70e9  # Young's modulus (Pa)
            nu = 0.3  # Poisson's ratio
            ys = 464.0e6  # yield stress

            # Plate geometry
            tplate = 0.005  # 5 mm

            # Set up property model
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set up constitutive model
            con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dv_num)
            transform = None
            # Set up element
            elem = elements.Quad4Shell(transform, con)
            scale = [100.0]
            return elem, scale

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        tacs_probs = []

        # Create two tacs problems with arbitrary loads
        # Static problem
        sp = fea_assembler.createStaticProblem(name="zero_load_static")
        F = np.array([0.0, 0.0, 1e6, 0.0, 0.0, 0.0])
        compIDs = fea_assembler.selectCompIDs(include="PLATE.00")
        sp.addLoadToComponents(compIDs, F)
        P = 100e5  # Pa
        compIDs = fea_assembler.selectCompIDs(include="PLATE.01")
        sp.addPressureToComponents(compIDs, P)
        tacs_probs.append(sp)

        # Transient problem
        tp = fea_assembler.createTransientProblem(
            name="zero_load_transient", tInit=0.0, tFinal=1.0, numSteps=5
        )
        time_steps = tp.getTimeSteps()
        for time_index, t in enumerate(time_steps):
            all_comps = fea_assembler.selectCompIDs()
            tp.addLoadToComponents(time_index, all_comps, t * F)
            tp.addPressureToComponents(time_index, all_comps, t * P)
        tacs_probs.append(tp)

        # Now loop through the problems we just created and zero out the loads
        # We really just want to test zero loads works as expected
        for prob in tacs_probs:
            prob.zeroLoads()

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction(
                "ks_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[100.0, 100.0, 100.0],
            )
            # Evaluate compliance of entire plate
            problem.addFunction("compliance", functions.Compliance)

        return tacs_probs, fea_assembler
