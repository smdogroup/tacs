"""
==============================================================================

==============================================================================
@File    :   test_multi_beam_frame.py
@Date    :   2025/10/31
@Author  :   Alasdair Christison Gray
@Description : This test involves a simple frame structure made of multiple beams as shown below:

    (5)        (3)       (6)
6 ------- 3 -------- 4 ------- 7
          |          |
          |          |
      (2) |          |
          |          |
          2          | (4)
          |          |
          |          |
      (1) |          |
          |          |
          1          5

The elements are grouped into components as follows:
- Component 0: Elements 1 and 2
- Component 1: Elements 3 and 6
- Component 2: Element 4
- Component 3: Element 5

Nodes 1 and 5 are fixed in all degrees of freedom, and forces are applied at nodes 6 and 7

This test is mostly designed to verify that adjacency constraints are correctly
set up across multiple beams with an unstructured connectivity.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os
import unittest

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, TACS, functions, elements, constitutive

TACS_IS_COMPLEX = TACS.dtype == complex

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/multi_beam_frame.bdf")

ksweight = 10.0

# Material properties
rho = 2700.0  # density kg/m^3
E = 70.0e9  # Young's modulus (Pa)
nu = 0.3  # Poisson's ratio
ys = 300e6  # yield stress


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "load_set_001_compliance": 7.418744658419356,
        "load_set_001_ks_vmfailure": 0.22143528124961626,
        "load_set_001_mass": 23.625000000000004,
        "load_set_001_x_disp": 0.1259149590010181,
        "load_set_001_y_disp": 0.12512170175459475,
        "load_set_001_z_disp": 0.10296855090569572,
    }

    def setup_tacs_problems(self, comm):
        FEAssembler = pytacs.pyTACS(bdf_file, comm)

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            t = 25e-3  # m
            w = 0.1  # m
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(
                prop, t=t, tNum=dv_num, w=w, wNum=dv_num + 1
            )
            refAxis = np.array([0.0, 0.0, 1.0])
            transform = elements.BeamRefAxisTransform(refAxis)
            # pass back the appropriate tacs element object
            elem = elements.Beam2(transform, con)
            return elem

        FEAssembler.initialize(elemCallBack=elem_call_back)

        tacsProbs = FEAssembler.createTACSProbsFromBDF()
        tacsProbs = list(tacsProbs.values())

        # Set convergence to be tight for test
        for problem in tacsProbs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacsProbs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction(
                "x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[10.0, 0.0, 0.0],
            )
            problem.addFunction(
                "y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 10.0, 0.0],
            )
            problem.addFunction(
                "z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, 10.0],
            )

        # Add adjacency constraints
        constr = FEAssembler.createAdjacencyConstraint("adjcon")
        allComps = FEAssembler.selectCompIDs()
        constr.addConstraint(
            "beamThicknessAdjCon",
            dvIndex=0,
            compIDs=allComps,
        )
        constr.addConstraint(
            "beamThicknessWidthCon",
            dvIndex=1,
            compIDs=allComps,
        )
        tacsProbs.append(constr)

        return tacsProbs, FEAssembler

    def test_adjacency(self):
        adjCon = self.tacs_probs[-1]
        expectedAdjacency = [(0, 1), (0, 3), (1, 2), (1, 3)]
        self.assertEqual(sorted(expectedAdjacency), sorted(adjCon.adjacentComps))


if __name__ == "__main__":
    unittest.main()
