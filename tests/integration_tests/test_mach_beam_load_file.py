"""
This tests TACS structural optimization capabilities.
The beam model that we will be using for this problem is a rectangular beam,
cantilevered, with a shear load applied at the tip. The beam is discretized using
160 shell elements along it's span and depth.

This tests the MACH StructProblem object's ability to load forces from a file.
The problem is identical to test_mach_beam_dvgeo, so reference values should match that test.
"""

import os
import numpy as np
from mpi4py import MPI
import unittest

from tacs import pyTACS
from tacs import elements, constitutive, functions
from tacs.mach import StructProblem
from mach_struct_problem_base_test import MACHStructProblemTestCase

# Get file paths
base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/coarse_beam.bdf")
load_file = os.path.join(base_dir, "./input_files/beam_load.dat")

from test_mach_beam_dvgeo import ksweight, TestMACHBeamExample as TMBE


class TestMACHBeamExample(MACHStructProblemTestCase.MACHStructProblemTest):
    """
    Test case for MACH StructProblem using the beam shape optimization example.
    """

    N_PROCS = 2

    # Set reference functions to match unrotated case
    FUNC_REFS = TMBE.FUNC_REFS

    def setup_struct_problems(self, comm):
        """
        Setup MACH StructProblem objects for testing.
        """

        # Beam properties
        t = 0.01  # m

        # Material properties
        rho = 2780.0  # kg /m^3
        E = 70.0e9
        nu = 0.0
        ys = 420.0e6

        # Callback function used to setup TACS element objects and DVs
        def element_callback(
            dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoShellConstitutive(prop, t=t, tNum=0)
            refAxis = np.array([1.0, 0.0, 0.0])
            transform = elements.ShellRefAxisTransform(refAxis)
            elem = elements.Quad4Shell(transform, con)
            return elem

        # Create pyTACS assembler
        FEAAssembler = pyTACS(bdf_file)
        FEAAssembler.initialize(element_callback)

        # Create static problem
        staticProb = FEAAssembler.createStaticProblem("tip_shear")
        # Add TACS Functions
        staticProb.addFunction("mass", functions.StructuralMass)
        staticProb.addFunction(
            "ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=ksweight
        )
        # Set convergence to be tight for test
        staticProb.setOption("L2Convergence", 1e-20)
        staticProb.setOption("L2ConvergenceRel", 1e-20)

        # Create MACH StructProblem
        structProb = StructProblem(staticProb, FEAAssembler, loadFile=load_file)

        return [structProb]
