"""
This tests TACS structural optimization capabilities.
The beam model that we will be using for this problem is a rectangular beam,
cantilevered, with a shear load applied at the tip. The beam is discretized using
160 shell elements along it's span and depth.

This tests the MACH StructProblem object's DVGeo and design variable sensitivities.
"""

import os
import numpy as np
from mpi4py import MPI
import unittest

from tacs import pyTACS
from tacs import elements, constitutive, functions
from tacs.mach import StructProblem
from mach_struct_problem_base_test import MACHStructProblemTestCase

try:
    from pygeo import DVGeometry
except ImportError:
    DVGeometry = None

# Get file paths
base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/coarse_beam.bdf")
ffd_file = os.path.join(base_dir, "./input_files/ffd_8_linear.fmt")


@unittest.skipIf(DVGeometry is None, "pygeo is not installed")
class TestMACHBeamExample(MACHStructProblemTestCase.MACHStructProblemTest):
    """
    Test case for MACH StructProblem using the beam shape optimization example.
    """

    N_PROCS = 2

    # Reference values for regression testing
    FUNC_REFS = {
        "tip_shear_ks_vmfailure": 2.492124644882399,
        "tip_shear_mass": 2.7799999999999963,
    }

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

        # Shear force applied at tip
        V = 2.5e4

        # Callback function used to setup TACS element objects and DVs
        def element_callback(
            dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoShellConstitutive(prop, t=t, tNum=0)
            # TACS shells are sometimes a little overly-rigid in shear
            # We can reduce this effect by decreasing the drilling regularization
            con.setDrillingRegularization(0.1)
            refAxis = np.array([1.0, 0.0, 0.0])
            transform = elements.ShellRefAxisTransform(refAxis)
            elem = elements.Quad4Shell(transform, con)
            return elem

        # Create pyTACS assembler
        FEAAssembler = pyTACS(bdf_file)
        FEAAssembler.initialize(element_callback)

        # Create DVGeometry
        DVGeo = DVGeometry(fileName=ffd_file, isComplex=self.dtype == complex)
        # Create reference axis
        nRefAxPts = DVGeo.addRefAxis(name="centerline", alignIndex="i", yFraction=0.5)

        # Set up global design variables
        def depth(val, geo):
            for i in range(nRefAxPts):
                geo.scale_y["centerline"].coef[i] = val[i]

        DVGeo.addGlobalDV(
            dvName="depth",
            value=np.ones(nRefAxPts),
            func=depth,
            lower=1e-3,
            upper=10.0,
            scale=20.0,
        )

        # Create static problem
        staticProb = FEAAssembler.createStaticProblem("tip_shear")
        # Add TACS Functions
        staticProb.addFunction("mass", functions.StructuralMass)
        staticProb.addFunction(
            "ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=10.0
        )
        # Add forces to static problem
        staticProb.addLoadToNodes(
            206, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True
        )
        # Set convergence to be tight for test
        staticProb.setOption("L2Convergence", 1e-20)
        staticProb.setOption("L2ConvergenceRel", 1e-20)

        # Create MACH StructProblem
        structProb = StructProblem(staticProb, FEAAssembler, DVGeo=DVGeo)

        return [structProb]
