"""
This tests TACS structural optimization capabilities.
The beam model that we will be using for this problem is a rectangular beam,
cantilevered, with a shear load applied at the tip. The beam is discretized using
160 shell elements along it's span and depth.

This tests the MACH StructProblem object's DVGeo and design variable sensitivities.
"""

import glob
import os
import shutil
import tempfile
import numpy as np
from mpi4py import MPI
import unittest

import pyNastran.bdf as pn

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
load_file = os.path.join(base_dir, "./input_files/beam_load_axial.dat")

ksweight = 10.0
# Shear force applied at tip
V = 2.5e4
# Pressure applied on face of beam
P = 1e5
# Gravity vector
G = np.array([0.0, -9.81, 0.0]) * 1000


@unittest.skipIf(DVGeometry is None, "pygeo is not installed")
class TestMACHBeamExample(MACHStructProblemTestCase.MACHStructProblemTest):
    """
    Test case for MACH StructProblem using the beam shape optimization example.
    """

    N_PROCS = 2

    # Reference values for regression testing
    FUNC_REFS = {
        "tip_shear_ks_vmfailure": 2.492124644887184,
        "tip_shear_mass": 2.7799999999999963,
        "pressure_ks_vmfailure": 6.37743474872165,
        "pressure_mass": 2.7799999999999963,
        "gravity_ks_vmfailure": 1.037787648560678,
        "gravity_mass": 2.7799999999999963,
        "combined_ks_vmfailure": 7.702671860961564,
        "combined_mass": 2.7799999999999963,
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

        staticProblems = []
        structProblems = []

        # Create static problem
        tipShear = FEAAssembler.createStaticProblem("tip_shear")
        # Add forces to static problem
        tipShear.addLoadToNodes(206, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)
        staticProblems.append(tipShear)

        # Create static problem
        pressure = FEAAssembler.createStaticProblem("pressure")
        # Add pressure to static problem
        allCompIDs = FEAAssembler.selectCompIDs()
        pressure.addPressureToComponents(allCompIDs, P)
        staticProblems.append(pressure)

        # Create static problem
        gravity = FEAAssembler.createStaticProblem("gravity")
        # Add gravity to static problem
        gravity.addInertialLoad(G)
        staticProblems.append(gravity)

        # Create static problem
        combined = FEAAssembler.createStaticProblem("combined")
        # Add forces to static problem
        combined.addLoadToNodes(206, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)
        # Add pressure to static problem
        combined.addPressureToComponents(allCompIDs, P)
        # Add gravity to static problem
        combined.addInertialLoad(G)
        staticProblems.append(combined)

        # Set convergence to be tight for test
        for problem in staticProblems:
            # Add TACS Functions
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction(
                "ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=ksweight
            )

            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

            # Create MACH StructProblem
            structProblems.append(StructProblem(problem, FEAAssembler, DVGeo=DVGeo))

        # Add external force file to combined problem that mimics aero load
        # should contain an axial force at tip in this case
        combinedStructProblem = structProblems[-1]
        combinedStructProblem.readExternalForceFile(load_file)

        return structProblems

    def test_write_external_force_file(self):
        """
        Test that writeExternalForceFile writes only the external (aero) loads,
        not the internal structural loads (shear, pressure, gravity). The combined
        problem has an axial force of 2.5e4 in the x direction set as its external
        force via readExternalForceFile. The written file should contain exactly one
        FORCE card with that value.
        """
        combined_prob = self.struct_probs[-1]

        # Create temp dir on rank 0 and broadcast so all MPI ranks share the same path
        if self.comm.rank == 0:
            tmpdir = tempfile.mkdtemp()
        else:
            tmpdir = None
        tmpdir = self.comm.bcast(tmpdir, root=0)

        try:
            combined_prob.writeExternalForceFile(outputDir=tmpdir)

            # Only rank 0 writes the file, so only rank 0 reads and validates it
            if self.comm.rank == 0:
                dat_files = glob.glob(os.path.join(tmpdir, "*.dat"))
                self.assertEqual(len(dat_files), 1)

                bdf = pn.bdf.read_bdf(
                    dat_files[0], validate=False, xref=False, debug=False, punch=True
                )

                all_forces = [card for cards in bdf.loads.values() for card in cards]
                self.assertEqual(len(all_forces), 1)

                force = all_forces[0]
                self.assertAlmostEqual(force.mag * force.xyz[0], V, places=1)
                self.assertAlmostEqual(force.mag * force.xyz[1], 0.0, places=1)
                self.assertAlmostEqual(force.mag * force.xyz[2], 0.0, places=1)
        finally:
            # Barrier ensures all ranks finish before rank 0 deletes the temp dir
            self.comm.barrier()
            if self.comm.rank == 0:
                shutil.rmtree(tmpdir)
