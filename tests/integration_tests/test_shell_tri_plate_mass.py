import os
import unittest

from mpi4py import MPI

from tacs import constitutive, elements, functions, pytacs

"""
The 4m x 3m x 20mm flat plate from test_shell_tri_buckling_axial.py, meshed with
three-node triangles. A constant-thickness flat plate's mass is integrated exactly by
any correct quadrature rule, so every Tri3 shell element type must give rho * t * area.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/ss_plate_tri.bdf")

DENSITY = 2500.0  # kg/m^3
THICKNESS = 0.020  # m
PLATE_LENGTH = 4.0  # m
PLATE_WIDTH = 3.0  # m

TRI3_SHELL_CLASSES = [
    elements.Tri3Shell,
    elements.Tri3NonlinearShell,
    elements.Tri3ShellModRot,
    elements.Tri3NonlinearShellModRot,
    elements.Tri3ShellQuaternion,
    elements.Tri3NonlinearShellQuaternion,
    elements.Tri3ThermalShell,
    elements.Tri3NonlinearThermalShell,
]


def compute_plate_mass(element_class):
    """Compute the structural mass of the plate meshed with the given Tri3 shell element type"""
    fea_assembler = pytacs.pyTACS(bdf_file, MPI.COMM_WORLD, options={"printLevel": 0})

    def elem_call_back(
        dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
    ):
        prop = constitutive.MaterialProperties(rho=DENSITY, E=205e9, nu=0.3, ys=464.0e6)
        con = constitutive.IsoShellConstitutive(prop, t=THICKNESS, tNum=dv_num)
        return element_class(None, con)

    fea_assembler.initialize(elem_call_back)
    problem = fea_assembler.createStaticProblem("plate")
    problem.addFunction("mass", functions.StructuralMass)
    funcs = {}
    problem.evalFunctions(funcs)
    return funcs["plate_mass"]


class Tri3ShellPlateMassTest(unittest.TestCase):
    def test_plate_mass(self):
        """Assert that every Tri3 shell element type computes the exact mass of a flat plate"""
        exact_mass = DENSITY * THICKNESS * PLATE_LENGTH * PLATE_WIDTH
        for element_class in TRI3_SHELL_CLASSES:
            with self.subTest(element=element_class.__name__):
                mass = compute_plate_mass(element_class)
                self.assertAlmostEqual(mass, exact_mass, delta=1e-10 * exact_mass)


if __name__ == "__main__":
    unittest.main()
