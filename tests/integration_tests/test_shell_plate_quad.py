import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
The nominal case is a 1m x 1m flat plate under three load cases:
a 10 kN point force at center, a 100kPa pressure applied to the surface, and a 100G gravity load. The
perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test KSFailure, StructuralMass, CenterOfMass, MomentOfInertia,
and Compliance functions and sensitivities. Finally, a modal analysis is performed
and the eigenvalues tested.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/plate.bdf")

# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "point_load_Ixx": 1.0416927083333256,
        "point_load_Ixy": 5.329070518200751e-15,
        "point_load_Ixz": 0.0,
        "point_load_Iyy": 1.0416927083333274,
        "point_load_Iyz": 0.0,
        "point_load_Izz": 2.0833333333333224,
        "point_load_cgx": 0.5000000000000003,
        "point_load_cgy": 0.5000000000000007,
        "point_load_cgz": 0.0,
        "point_load_compliance": 683.8571611645635,
        "point_load_ks_vmfailure": 0.5757488025912149,
        "point_load_ks_vmfailure_discrete": 1.1748952573020133,
        "point_load_ks_vmfailure_discrete_average": 1.0273708117809457,
        "point_load_mass": 12.5,
        "pressure_Ixx": 1.0416927083333256,
        "pressure_Ixy": 5.329070518200751e-15,
        "pressure_Ixz": 0.0,
        "pressure_Iyy": 1.0416927083333274,
        "pressure_Iyz": 0.0,
        "pressure_Izz": 2.0833333333333224,
        "pressure_cgx": 0.5000000000000003,
        "pressure_cgy": 0.5000000000000007,
        "pressure_cgz": 0.0,
        "pressure_compliance": 4679.345460332448,
        "pressure_ks_vmfailure": 1.2938623156883402,
        "pressure_ks_vmfailure_discrete": 1.8930087703991383,
        "pressure_ks_vmfailure_discrete_average": 1.7458488121120046,
        "pressure_mass": 12.5,
        "gravity_Ixx": 1.0416927083333256,
        "gravity_Ixy": 5.329070518200751e-15,
        "gravity_Ixz": 0.0,
        "gravity_Iyy": 1.0416927083333274,
        "gravity_Iyz": 0.0,
        "gravity_Izz": 2.0833333333333224,
        "gravity_cgx": 0.5000000000000003,
        "gravity_cgy": 0.5000000000000007,
        "gravity_cgz": 0.0,
        "gravity_compliance": 70.36280588360437,
        "gravity_ks_vmfailure": 0.11707320009756256,
        "gravity_ks_vmfailure_discrete": 0.7162196548083607,
        "gravity_ks_vmfailure_discrete_average": 0.5772792866164124,
        "gravity_mass": 12.5,
        "modal_eigsm.0": 87437.50645925231,
        "modal_eigsm.1": 396969.8881662692,
        "modal_eigsm.2": 396969.8881667623,
        "modal_eigsm.3": 866727.671481833,
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
        struct_options = {
            "writeCoordinateFrame": True,
        }

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

        # Add point force to node 81 (center of plate)
        sp = fea_assembler.createStaticProblem(name="point_load")
        F = np.array([0.0, 0.0, 1e4, 0.0, 0.0, 0.0])
        sp.addLoadToNodes(81, F, nastranOrdering=True)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name="pressure")
        P = 100e3  # Pa
        compIDs = fea_assembler.selectCompIDs(include="PLATE")
        sp.addPressureToComponents(compIDs, P)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name="gravity")
        g = np.array([0.0, 0.0, -981.0], dtype=self.dtype)
        sp.addInertialLoad(g)
        tacs_probs.append(sp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction(
                "ks_vmfailure_discrete",
                functions.KSFailure,
                ksWeight=ksweight,
                ftype="discrete",
            )
            problem.addFunction(
                "ks_vmfailure_discrete_average",
                functions.KSFailure,
                ksWeight=ksweight,
                ftype="discrete-average",
            )
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "cgx", functions.CenterOfMass, direction=[1.0, 0.0, 0.0]
            )
            problem.addFunction(
                "cgy", functions.CenterOfMass, direction=[0.0, 1.0, 0.0]
            )
            problem.addFunction(
                "cgz", functions.CenterOfMass, direction=[0.0, 0.0, 1.0]
            )
            problem.addFunction(
                "Ixx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Ixy",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Ixz",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Iyy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Iyz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Izz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )

        # Add modal problem
        mp = fea_assembler.createModalProblem("modal", 2e5, 10)
        mp.setOption("printLevel", 2)
        tacs_probs.append(mp)

        return tacs_probs, fea_assembler
