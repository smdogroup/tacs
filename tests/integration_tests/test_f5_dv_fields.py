"""Tests for the named design variable fields written to F5 output files."""

import os
import shutil
import tempfile
import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements, pyTACS

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARTITIONED_PLATE_BDF = os.path.join(BASE_DIR, "input_files", "partitioned_plate.bdf")

# Components PLATE.00/PLATE.02 get an isotropic shell (scalar group "t"), PLATE.01 gets a smeared composite shell (scalar "thickness" plus array "ply_fractions"), and PLATE.03 gets a blade-stiffened shell with inactive DVs and derived outputs (e.g "effectiveThickness").
# The union contains columns that are NaN on some components and real on others.
ISO_COMPONENTS = ["PLATE.00", "PLATE.02"]
ISO_THICKNESSES = {"PLATE.00": 0.010, "PLATE.02": 0.014}
# PLATE.02 deliberately has NO active DV: inactive groups must still be written
ACTIVE_ISO_COMPONENTS = ["PLATE.00"]

BLADE_COMPONENT = "PLATE.03"
BLADE_PANEL_LENGTH = 0.5
BLADE_STIFFENER_PITCH = 0.2
BLADE_PANEL_THICK = 0.02
BLADE_STIFFENER_HEIGHT = 0.05
BLADE_STIFFENER_THICK = 0.008
BLADE_FLANGE_FRACTION = 0.8


def makeIsoCon(t, tNum):
    prop = constitutive.MaterialProperties(rho=2780.0, E=73.1e9, nu=0.33, ys=324e6)
    return constitutive.IsoShellConstitutive(prop, t=t, tNum=tNum)


def makeSmearedCon(dvNum):
    prop = constitutive.MaterialProperties(
        rho=1550.0,
        specific_heat=921.096,
        E1=54e9,
        E2=18e9,
        nu12=0.25,
        G12=9e9,
        G13=9e9,
        G23=9e9,
        Xt=2410.0e6,
        Xc=1040.0e6,
        Yt=73.0e6,
        Yc=173.0e6,
        S12=71.0e6,
        alpha=24.0e-6,
        kappa=230.0,
    )
    ply = constitutive.OrthotropicPly(0.1, prop)
    plyAngles = np.deg2rad([0.0, 45.0, 90.0]).astype(TACS.dtype)
    plyFractions = np.array([0.5, 0.3, 0.2], dtype=TACS.dtype)
    plyFractionDVNums = np.array([dvNum + 1, dvNum + 2, dvNum + 3], dtype=np.intc)
    return constitutive.SmearedCompositeShellConstitutive(
        [ply] * 3,
        0.01,
        plyAngles,
        plyFractions,
        thickness_dv_num=dvNum,
        ply_fraction_dv_nums=plyFractionDVNums,
    )


def makeBladeCon():
    prop = constitutive.MaterialProperties(
        rho=1550.0,
        E1=54e9,
        E2=18e9,
        nu12=0.25,
        G12=9e9,
        G13=9e9,
        G23=9e9,
        Xt=2410.0e6,
        Xc=1040.0e6,
        Yt=73.0e6,
        Yc=173.0e6,
        S12=71.0e6,
    )
    ply = constitutive.OrthotropicPly(1e-3, prop)
    panelPlyAngles = np.deg2rad([0.0, 45.0, 90.0]).astype(TACS.dtype)
    panelPlyFracs = np.array([0.5, 0.3, 0.2], dtype=TACS.dtype)
    stiffenerPlyAngles = np.deg2rad([0.0, 60.0]).astype(TACS.dtype)
    stiffenerPlyFracs = np.array([0.7, 0.3], dtype=TACS.dtype)
    # All DVs left inactive: groups and derived outputs must be written anyway
    return constitutive.BladeStiffenedShellConstitutive(
        ply,
        ply,
        BLADE_PANEL_LENGTH,
        BLADE_STIFFENER_PITCH,
        BLADE_PANEL_THICK,
        panelPlyAngles,
        panelPlyFracs,
        BLADE_STIFFENER_HEIGHT,
        BLADE_STIFFENER_THICK,
        stiffenerPlyAngles,
        stiffenerPlyFracs,
        flangeFraction=BLADE_FLANGE_FRACTION,
    )


def entryValues(groupValues):
    """Flatten a component's group-value dict into {entry name: float}.
    Mirrors the C++ naming rule: scalars keep the group name, arrays get _i suffixes.
    """
    entries = {}
    for name, value in groupValues.items():
        if isinstance(value, np.ndarray):
            for i in range(len(value)):
                entries[f"{name}_{i}"] = float(np.real(value[i]))
        else:
            entries[name] = float(np.real(value))
    return entries


def expectedEntries(fea, consByComp, includeGroups=True, includeDerived=True):
    """Build {component: {field name: value}} the way the writer does."""
    entries = {descript: {} for descript in consByComp}
    if includeGroups:
        # getComponentDesignVars is now a problem-level method; a throwaway
        # problem reads back the same (unmodified) design vars this test uses.
        problem = fea.createStaticProblem("expectedEntries")
        compDVs = problem.getComponentDesignVars()
        for descript, groups in compDVs.items():
            entries[descript].update(entryValues(groups))
    if includeDerived:
        for descript, con in consByComp.items():
            for name, value in con.getDerivedOutputs().items():
                entries[descript][name] = float(np.real(value))
    return entries


class F5DesignVarFieldTest(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    def setupMixedPlate(self, writeDesignVars=True, writeDerivedOutputs=True):
        fea = pyTACS(
            PARTITIONED_PLATE_BDF,
            options={
                "writeDesignVars": writeDesignVars,
                "writeDerivedOutputs": writeDerivedOutputs,
            },
            comm=self.comm,
        )
        consByComp = {}

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            if compDescript in ISO_COMPONENTS:
                if compDescript in ACTIVE_ISO_COMPONENTS:
                    tNum = dvNum
                else:
                    tNum = -1
                con = makeIsoCon(ISO_THICKNESSES[compDescript], tNum)
            elif compDescript == BLADE_COMPONENT:
                con = makeBladeCon()
            else:
                con = makeSmearedCon(dvNum)
            consByComp[compDescript] = con
            return elements.Quad4Shell(None, con)

        fea.initialize(elemCallBack)
        return fea, consByComp

    def writeAndLoad(self, fea):
        """Write an f5 on all ranks, read the dv zone back on the root rank, and broadcast the results so every rank can run identical assertions."""
        if self.comm.rank == 0:
            tmpDir = tempfile.mkdtemp()
        else:
            tmpDir = None
        tmpDir = self.comm.bcast(tmpDir, root=0)

        problem = fea.createStaticProblem("dvfields")
        # Pass an explicit solution number because the problem is never solved, so the internal call counter is still -1.
        problem.writeSolution(outputDir=tmpDir, baseName="dvfields", number=0)
        self.comm.barrier()

        if self.comm.rank == 0:
            loader = TACS.FH5Loader()
            loader.loadData(os.path.join(tmpDir, "dvfields_000.f5"))
            names, data = loader.getDesignData()
            compNums, _, _, _ = loader.getConnectivity()
            compNames = [
                loader.getComponentName(k) for k in range(loader.getNumComponents())
            ]
            payload = (names, np.array(data), np.array(compNums), compNames)
        else:
            payload = None
        payload = self.comm.bcast(payload, root=0)

        self.comm.barrier()
        if self.comm.rank == 0:
            shutil.rmtree(tmpDir)
        return payload

    def test_dv_field_names_and_values(self):
        fea, consByComp = self.setupMixedPlate()
        names, data, compNums, compNames = self.writeAndLoad(fea)

        # The expected union comes from the (rank-identical) component DV dict
        # and the constitutive derived outputs, flattened with the same naming
        # rule the writer uses
        compEntries = expectedEntries(fea, consByComp)
        expectedUnion = sorted(set().union(*[e.keys() for e in compEntries.values()]))
        self.assertEqual(names.split(","), expectedUnion)

        self.assertEqual(data.shape, (len(compNums), len(expectedUnion)))

        # Every element row must hold its component's values, NaN elsewhere
        for e in range(len(compNums)):
            descript = compNames[compNums[e]]
            expected = compEntries[descript]
            for j, name in enumerate(expectedUnion):
                if name in expected:
                    np.testing.assert_allclose(
                        data[e, j],
                        expected[name],
                        rtol=1e-12,
                        err_msg=f"element {e} ({descript}), field {name}",
                    )
                else:
                    self.assertTrue(
                        np.isnan(data[e, j]),
                        msg=f"element {e} ({descript}), field {name} should be NaN",
                    )

        # Spot-check the blade component's effective thickness against the formula
        effTCol = expectedUnion.index("effectiveThickness")
        stiffenerArea = (
            (1.0 + BLADE_FLANGE_FRACTION)
            * BLADE_STIFFENER_HEIGHT
            * BLADE_STIFFENER_THICK
        )
        expectedEffT = BLADE_PANEL_THICK + stiffenerArea / BLADE_STIFFENER_PITCH
        bladeRows = [
            e for e in range(len(compNums)) if compNames[compNums[e]] == BLADE_COMPONENT
        ]
        self.assertGreater(len(bladeRows), 0)
        for e in bladeRows:
            np.testing.assert_allclose(data[e, effTCol], expectedEffT, rtol=1e-12)
        # Non-blade components must hold NaN in the derived column
        for e in range(len(compNums)):
            if compNames[compNums[e]] != BLADE_COMPONENT:
                self.assertTrue(np.isnan(data[e, effTCol]))

    def test_inactive_groups_are_written(self):
        fea, _ = self.setupMixedPlate()
        names, data, compNums, compNames = self.writeAndLoad(fea)
        nameList = names.split(",")
        tCol = nameList.index("t")
        # PLATE.02 has tNum=-1, but its thickness must still appear in the file
        inactiveRows = [
            e for e in range(len(compNums)) if compNames[compNums[e]] == "PLATE.02"
        ]
        self.assertGreater(len(inactiveRows), 0)
        for e in inactiveRows:
            np.testing.assert_allclose(data[e, tCol], 0.014, rtol=1e-12)

    def test_derived_outputs_only(self):
        fea, consByComp = self.setupMixedPlate(writeDesignVars=False)
        names, data, compNums, compNames = self.writeAndLoad(fea)
        expected = expectedEntries(fea, consByComp, includeGroups=False)
        expectedUnion = sorted(set().union(*[e.keys() for e in expected.values()]))
        self.assertEqual(names.split(","), expectedUnion)
        self.assertNotIn("t", names.split(","))

    def test_design_vars_only(self):
        fea, consByComp = self.setupMixedPlate(writeDerivedOutputs=False)
        names, data, compNums, compNames = self.writeAndLoad(fea)
        expected = expectedEntries(fea, consByComp, includeDerived=False)
        expectedUnion = sorted(set().union(*[e.keys() for e in expected.values()]))
        self.assertEqual(names.split(","), expectedUnion)
        self.assertNotIn("effectiveThickness", names.split(","))

    def test_no_dv_zone_when_disabled(self):
        fea, _ = self.setupMixedPlate(writeDesignVars=False, writeDerivedOutputs=False)
        if self.comm.rank == 0:
            tmpDir = tempfile.mkdtemp()
        else:
            tmpDir = None
        tmpDir = self.comm.bcast(tmpDir, root=0)

        problem = fea.createStaticProblem("dvfields")
        # Pass an explicit solution number because the problem is never solved, so the internal call counter is still -1.
        problem.writeSolution(outputDir=tmpDir, baseName="dvfields", number=0)
        self.comm.barrier()

        raised = False
        if self.comm.rank == 0:
            loader = TACS.FH5Loader()
            loader.loadData(os.path.join(tmpDir, "dvfields_000.f5"))
            try:
                loader.getDesignData()
            except ValueError:
                raised = True
        raised = self.comm.bcast(raised, root=0)

        self.comm.barrier()
        if self.comm.rank == 0:
            shutil.rmtree(tmpDir)
        self.assertTrue(raised)


if __name__ == "__main__":
    unittest.main()
