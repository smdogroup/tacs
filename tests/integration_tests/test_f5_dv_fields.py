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

# Components PLATE.00/PLATE.02 get an isotropic shell (scalar group "t"), the other two
# get a smeared composite shell (scalar "thickness" plus array "ply_fractions"), so the
# union contains columns that are NaN on some components and real on others.
ISO_COMPONENTS = ["PLATE.00", "PLATE.02"]
ISO_THICKNESSES = {"PLATE.00": 0.010, "PLATE.02": 0.014}
# PLATE.02 deliberately has NO active DV: inactive groups must still be written
ACTIVE_ISO_COMPONENTS = ["PLATE.00"]


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


class F5DesignVarFieldTest(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    def setupMixedPlate(self, writeDesignVars=True):
        fea = pyTACS(
            PARTITIONED_PLATE_BDF,
            options={"writeDesignVars": writeDesignVars},
            comm=self.comm,
        )

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            if compDescript in ISO_COMPONENTS:
                if compDescript in ACTIVE_ISO_COMPONENTS:
                    tNum = dvNum
                else:
                    tNum = -1
                con = makeIsoCon(ISO_THICKNESSES[compDescript], tNum)
            else:
                con = makeSmearedCon(dvNum)
            return elements.Quad4Shell(None, con)

        fea.initialize(elemCallBack)
        return fea

    def writeAndLoad(self, fea):
        """Write an f5 on all ranks, read the dv zone back on the root rank, and
        broadcast the results so every rank can run identical assertions.
        """
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
            names, data = loader.getDesignVarData()
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
        fea = self.setupMixedPlate()
        names, data, compNums, compNames = self.writeAndLoad(fea)

        # The expected union comes from the (rank-identical) component DV dict,
        # flattened with the same naming rule the writer uses
        compDVs = fea.getComponentDesignVars()
        compEntries = {
            descript: entryValues(groups) for descript, groups in compDVs.items()
        }
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

    def test_inactive_groups_are_written(self):
        fea = self.setupMixedPlate()
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

    def test_no_dv_zone_when_disabled(self):
        fea = self.setupMixedPlate(writeDesignVars=False)
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
                loader.getDesignVarData()
            except ValueError:
                raised = True
        raised = self.comm.bcast(raised, root=0)

        self.comm.barrier()
        if self.comm.rank == 0:
            shutil.rmtree(tmpDir)
        self.assertTrue(raised)


if __name__ == "__main__":
    unittest.main()
