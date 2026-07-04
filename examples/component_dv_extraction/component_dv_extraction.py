"""
This example demonstrates how to extract a complete description of a model's sizing variables with pyTACS and use it to recreate the same sizing state in a separate TACS execution.

The model is a flat plate split into four components (PLATE.00 through PLATE.03), each with its own thickness.
Two components have active thickness design variables, and two are held fixed.
We perturb the active design variables (standing in for the result of an optimization), extract the component design variable dictionary, save it to disk with pickle, and then rebuild a second, identical model from the file.
"""

# [docs:imports-start]
import os
import pickle

import numpy as np

from tacs import constitutive, elements, functions, pyTACS

# [docs:imports-end]

# [docs:parameters-start]
bdf_file = os.path.join(os.path.dirname(__file__), "partitioned_plate.bdf")
sizing_file = os.path.join(os.path.dirname(__file__), "sizing.pkl")

# Constructor thickness for each component; only the components listed in
# active_components are given active thickness design variables
thicknesses = {
    "PLATE.00": 0.010,
    "PLATE.01": 0.012,
    "PLATE.02": 0.014,
    "PLATE.03": 0.016,
}
active_components = ["PLATE.00", "PLATE.02"]
# [docs:parameters-end]


# [docs:element-callback-start]
def element_callback(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    prop = constitutive.MaterialProperties(rho=2780.0, E=73.1e9, nu=0.33, ys=324e6)
    if compDescript in active_components:
        tNum = dvNum
    else:
        tNum = -1
    con = constitutive.IsoShellConstitutive(
        prop, t=thicknesses[compDescript], tNum=tNum
    )
    return elements.Quad4Shell(None, con)


# [docs:element-callback-end]

# [docs:setup-start]
FEAAssembler = pyTACS(bdf_file)
FEAAssembler.initialize(element_callback)

problem = FEAAssembler.createStaticProblem("gravity")
problem.addInertialLoad(np.array([0.0, 0.0, -9.81]))
problem.addFunction("mass", functions.StructuralMass)
problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=100.0)
# [docs:setup-end]

# [docs:perturb-start]
# Perturb the active design variables, standing in for the result of an optimization
x = problem.getDesignVars()
problem.setDesignVars(1.5 * x)
problem.solve()
funcs = {}
problem.evalFunctions(funcs)
# [docs:perturb-end]

# [docs:extract-start]
# Extract the sizing state. Every component appears in the dictionary, and every design
# variable group is included whether its entries are active design variables or not.
# Active entries hold the assembler's current values; inactive entries hold the values
# the constitutive objects were constructed with.
componentDVs = FEAAssembler.getComponentDesignVars()
if FEAAssembler.comm.rank == 0:
    print("Extracted component design variables:")
    for compDescript, groups in componentDVs.items():
        print(f"  {compDescript}: {groups}")
    with open(sizing_file, "wb") as f:
        pickle.dump(componentDVs, f)
# Make sure the file is written before any other processor tries to read it
FEAAssembler.comm.barrier()
# [docs:extract-end]

# [docs:restore-start]
# In a separate TACS execution, load the sizing file and use it inside the element
# callback. The group names and value types match the constitutive constructor keyword
# arguments, so the values can be passed straight through.
with open(sizing_file, "rb") as f:
    sizing = pickle.load(f)


def element_callback_restored(
    dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
):
    prop = constitutive.MaterialProperties(rho=2780.0, E=73.1e9, nu=0.33, ys=324e6)
    con = constitutive.IsoShellConstitutive(
        prop, t=sizing[compDescript]["t"], tNum=dvNum
    )
    return elements.Quad4Shell(None, con)


FEAAssembler2 = pyTACS(bdf_file)
FEAAssembler2.initialize(element_callback_restored)
problem2 = FEAAssembler2.createStaticProblem("gravity")
problem2.addInertialLoad(np.array([0.0, 0.0, -9.81]))
problem2.addFunction("mass", functions.StructuralMass)
problem2.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=100.0)
problem2.solve()
funcs2 = {}
problem2.evalFunctions(funcs2)

if FEAAssembler2.comm.rank == 0:
    print("Original model functions: ", funcs)
    print("Restored model functions: ", funcs2)
# [docs:restore-end]
