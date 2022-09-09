import numpy as np
from tacs import constitutive, elements

# Material properties
rho = 2780.0  # density kg/m^3
E = 73.1e9  # Young's modulus (Pa)
nu = 0.33  # Poisson's ratio
ys = 324.0e6  # yield stress

# Shell thickness
t = 0.005  # m
tMin = 0.002  # m
tMax = 0.05  # m

# Callback function used to setup TACS element objects and DVs
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []

    if (
        "MAIN_FLOOR" in compDescript
        or "CARGO_FLOOR" in compDescript
        or "SKIN" in compDescript
    ):
        ref_axis = np.array([0.0, 0.0, 1.0])
    else:  # Bulkheads, floor beams, and frames
        ref_axis = np.array([0.0, 1.0, 0.0])
    transform = elements.ShellRefAxisTransform(ref_axis)

    for elemDescript in elemDescripts:
        if elemDescript in ["CQUAD4", "CQUADR"]:
            elem = elements.Quad4Shell(transform, con)
        elif elemDescript in ["CTRIA3", "CTRIAR"]:
            elem = elements.Tri3Shell(transform, con)
        else:
            print("Uh oh, '%s' not recognized" % (elemDescript))
        elemList.append(elem)

    # Add scale for thickness dv
    scale = [100.0]
    return elemList, scale
