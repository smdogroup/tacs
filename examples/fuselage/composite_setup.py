import numpy as np
from tacs import constitutive, elements

# Create the quasi-isotropic layup with 48 plies
rho = 1550.0
specific_heat = 921.096
E1 = 54e9
E2 = 18e9
nu12 = 0.25
G12 = 9e9
G13 = 9e9
Xt = 2410.0e6
Xc = 1040.0e6
Yt = 73.0e6
Yc = 173.0e6
S12 = 71.0e6
cte = 24.0e-6
kappa = 230.0

ply_thickness = 0.125e-3
nplies = 48

# Quasi-iso layup
# NOTE: Angles must be in radians
qi_angles = np.deg2rad(np.array([0.0, -45.0, 90.0, 45.0]))
nrepeats = nplies // len(qi_angles)

# Callback function used to setup TACS element objects and DVs
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup (orthotropic) property and constitutive objects
    ortho_prop = constitutive.MaterialProperties(
        rho=rho,
        specific_heat=specific_heat,
        E1=E1,
        E2=E2,
        nu12=nu12,
        G12=G12,
        G13=G13,
        G23=G13,
        Xt=Xt,
        Xc=Xc,
        Yt=Yt,
        Yc=Yc,
        S12=S12,
        cte=cte,
        kappa=kappa,
    )

    ortho_ply = constitutive.OrthotropicPly(ply_thickness, ortho_prop)

    # Setup layup properties
    # copy each list until we have entry for every ply
    ortho_layup = [ortho_ply] * nplies
    ply_thicknesses = np.array([ply_thickness] * nplies)
    ply_angles = np.tile(qi_angles, nrepeats)

    # Setup composite constitutive object
    # NOTE: The CompositeShellConstitutive has no design variables
    con = constitutive.CompositeShellConstitutive(
        ortho_layup, ply_thicknesses, ply_angles
    )

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []

    # Set reference axis for plies, this determines what direction aligns with 0 degrees for the laminate
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

    return elemList
