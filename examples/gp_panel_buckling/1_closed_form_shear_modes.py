"""
Check axial closed-form modes in GPBladeStiffenedShellConstitutive class are reasonable
@Author Sean Engelstad
@Date 05/16/2024
"""

import numpy as np
import matplotlib.pyplot as plt
import niceplots
from tacs import TACS, constitutive

DEG2RAD = np.pi / 180.0

dtype = TACS.dtype

# Create the orthotropic layup
ortho_prop = constitutive.MaterialProperties(
    rho=1550,
    specific_heat=921.096,
    E1=54e3,
    E2=18e3,
    nu12=0.25,
    G12=9e3,
    G13=9e3,
    G23=9e3,
    Xt=2410.0,
    Xc=1040.0,
    Yt=73.0,
    Yc=173.0,
    S12=71.0,
    alpha=24.0e-6,
    kappa=230.0,
)
ortho_ply = constitutive.OrthotropicPly(1e-3, ortho_prop)

# don't put in any GP models (so using closed-form solutions rn)
con = constitutive.GPBladeStiffenedShellConstitutive(
    panelPly=ortho_ply,
    stiffenerPly=ortho_ply,
    panelLength=2.0,
    stiffenerPitch=0.2,
    panelThick=1.5e-2,
    panelPlyAngles=np.array([0.0, 45.0, 90.0], dtype=dtype) * DEG2RAD,
    panelPlyFracs=np.array([0.5, 0.3, 0.2], dtype=dtype),
    stiffenerHeight=0.075,
    stiffenerThick=1e-2,
    stiffenerPlyAngles=np.array([0.0, 60.0], dtype=dtype) * DEG2RAD,
    stiffenerPlyFracs=np.array([0.6, 0.4], dtype=dtype),
    panelWidth=1.0,
    flangeFraction=0.8,
)
# Set the KS weight really low so that all failure modes make a
# significant contribution to the failure function derivatives
con.setKSWeight(20.0)

# get the axial loads in nondimensional space w.r.t. rho_0
n = 500
plt.style.use(niceplots.get_style())
rho0_vec = np.linspace(0.5, 10.0, n)
N12cr_vec = np.zeros((n,), dtype=dtype)
for gamma in [0.0, 0.1, 0.5, 1.0]:
    for i, rho0 in enumerate(rho0_vec):
        N12cr_vec[i] = con.nondimCriticalGlobalShearLoad(rho0, 0.4, gamma)
    plt.plot(rho0_vec, N12cr_vec, label=f"gamma={gamma:.2f}")

# plot it
plt.margins(x=0.05, y=0.05)
plt.xlabel(r"$\rho_0$")
plt.ylabel(r"$N_{12,cr}^*$")
plt.legend()
plt.savefig("1-verify-shear.png", dpi=400)
