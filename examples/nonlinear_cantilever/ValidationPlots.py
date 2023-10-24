"""
==============================================================================
Nonlinear cantilever under tip force - Validation plots
==============================================================================
@File    :   ValidationPlots.py
@Date    :   2023/06/08
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import pickle

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
# import niceplots
import matplotlib.pyplot as plt

# ==============================================================================
# Extension modules
# ==============================================================================
from analytic_displacement import analyticCantileverDisplacement
from TipForceReferenceResults import AlphaS4R, WtipS4R

LENGTH = 10.0

alphaTest = np.arange(0.05, 1.01, 0.05) * 4
alphaPlot = np.linspace(0.0, 1.0, 101) * 4

thetaTip, XDispNorm, ZDispNorm = analyticCantileverDisplacement(alphaPlot)

# ==============================================================================
# Z displacement plot
# ==============================================================================

fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")

ax.plot(alphaPlot, ZDispNorm, label="Analytic", clip_on=False, c="gray", lw=4)
ax.plot(
    AlphaS4R, WtipS4R / LENGTH, "o", label="Abaqus S4R", clip_on=False, markersize=10
)

for strainType in ["nonlinear"]:  # "linear"
    for rotationType in ["linear", "quadratic"]:
        try:
            with open(f"TACS-Disps-{strainType}_{rotationType}.pkl", "rb") as f:
                TACSDisps = pickle.load(f)
            ax.plot(
                TACSDisps["tipForce"] * 4,
                TACSDisps["zDisp"] / LENGTH,
                "o",
                label=f"TACS, {strainType} strain, {rotationType} rotation",
                clip_on=False,
                markersize=10,
            )
        except FileNotFoundError:
            pass

ax.set_xlabel(r"$\frac{F_{tip}L^2}{EI}$", fontsize=24)
ax.set_ylabel(
    r"$\frac{\Delta Z_{tip}}{L}$", rotation="horizontal", ha="right", fontsize=24
)
ax.set_xticks(np.arange(0, 5, 1))
ax.set_xticklabels([str(int(s)) for s in np.arange(0, 5, 1)])

ax.set_yticks(np.arange(0, 0.7, 0.1))
ax.set_yticklabels([f"{str(int(s*100))} %" for s in np.arange(0, 0.7, 0.1)])
ax.legend(labelcolor="linecolor")

plt.savefig("TipForceZDispValidation.png", bbox_inches="tight")
plt.show()
