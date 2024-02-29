import pandas as pd
import matplotlib.pyplot as plt
import niceplots
import numpy as np

df = pd.read_csv("compare-cases.csv")
AR_list = df["AR"].to_numpy()
kx0_FEA_list = df["kx0_FEA"].to_numpy()
BC = df["BC"].to_numpy()
load = df["load"].to_numpy()

# plot the Closed-Form versus the FEA for affine equations
plt.style.use(niceplots.get_style())
plt.margins(x=0.02, y=0.02)
plt.figure("affine")
for _BC in ["CL", "SS"]:
    for _load in ["axial", "shear"]:
        BC_mask = BC == _BC
        load_mask = load == _load
        mask = np.logical_and(BC_mask, load_mask)
        plt.plot(AR_list[mask], kx0_FEA_list[mask], "o", label=_load + "-" + _BC)
plt.xlabel(r"$a_0/b_0$")
plt.ylabel(r"$k$")
plt.ylim(0.0, 20.0)
plt.xlim(0.0, 5.0)
plt.legend()
plt.savefig("compare-cases.png", dpi=400)
