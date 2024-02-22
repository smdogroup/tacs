import pandas as pd
import matplotlib.pyplot as plt
import niceplots

df = pd.read_csv("axial-SS1.csv")
AR_list = df["AR"].to_numpy()
kx0_CF_list = df["kx0_CF"].to_numpy()
kx0_FEA_list = df["kx0_FEA"].to_numpy()

# plot the Closed-Form versus the FEA for affine equations
plt.style.use(niceplots.get_style())
plt.margins(x=0.02, y=0.02)
plt.figure("affine")
plt.vlines(x=0.5, ymin=0.0, ymax=7.5, linestyles="--", colors="k", linewidth=1)
plt.plot(AR_list, kx0_CF_list, label="closed-form")
plt.plot(AR_list, kx0_FEA_list, "o", label="FEA")
plt.xlabel(r"$a_0/b_0$")
plt.ylabel(r"$k_{x_0}$")
plt.ylim(0.0, 7.5)
plt.xlim(0.0, 5.2)
plt.legend()
plt.savefig("affine-AR-kx0.png", dpi=400)
