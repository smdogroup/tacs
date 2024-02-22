import numpy as np
import matplotlib.pyplot as plt
import niceplots, pandas as pd, os

df = pd.read_csv("Nxcrit.csv")
Dstar = df["Dstar"].to_numpy()
affine_AR = df["a0/b0"].to_numpy()
slenderness = df["b/h"].to_numpy()
kx0 = df["kx_0"].to_numpy()

# make folder for the slender bins
slender_folder = os.path.join(os.getcwd(), "slender-bin")
if not os.path.exists(slender_folder):
    os.mkdir(slender_folder)

# loop over different slenderness bins
slender_bins = [[5.0, 10.0], [10.0, 20.0], [20.0, 50.0], [50.0, 100.0]]

Dstar_bins = [[0.0, 0.2], [0.2, 0.4], [0.4, 0.6], [0.6, 0.8], [0.8, 1.0]]

plt.style.use(niceplots.get_style())

for ibin, bin in enumerate(slender_bins):
    if ibin < len(slender_bins) - 1:
        mask1 = np.logical_and(bin[0] <= slenderness, slenderness < bin[1])
    else:
        mask1 = np.logical_and(bin[0] <= slenderness, slenderness <= bin[1])
    if np.sum(mask1) == 0:
        continue

    plt.figure(f"islender{ibin}")
    plt.margins(x=0.05, y=0.05)

    for iDstar, Dstar_bin in enumerate(Dstar_bins):
        if iDstar < len(Dstar_bins) - 1:
            mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar < Dstar_bin[1])
        else:
            mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar <= Dstar_bin[1])

        mask = np.logical_and(mask1, mask2)
        if np.sum(mask) == 0:
            continue

        plt.plot(
            affine_AR[mask],
            kx0[mask],
            "o",
            label=f"D* in [{Dstar_bin[0]:.1f},{Dstar_bin[1]:.1f}]",
        )

    plt.legend()
    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$k_{x_0}$")
    plt.title(f"slender bin b/h in [{bin[0]:.1f},{bin[1]:.1f}]")
    plt.savefig(
        os.path.join(slender_folder, f"slender_bin_{bin[0]:.1f}_{bin[1]:.1f}.png"),
        dpi=400,
    )
    plt.close(f"islender{ibin}")
