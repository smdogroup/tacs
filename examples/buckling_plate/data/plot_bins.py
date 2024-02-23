import numpy as np
import matplotlib.pyplot as plt
import niceplots, pandas as pd, os

df = pd.read_csv("Nxcrit.csv")
Dstar = df["Dstar"].to_numpy()
affine_AR = df["a0/b0"].to_numpy()
slenderness = df["b/h"].to_numpy()
kx0 = df["kx_0"].to_numpy()

# a = df["a"].to_numpy()
# b = df["b"].to_numpy()
# AR = a / b

# make folder for the slender bins
slender_folder = os.path.join(os.getcwd(), "slender-bin")
if not os.path.exists(slender_folder):
    os.mkdir(slender_folder)

Dstar_folder = os.path.join(os.getcwd(), "Dstar-bin")
if not os.path.exists(Dstar_folder):
    os.mkdir(Dstar_folder)

# loop over different slenderness bins
slender_bins = [[5.0, 10.0], [10.0, 20.0], [20.0, 50.0], [50.0, 100.0], [100.0, 200.0]]

Dstar_bins = [[0.25 * i, 0.25 * (i + 1)] for i in range(6)]

plt.style.use(niceplots.get_style())

mask0 = kx0 <= 20.0

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
        mask = np.logical_and(mask, mask0)
        if np.sum(mask) == 0:
            continue

        plt.plot(
            affine_AR[mask],
            kx0[mask],
            "o",
            label=f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}]",
        )

    plt.legend()
    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$k_{x_0}$")
    plt.xlim(0.0, 5.0)
    plt.title(f"slender bin b/h in [{bin[0]:.1f},{bin[1]:.1f}]")
    plt.savefig(
        os.path.join(slender_folder, f"slender_bin_{bin[0]:.1f}_{bin[1]:.1f}.png"),
        dpi=400,
    )
    plt.close(f"islender{ibin}")


# now plot with Dstar_bins folder and those as the different plots
for iDstar, Dstar_bin in enumerate(Dstar_bins):
    if iDstar < len(Dstar_bins) - 1:
        mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar < Dstar_bin[1])
    else:
        mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar <= Dstar_bin[1])

    plt.figure(f"iDstar{iDstar}")
    plt.margins(x=0.05, y=0.05)

    if np.sum(mask2) == 0:
        continue

    for ibin, bin in enumerate(slender_bins):
        if ibin < len(slender_bins) - 1:
            mask1 = np.logical_and(bin[0] <= slenderness, slenderness < bin[1])
        else:
            mask1 = np.logical_and(bin[0] <= slenderness, slenderness <= bin[1])

        mask = np.logical_and(mask1, mask2)
        mask = np.logical_and(mask, mask0)
        if np.sum(mask) == 0:
            continue

        plt.plot(
            affine_AR[mask],
            kx0[mask],
            "o",
            label=f"b/h - [{bin[0]:.1f},{bin[1]:.1f}]",
        )

    plt.legend()
    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$k_{x_0}$")
    plt.xlim(0.0, 5.0)
    plt.title(f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}]")
    plt.savefig(
        os.path.join(
            Dstar_folder, f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}].png"
        ),
        dpi=400,
    )
    plt.close(f"iDstar{iDstar}")
