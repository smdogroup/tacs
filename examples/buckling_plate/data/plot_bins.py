import numpy as np
import matplotlib.pyplot as plt
import niceplots, pandas as pd, os
import shutil
from mpi4py import MPI
import time
import argparse

# parse the arguments
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--load', type=str)
parent_parser.add_argument('--BC', type=str)
parent_parser.add_argument('--log', type=bool, default=False)

args = parent_parser.parse_args()

assert args.load in ["Nx", "Nxy", "axial", "shear"]
assert args.BC in ["SS", "CL"]

if args.load in ["Nx", "axial"]:
    loading = "Nx"
else:
    loading = "Nxy"
BC = args.BC
use_logs = args.log
# print(f"use logs = {use_logs}", flush=True)
# exit()

csv_filename = loading + "crit_" + BC + ".csv"

comm = MPI.COMM_WORLD

df = pd.read_csv(csv_filename)
Dstar = df["Dstar"].to_numpy()
affine_AR = df["a0/b0"].to_numpy()
slenderness = df["b/h"].to_numpy()
kx0 = df["kmin"].to_numpy()
kmodes = [df[f"k_{imode+1}"].to_numpy() for imode in range(20)]

# a = df["a"].to_numpy()
# b = df["b"].to_numpy()
# AR = a / b

case_folder = os.path.join(os.getcwd(), csv_filename.split(".")[0])
if not os.path.exists(case_folder):
    os.mkdir(case_folder)

# make folder for the slender bins
slender_folder = os.path.join(case_folder, "slender-bin")
if os.path.exists(slender_folder):
    shutil.rmtree(slender_folder)
if not os.path.exists(slender_folder):
    os.mkdir(slender_folder)

Dstar_folder = os.path.join(case_folder, "Dstar-bin")
if os.path.exists(Dstar_folder):
    shutil.rmtree(Dstar_folder)
if not os.path.exists(Dstar_folder):
    os.mkdir(Dstar_folder)

comm.Barrier()
time.sleep(3)

# loop over different slenderness bins
slender_bins = [[5.0, 10.0], [10.0, 20.0], [20.0, 50.0], [50.0, 100.0], [100.0, 200.0]]

Dstar_bins = [[0.25 * i, 0.25 * (i + 1)] for i in range(7)]

plt.style.use(niceplots.get_style())

mask0 = kx0 <= 20.0

figsize=(8,6)

for ibin, bin in enumerate(slender_bins):
    if ibin < len(slender_bins) - 1:
        mask1 = np.logical_and(bin[0] <= slenderness, slenderness < bin[1])
    else:
        mask1 = np.logical_and(bin[0] <= slenderness, slenderness <= bin[1])
    if np.sum(mask1) == 0:
        continue

    # slender bin folders
    c_slender_folder = os.path.join(slender_folder, f"b_h_{bin[0]:.1f}_{bin[1]:.1f}")
    if not os.path.exists(c_slender_folder):
        os.mkdir(c_slender_folder)

    # plot the kmin values
    # -------------------------------------------------------------------------
    plt.figure(f"islender{ibin}", figsize=figsize)
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
    if use_logs: 
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.xlim(0.0, 5.0)
        plt.ylim(0.0, 20.0)
    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$\lambda_{min}$")
    plt.title(f"slender bin b/h in [{bin[0]:.1f},{bin[1]:.1f}]")
    plt.savefig(
        os.path.join(c_slender_folder, "_kmin.png"),
        dpi=400,
    )
    plt.close(f"islender{ibin}")

    # plot the individual \lambda_i modes
    # -------------------------------------------------------------------------
    for imode in range(20):
        has_values = False

        for iDstar, Dstar_bin in enumerate(Dstar_bins):
            if iDstar < len(Dstar_bins) - 1:
                mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar < Dstar_bin[1])
            else:
                mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar <= Dstar_bin[1])

            mask = np.logical_and(mask1, mask2)
            mask = np.logical_and(mask, mask0)
            mask = np.logical_and(mask, kmodes[imode] <= 20.0)

            if np.sum(mask) == 0:
                continue

            if any(kmodes[imode][mask]) is not None and any(kmodes[imode][mask]) < 20.0:
                has_values = True

        if not has_values:
            continue

        plt.figure(f"islender{ibin}-mode{imode}", figsize=figsize)
        plt.margins(x=0.05, y=0.05)

        for iDstar, Dstar_bin in enumerate(Dstar_bins):
            if iDstar < len(Dstar_bins) - 1:
                mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar < Dstar_bin[1])
            else:
                mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar <= Dstar_bin[1])

            mask = np.logical_and(mask1, mask2)
            mask = np.logical_and(mask, mask0)
            mask = np.logical_and(mask, kmodes[imode] <= 20.0)
            if np.sum(mask) == 0:
                continue

            plt.plot(
                affine_AR[mask],
                kmodes[imode][mask],
                "o",
                label=f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}]",
            )

        plt.legend()
        plt.xlabel(r"$a_0/b_0$")
        plt.ylabel(r"$\lambda_{" + str(imode + 1) + r"}$")
        if use_logs: 
            plt.xscale('log')
            plt.yscale('log')
        else:
            plt.xlim(0.0, 5.0)
            plt.ylim(0.0, 20.0)
        plt.title(f"slender bin b/h in [{bin[0]:.1f},{bin[1]:.1f}]")
        plt.savefig(
            os.path.join(c_slender_folder, f"kmode_{imode+1}.png"),
            dpi=400,
        )
        plt.close(f"islender{ibin}-mode{imode}")

    # show each of the modes (w/o D* slicing)
    plt.figure("all-modes", figsize=(8, 6))
    ax = plt.subplot(111)
    for imode in range(20):
        mask = np.logical_and(mask1, mask0)
        if any(kmodes[imode][mask]) is not None:
            pass
        else:
            continue

        ax.plot(
            affine_AR[mask],
            kmodes[imode][mask],
            "o",
            label=r"$\lambda_{" + str(imode + 1) + r"}$",
        )

    # Shrink current axis by 20%
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # maybe in log-log scale the modes will become nearly bilinear?
    if use_logs: 
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.xlim(0.0, 5.0)
        plt.ylim(0.0, 20.0)

    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$\lambda$")
    plt.title(f"slender bin b/h in [{bin[0]:.1f},{bin[1]:.1f}]")
    plt.savefig(
        os.path.join(c_slender_folder, f"_kmodes.png"),
        dpi=400,
    )
    plt.close("all-modes")


# now plot with Dstar_bins folder and those as the different plots
for iDstar, Dstar_bin in enumerate(Dstar_bins):
    if iDstar < len(Dstar_bins) - 1:
        mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar < Dstar_bin[1])
    else:
        mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar <= Dstar_bin[1])

    if np.sum(mask2) == 0:
        continue

    # Dstar folders
    c_dstar_folder = os.path.join(
        Dstar_folder, f"D*-[{Dstar_bin[0]:.2f}-{Dstar_bin[1]:.2f}]"
    )
    if not os.path.exists(c_dstar_folder):
        os.mkdir(c_dstar_folder)

    # plot kmin
    # ------------------------------------------------------
    plt.figure(f"iDstar{iDstar}", figsize=figsize)
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
    if use_logs: 
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.xlim(0.0, 5.0)
        plt.ylim(0.0, 20.0)
    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$\lambda_{min}$")
    plt.title(f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}]")
    plt.savefig(
        os.path.join(c_dstar_folder, "_kmin.png"),
        dpi=400,
    )
    plt.close(f"iDstar{iDstar}")

    # plot the individual \lambda_i modes
    # -------------------------------------------------------------------------
    for imode in range(20):
        has_values = False

        for ibin, bin in enumerate(slender_bins):
            if ibin < len(slender_bins) - 1:
                mask1 = np.logical_and(bin[0] <= slenderness, slenderness < bin[1])
            else:
                mask1 = np.logical_and(bin[0] <= slenderness, slenderness <= bin[1])

            mask = np.logical_and(mask1, mask2)
            mask = np.logical_and(mask, mask0)
            mask = np.logical_and(mask, kmodes[imode] <= 20.0)
            if np.sum(mask) == 0:
                continue

            if any(kmodes[imode][mask]) is not None and any(kmodes[imode][mask]) < 20.0:
                has_values = True

        if not has_values:
            continue

        plt.figure(f"iDstar{iDstar}-imode{imode+1}", figsize=figsize)
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
                kmodes[imode][mask],
                "o",
                label=f"b/h - [{bin[0]:.1f},{bin[1]:.1f}]",
            )

        plt.legend()
        if use_logs: 
            plt.xscale('log')
            plt.yscale('log')
        else:
            plt.xlim(0.0, 5.0)
            plt.ylim(0.0, 20.0)
        plt.xlabel(r"$a_0/b_0$")
        plt.ylabel(r"$\lambda_{min}$")
        plt.title(f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}]")
        plt.savefig(
            os.path.join(c_dstar_folder, f"kmode_{imode+1}.png"),
            dpi=400,
        )
        plt.close(f"iDstar{iDstar}-imode{imode+1}")

    # show each of the modes (w/o D* slicing)
    plt.figure("all-modes", figsize=(8, 6))
    ax = plt.subplot(111)
    for imode in range(20):
        mask = np.logical_and(mask2, mask0)
        if any(kmodes[imode][mask]) is not None:
            pass
        else:
            continue

        ax.plot(
            affine_AR[mask],
            kmodes[imode][mask],
            "o",
            label=r"$\lambda_{" + str(imode + 1) + r"}$",
        )

    # Shrink current axis by 20%
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # maybe in log-log scale the modes will become nearly bilinear?
    if use_logs: 
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.xlim(0.0, 5.0)
        plt.ylim(0.0, 20.0)

    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$\lambda$")
    plt.title(f"D* in [{Dstar_bin[0]:.2f},{Dstar_bin[1]:.2f}]")
    plt.savefig(
        os.path.join(c_dstar_folder, f"_kmodes.png"),
        dpi=400,
    )
    plt.close("all-modes")
