import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import niceplots, scipy, time, os
import argparse
from mpl_toolkits import mplot3d
from matplotlib import cm
import shutil

"""
This time I'll try a Gaussian Process model to fit the axial critical load surrogate model
Inputs: D*, a0/b0, ln(b/h)
Output: k_x0
"""
# parse the arguments
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--load', type=str)
parent_parser.add_argument('--BC', type=str)

args = parent_parser.parse_args()

assert args.load in ["Nx", "Nxy", "axial", "shear"]
assert args.BC in ["SS", "CL"]

print(f"args.load = {args.load}")
if args.load in ["Nx", "axial"]:
    load = "Nx"
else:
    load = "Nxy"
BC = args.BC

# load the Nxcrit dataset
load_prefix = "Nxcrit" if load == "Nx" else "Nxycrit"
csv_filename = f"{load_prefix}_{BC}"
print(f"csv filename = {csv_filename}")
df = pd.read_csv("data/" + csv_filename + ".csv")

# extract only the model columns
# TODO : if need more inputs => could maybe try adding log(E11/E22) in as a parameter?
# or also log(E11/G12)
X = df[["Dstar", "a0/b0", "b/h"]].to_numpy()
Y = df["kmin"].to_numpy()
Y = np.reshape(Y, newshape=(Y.shape[0], 1))

print(f"Monte Carlo #data = {X.shape[0]}")
N_data = X.shape[0]

n_train = int(0.9 * N_data)

# REMOVE THE OUTLIERS in local 4d regions
# loop over different slenderness bins
slender_bins = [
    [10.0, 20.0],
    [20.0, 50.0],
    [50.0, 100.0],
    [100.0, 200.0],
]  # [5.0, 10.0],
Dstar_bins = [[0.25 * i, 0.25 * (i + 1)] for i in range(7)]
# added smaller and larger bins here cause we miss some of the outliers near the higher a0/b0 with less data
aff_AR_bins = (
    [[0.5 * i, 0.5 * (i + 1)] for i in range(4)]
    + [[1.0 * i, 1.0 * (i + 1)] for i in range(2, 5)]
    + [[2.0, 10.0]]
)

# make a folder for the model fitting
data_folder = os.path.join(os.getcwd(), "data")
sub_data_folder = os.path.join(data_folder, csv_filename)
wo_outliers_folder = os.path.join(sub_data_folder, "model-no-outliers")
w_outliers_folder = os.path.join(sub_data_folder, "model-w-outliers")
GP_folder = os.path.join(sub_data_folder, "GP")
for ifolder,folder in enumerate([
    data_folder,
    sub_data_folder,
    wo_outliers_folder,
    w_outliers_folder,
    GP_folder,
]):
    if ifolder >= 2 and os.path.exists(folder):
        shutil.rmtree(folder)
    if not os.path.exists(folder):
        os.mkdir(folder)

plt.style.use(niceplots.get_style())
Dstar = X[:, 0]
affine_AR = X[:, 1]
slenderness = X[:, 2]
kx0 = Y[:, 0]

n_data = Dstar.shape[0]
# print(f"n data = {n_data}")
global_outlier_mask = np.full((n_data,), False, dtype=bool)

plt_ct = 0

for ibin, bin in enumerate(slender_bins):
    if ibin < len(slender_bins) - 1:
        mask1 = np.logical_and(bin[0] <= slenderness, slenderness < bin[1])
    else:
        mask1 = np.logical_and(bin[0] <= slenderness, slenderness <= bin[1])
    if np.sum(mask1) == 0:
        continue

    for iDstar, Dstar_bin in enumerate(Dstar_bins):
        if iDstar < len(Dstar_bins) - 1:
            mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar < Dstar_bin[1])
        else:
            mask2 = np.logical_and(Dstar_bin[0] <= Dstar, Dstar <= Dstar_bin[1])

        # limit to local regions in D*, slenderness
        layer2_mask = np.logical_and(mask1, mask2)

        for iAR, AR_bin in enumerate(aff_AR_bins):
            if iAR < len(aff_AR_bins) - 1:
                mask3 = np.logical_and(AR_bin[0] <= affine_AR, affine_AR < AR_bin[1])
            else:
                mask3 = np.logical_and(AR_bin[0] <= affine_AR, affine_AR <= AR_bin[1])

            mask = np.logical_and(layer2_mask, mask3)
            # and also now in affine_AR
            N = np.sum(mask)

            # print(f"aspect ratio bin = {AR_bin}, N = {N}")

            if N < 5:
                continue

            local_indices = np.where(mask)[0]

            # compute a local polynomial fit - mean and covariance
            # eliminate local outliers
            X_local = X[mask, :]
            t_local = Y[mask, :]

            loc_AR = X_local[:, 1:2]
            X_fit = np.concatenate([np.ones((N, 1)), loc_AR, loc_AR**2], axis=1)
            # print(f"X_fit shape = {X_fit.shape}")
            # print(f"t_local shape = {t_local.shape}")
            w_hat = np.linalg.solve(X_fit.T @ X_fit, X_fit.T @ t_local)
            # print(f"w_hat shape = {w_hat.shape}")

            # compute the local noise
            t_pred = X_fit @ w_hat
            resids = t_local - t_pred
            variance = float(1 / N * resids.T @ resids)
            sigma = np.sqrt(variance)

            _plot_outliers = False

            # compute plotted model
            if _plot_outliers:
                plot_AR = np.linspace(AR_bin[0], AR_bin[1], 20).reshape((20, 1))
                X_plot = np.concatenate(
                    [np.ones((20, 1)), plot_AR, plot_AR**2], axis=1
                )
                t_plot = X_plot @ w_hat

                # plot the local polynomial model and its variance range
                plt.figure("temp")
                ax = plt.subplot(111)
                plt.margins(x=0.05, y=0.05)
                plt.title(
                    f"b/h - [{bin[0]},{bin[1]}], D* - [{Dstar_bin[0]},{Dstar_bin[1]}]"
                )
                plt.plot(loc_AR, t_local, "ko")
                plt.plot(plot_AR, t_plot - 2 * sigma, "r--", label=r"$\mu-2 \sigma$")
                plt.plot(plot_AR, t_plot, "b-", label=r"$\mu$")
                plt.plot(plot_AR, t_plot + 2 * sigma, "b--", label=r"$\mu+2 \sigma$")
                plt.legend()
                plt.xlabel(r"$a_0/b_0$")
                plt.ylabel(r"$k_{x_0}$")
                plt.savefig(
                    os.path.join(
                        wo_outliers_folder, f"slender{ibin}_Dstar{iDstar}_AR{iAR}.png"
                    ),
                    dpi=400,
                )
                plt_ct += 1
                plt.close("temp")

            # get all of the datapoints that lie outside the +- 2 sigma bounds
            z_scores = abs(resids[:, 0]) / sigma
            # print(f"zscores = {z_scores}")
            outlier_mask = z_scores >= 2.0
            if np.sum(outlier_mask) == 0:
                continue  # no outliers found so exit
            X_outliers = X_local[outlier_mask, :]

            # now show plot with outliers indicated for sanity check
            if _plot_outliers:
                loc_AR_noout = loc_AR[np.logical_not(outlier_mask), :]
                t_local_noout = t_local[np.logical_not(outlier_mask), :]
                loc_AR_out = loc_AR[outlier_mask, :]
                t_local_out = t_local[outlier_mask, :]
                plt.figure("temp", figsize=(8, 6))
                ax = plt.subplot(111)
                plt.margins(x=0.05, y=0.05)
                plt.title(
                    f"b/h - [{bin[0]},{bin[1]}], D* - [{Dstar_bin[0]},{Dstar_bin[1]}]"
                )
                ax.fill_between(
                    plot_AR[:, 0],
                    t_plot[:, 0] - 2 * sigma,
                    t_plot[:, 0] + 2 * sigma,
                    color="g",
                )
                ax.plot(plot_AR, t_plot + 2 * sigma, "b--", label=r"$\mu+2 \sigma$")
                ax.plot(plot_AR, t_plot, "b-", label=r"$\mu$")
                ax.plot(plot_AR, t_plot - 2 * sigma, "r--", label=r"$\mu-2 \sigma$")
                ax.plot(loc_AR_noout, t_local_noout, "ko", label="data")
                ax.plot(loc_AR_out, t_local_out, "ro", label="outlier")
                plt.legend()
                plt.xlabel(r"$a_0/b_0$")
                plt.ylabel(r"$k_{x_0}$")
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
                plt.savefig(
                    os.path.join(
                        w_outliers_folder, f"slender{ibin}_Dstar{iDstar}_AR{iAR}.png"
                    ),
                    dpi=400,
                )
                plt_ct += 1
                plt.close("temp")

            # now record the outliers to the global outlier indices and remove them
            n_outliers = np.sum(outlier_mask)
            # print(f"local indices = {local_indices[:3]} type {type(local_indices)}")
            # print(f"outlier mask = {outlier_mask} type {type(outlier_mask)} shape = {outlier_mask.shape}")
            _outlier_indices = local_indices[outlier_mask]
            # print(f"local outlier indices = {_outlier_indices} shape = {_outlier_indices.shape}")
            for _outlier in _outlier_indices:
                global_outlier_mask[_outlier] = True

# print(f"global outlier mask = {global_outlier_mask}")
print(f"num outliers = {np.sum(global_outlier_mask)}")
# exit()
# remove the outliers from the dataset
keep_mask = np.logical_not(global_outlier_mask)
X = X[keep_mask, :]
Y = Y[keep_mask, :]

# convert a0/b0, b/h and kmin to log space
X[:, 1:] = np.log(X[:, 1:])
Y[:, 0:] = np.log(Y[:, 0:])

# split into training and test datasets
n_total = X.shape[0]
n_test = n_total - n_train
assert n_test > 100

X_train = X[:n_train, :]
X_test = X[n_train:, :]
Y_train = Y[:n_train, :]
Y_test = Y[n_train:, :]

# print(f"Y train = {Y_train.shape}")

# TRAIN THE MODEL HYPERPARAMETERS
# by maximizing the log marginal likelihood log p(y|X)
# or minimizing the negative log marginal likelihood
# ----------------------------------------------------

# y is the training set observations
y = Y_train

# update the local hyperparameter variables
# initial hyperparameter vector
# sigma_n, sigma_f, L1, L2, L3
theta0 = np.array([1e-1, 3e-1, 1.0,
                  0.2, 1.0, 1.0, 0.3])
sigma_n = 1e-3

def relu(x):
    return max([0.0, x])

def soft_relu(x,rho=10):
    return 1.0/rho * np.log(1 + np.exp(rho*x))

def kernel(xp, xq, theta):
    # xp, xq are Nx1,Mx1 vectors (D*, a0/b0, ln(b/h))
    vec = xp - xq

    S1 = theta[0]; S2 = theta[1]; c = theta[2]; L1 = theta[3]; S4 = theta[4]; S5 = theta[5]; L2 = theta[6]

    d1 = vec[1] # first two entries
    d2 = vec[2]

    # kernel in 
    kernel0 = S1**2 + S2**2 * (xp[0]-c) * (xq[0]-c)
    kernel1 = np.exp(-0.5 * (d1**2/L1**2)) * soft_relu(1-abs(xp[1])) * soft_relu(1-abs(xq[1])) + \
            S4 + S5 * soft_relu(-xp[1]) * soft_relu(-xq[1])
    kernel2 = np.exp(-0.5 * d2**2 / L2**2)
    return kernel0 * kernel1 * kernel2


# compute the training kernel matrix
K_y = np.array(
    [
        [kernel(X_train[i, :], X_train[j, :],theta0) for i in range(n_train)]
        for j in range(n_train)
    ]
) + sigma_n**2 * np.eye(n_train)

# compute the objective : maximize log marginal likelihood
# sign,log_detK = np.linalg.slogdet(K_y) # special numpy routine for log(det(K))
# print(f"\tlog detK = {log_detK}, sign = {sign}")
# _start = time.time()
alpha = np.linalg.solve(K_y, y)

# plot the model and some of the data near the model range in D*=1, AR from 0.5 to 5.0, b/h=100
# ---------------------------------------------------------------------------------------------
_plot = True
_plot_Dstar_2d = True
_plot_slender_2d = True
_plot_3d = True

if _plot:
    # get the available colors
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    plt.style.use(niceplots.get_style())

    if _plot_3d:
        # plot the 3D version of the GP curve for each slenderness value
        # iterate over the different slender,D* bins
        for ibin, bin in enumerate(slender_bins):
            slender_bin = [np.log(bin[0]), np.log(bin[1])]
            avg_log_slender = 0.5 * (slender_bin[0] + slender_bin[1])
            mask1 = np.logical_and(slender_bin[0] <= X[:, 2], X[:, 2] <= slender_bin[1])
            if np.sum(mask1) == 0:
                continue

            fig = plt.figure("3d-GP", figsize =(14, 9))
            ax = plt.axes(projection ='3d',computed_zorder=False)

            # plot data in certain range of the training set
            for iDstar, Dstar_bin in enumerate(Dstar_bins):
                mask2 = np.logical_and(Dstar_bin[0] <= X[:, 0], X[:, 0] <= Dstar_bin[1])
                avg_Dstar = 0.5 * (Dstar_bin[0] + Dstar_bin[1])

                mask = np.logical_and(mask1, mask2)
                if np.sum(mask) == 0:
                    continue
                X_in_range = X[mask, :]
                Y_in_range = Y[mask, :]
                ax.scatter(X_in_range[:,0], X_in_range[:,1], Y_in_range[:,0], s=40, edgecolors='black', zorder=1+iDstar)

            n_plot = 3000
            X_plot_mesh = np.zeros((30,100))
            X_plot = np.zeros((n_plot, 3))
            ct = 0
            Dstar_vec = np.linspace(0.25, 1.75, 30)
            AR_vec = np.log(np.linspace(0.1, 10.0, 100))
            for iDstar in range(30):
                for iAR in range(100):
                    X_plot[ct,:] = np.array([Dstar_vec[iDstar], AR_vec[iAR], avg_log_slender])
                    ct += 1

            Kplot = np.array(
                [
                    [kernel(X_train[i, :], X_plot[j, :],theta0) for i in range(n_train)]
                    for j in range(n_plot)
                ]
            )
            f_plot = Kplot @ alpha

            # make meshgrid of outputs
            DSTAR = np.zeros((30, 100))
            AR = np.zeros((30,100))
            KMIN = np.zeros((30,100))
            ct = 0
            for iDstar in range(30):
                for iAR in range(100):
                    DSTAR[iDstar,iAR] = Dstar_vec[iDstar]
                    AR[iDstar,iAR] = AR_vec[iAR]
                    KMIN[iDstar,iAR] = f_plot[ct]
                    ct += 1

            # plot the model curve            
            # Creating plot
            face_colors = cm.jet(KMIN/10.0)
            ax.plot_surface(DSTAR, AR, KMIN, antialiased=False, facecolors = face_colors, alpha=0.4, zorder=1)

            # save the figure
            ax.set_xlabel(r"$D^*$")
            ax.set_ylabel(r"$log(a_0/b_0)$")
            ax.set_zlabel(r"$log(k_{min})$")
            ax.set_ylim3d(np.log(0.1), np.log(10.0))
            ax.set_zlim3d(0.0, np.log(30.0))
            ax.view_init(elev=20, azim=20, roll=0)
            plt.gca().invert_xaxis()
            plt.title(f"b/h in [{bin[0]},{bin[1]}]")
            #plt.show()
            plt.savefig(os.path.join(GP_folder, f"3d-slender-{ibin}.png"), dpi=400)
            plt.close("3d-GP")

    if _plot_slender_2d:
        # iterate over the different slender,D* bins
        for ibin, bin in enumerate(slender_bins):
            slender_bin = [np.log(bin[0]), np.log(bin[1])]
            avg_log_slender = 0.5 * (slender_bin[0] + slender_bin[1])
            mask1 = np.logical_and(slender_bin[0] <= X[:, 2], X[:, 2] <= slender_bin[1])
            if np.sum(mask1) == 0:
                continue

            plt.figure("check model", figsize=(8, 6))
            plt.margins(x=0.05, y=0.05)
            plt.title(f"b/h in [{bin[0]},{bin[1]}]")
            ax = plt.subplot(111)

            for iDstar, Dstar_bin in enumerate(Dstar_bins):
                mask2 = np.logical_and(Dstar_bin[0] <= X[:, 0], X[:, 0] <= Dstar_bin[1])
                avg_Dstar = 0.5 * (Dstar_bin[0] + Dstar_bin[1])

                mask = np.logical_and(mask1, mask2)
                if np.sum(mask) == 0:
                    continue

                # predict the GP curve
                n_plot = 300
                X_plot = np.zeros((n_plot, 3))
                X_plot[:, 0] = avg_Dstar
                X_plot[:, 1] = np.linspace(0.1, 10.0, n_plot)
                X_plot[:, 2] = avg_log_slender

                Kplot = np.array(
                    [
                        [kernel(X_train[i, :], X_plot[j, :],theta0) for i in range(n_train)]
                        for j in range(n_plot)
                    ]
                )
                f_plot = Kplot @ alpha

                # plot data in certain range of the training set
                mask = np.logical_and(mask1, mask2)
                X_in_range = X[mask, :]
                Y_in_range = Y[mask, :]

                AR_in_range = X_in_range[:, 1]

                # plot the raw data and the model in this range
                ax.plot(
                    AR_in_range, Y_in_range, "o", color=colors[iDstar]
                )  # , label=f"D*-[{Dstar_bin[0]},{Dstar_bin[1]}]""
                ax.plot(
                    X_plot[:, 1],
                    f_plot[:, 0],
                    color=colors[iDstar],
                    label=f"D*-[{Dstar_bin[0]},{Dstar_bin[1]}]",
                )

            # outside of for loop save the plot
            plt.xlabel(r"$log(a_0/b_0)$")
            plt.ylabel(r"$log(k_{x_0})$")
            plt.legend()
            plt.xlim(np.log(0.1), np.log(20.0))
            plt.ylim(0.0, np.log(20.0))
            # box = ax.get_position()
            # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # # Put a legend to the right of the current axis
            # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(os.path.join(GP_folder, f"slender{ibin}.png"), dpi=400)
            plt.close("check model")

    if _plot_Dstar_2d:
        # then do the same but reverse D* then slender for loop
        # iterate over the different slender,D* bins
        for iDstar, Dstar_bin in enumerate(Dstar_bins):
            mask2 = np.logical_and(Dstar_bin[0] <= X[:, 0], X[:, 0] <= Dstar_bin[1])
            avg_Dstar = 0.5 * (Dstar_bin[0] + Dstar_bin[1])
            if np.sum(mask2) == 0:
                continue

            plt.figure("check model", figsize=(8, 6))
            plt.margins(x=0.05, y=0.05)
            plt.title(f"D* in [{Dstar_bin[0]},{Dstar_bin[1]}]")
            ax = plt.subplot(111)

            for ibin, bin in enumerate(slender_bins):
                slender_bin = [np.log(bin[0]), np.log(bin[1])]
                avg_log_slender = 0.5 * (slender_bin[0] + slender_bin[1])
                mask1 = np.logical_and(slender_bin[0] <= X[:, 2], X[:, 2] <= slender_bin[1])

                mask = np.logical_and(mask1, mask2)
                if np.sum(mask) == 0:
                    continue

                # predict the GP curve
                n_plot = 300
                X_plot = np.zeros((n_plot, 3))
                X_plot[:, 0] = avg_Dstar
                X_plot[:, 1] = np.linspace(0.1, 10.0, n_plot)
                X_plot[:, 2] = avg_log_slender

                Kplot = np.array(
                    [
                        [kernel(X_train[i, :], X_plot[j, :],theta0) for i in range(n_train)]
                        for j in range(n_plot)
                    ]
                )
                f_plot = Kplot @ alpha

                # plot data in certain range of the training set
                mask = np.logical_and(mask1, mask2)
                X_in_range = X[mask, :]
                Y_in_range = Y[mask, :]

                AR_in_range = X_in_range[:, 1]

                # plot the raw data and the model in this range
                ax.plot(
                    AR_in_range, Y_in_range, "o", color=colors[ibin]
                )  # , label=f"D*-[{Dstar_bin[0]},{Dstar_bin[1]}]""
                ax.plot(
                    X_plot[:, 1],
                    f_plot[:, 0],
                    color=colors[ibin],
                    label=f"b/h-[{bin[0]},{bin[1]}]",
                )

            # outside of for loop save the plot
            plt.xlabel(r"$log(a_0/b_0)$")
            plt.ylabel(r"$log(k_{x_0})$")
            plt.legend()
            plt.xlim(np.log(0.1), np.log(20.0))
            plt.ylim(0.0, np.log(20.0))
            # box = ax.get_position()
            # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # # Put a legend to the right of the current axis
            # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(os.path.join(GP_folder, f"Dstar{iDstar}.png"), dpi=400)
            plt.close("check model")

# predict against the test dataset
# --------------------------------------------------------------------------------------------------
# K(X*,X)
Kstar = np.array(
    [
        [kernel(X_train[i, :], X_test[j, :],theta0) for i in range(n_train)]
        for j in range(n_test)
    ]
)

f_test = Kstar @ alpha
# compute the RMSE on the test dataset
test_resid = Y_test - f_test
# print(f"test resid = {test_resid}")
RMSE = np.sqrt(1.0 / n_test * float(test_resid.T @ test_resid))

# norm RMSE
norm_resid = (Y_test - f_test) / Y_test
avg_rel_error = np.sum(np.abs(norm_resid)) / n_test

# write this to file in the GP folder
txt_file = os.path.join(GP_folder, "model-fit.txt")
hdl = open(txt_file, "w")
hdl.write(f"RMSE test = {RMSE}")
hdl.write(f"avg relative error = {avg_rel_error}")
hdl.close()
