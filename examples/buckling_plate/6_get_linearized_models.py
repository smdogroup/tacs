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

# ignore all AR data below 2.0
mask = np.logical_and(X[:,1] > 2.0, X[:,1] < 5.0)
X = X[mask,:]
Y = Y[mask,:]

# remove outliers
if load == "Nxy":
    mask = Y[:,0] > 2.0
    X = X[mask,:]
    Y = Y[mask,:]
elif load == "Nx":
    mask = Y[:,0] > 2.2
    X = X[mask,:]
    Y = Y[mask,:]

# convert to sigmoid(b/h) data
X[:,2] = 1.0 / (1.0 + np.exp(-1.0 * X[:,2] / 10) ) - 1.0

print(f"X = {X}")

n_data = X.shape[0]
n_train = int(0.8 * n_data)
n_test = int(0.2 * n_data)

include_quadratic = load == "Nxy"
include_quadratic = False

if include_quadratic:

    # linear model kmin = w_0 + Dstar * w_1 + Dstar**2 * w_2 + log(b/h) * w_3
    X_linear = np.concatenate((np.ones((n_train,1)), X[:n_train,0:1], X[:n_train,0:1]**2, X[:n_train,2:3]), axis=-1)
    Y_linear = Y[:n_train,:]

    w_linear = np.linalg.solve(X_linear.T @ X_linear, X_linear.T @ Y_linear)

    print(f"{load}_{BC}_crit norm = {w_linear[0,0]} + {w_linear[1,0]} * Dstar + {w_linear[2,0]} * Dstar**2 + {w_linear[3,0]} * (sigmoid(b/h/10)-1)")

else:


    # linear model kmin = w_0 + Dstar * w_1 + log(b/h) * w_2
    X_linear = np.concatenate((np.ones((n_train,1)), X[:n_train,0:1], X[:n_train,2:3]), axis=-1)
    Y_linear = Y[:n_train,:]

    print(f"X linear = {X_linear}")

    w_linear = np.linalg.solve(X_linear.T @ X_linear, X_linear.T @ Y_linear)

    print(f"{load}_{BC}_crit norm = {w_linear[0,0]} + {w_linear[1,0]} * Dstar + {w_linear[2,0]} * (sigmoid(b/h/10)-1)")

# plot the b/h data versus kmin
#mask = np.logical_and(X[:,0] < 1.1, X[:,0] > 0.9)
plt.plot(X[:,0:1], Y[:,:], "ko")
plt.savefig(load + "_Dstar.png", dpi=400)


# then also get the avg relative error of the model
X_test = X[n_train:,:]
Y_test = Y[n_train:,:]

if include_quadratic:
    X_linear_test = np.concatenate((np.ones((X_test.shape[0],1)), X_test[:,0:1], X_test[:,0:1]**2, X_test[:,2:3]), axis=-1)
else:
    X_linear_test = np.concatenate((np.ones((X_test.shape[0],1)), X_test[:,0:1], X_test[:,2:3]), axis=-1)
Y_pred = X_linear_test @ w_linear

n_pred = Y_pred.shape[0]
rel_err_sum = 0.0
for i in range(n_pred):
    rel_err = (Y_pred[i,0] - Y_test[i,0]) / Y_test[i,0]
    rel_err_sum += abs(rel_err)
avg_rel_err = rel_err_sum / n_pred

print(f"avg rel err model = {avg_rel_err}")

# now compare tacs closed-form
w_tacs = w_linear * 0.0
if load == "Nxy":
    Y_pred = np.zeros((X_test.shape[0],1))
    for i in range(X_test.shape[0]):
        xi = 1.0/X_test[i,0]
        #xi = X_test[i,0]
        if xi > 1.0:
            Y_pred[i,0] = 4.0/np.pi**2 * (8.125 + 5.045 / xi)
        else:
            Y_pred[i,0] = 4.0/np.pi**2 * xi**0.5 * (11.7 + 0.532 * xi + 0.938 * xi**2)
else:
    Y_pred = 2.0 * (1.0 + X_test[:,0])

Dstar = np.linspace(0.0, 1.5, 100)
k = np.zeros((100,))
for i in range(100):
    xi = 1.0/Dstar[i]
    #xi = Dstar[i]
    if xi > 1.0:
        k[i] = 4.0/np.pi**2 * (8.125 + 5.045 / xi)
    else:
        k[i] = 4.0/np.pi**2 / xi**0.5 * (11.7 + 0.532 * xi + 0.938 * xi**2)
plt.plot(Dstar, k, linewidth=2)
plt.ylim(0.0, 6.5)
plt.savefig("Dstar-plot.png", dpi=400)

n_pred = Y_pred.shape[0]
rel_err_sum = 0.0
for i in range(n_pred):
    rel_err = (Y_pred[i,0] - Y_test[i,0]) / Y_test[i,0]
    rel_err_sum += abs(rel_err)
avg_rel_err = rel_err_sum / n_pred

print(f"avg rel err tacs = {avg_rel_err}")