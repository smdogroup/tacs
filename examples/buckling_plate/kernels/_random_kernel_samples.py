import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import niceplots, scipy, time, os
import argparse

parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--case', type=int)

args = parent_parser.parse_args()

assert args.case in [1, 2, 3, 4, 5, 6, 7, 8]

def relu(x):
    return max([0, x])
if args.case == 1: # large length scale SE kernel
    def kernel(x,xp):
        return np.exp(-0.5 * (x-xp)**2)
elif args.case == 2: # short length scale SE kernel
    def kernel(x,xp):
        return np.exp(-0.5 * (x-xp)**2/0.2**2)
elif args.case == 3: # linear kernel (centered at c=0)
    def kernel(x,xp):
        return 1 + 0.2 * x * xp
elif args.case == 4:
    def kernel(x,xp):
        return 1 + 0.2 * max([x * xp, 0])
elif args.case == 5:
    def kernel(x,xp):
        return np.exp(-0.5 * (x-xp)**2/0.2**2) * relu(1-abs(x)) * relu(1 - abs(xp))  + 1 + 1.0 * relu(-x) * relu(-xp)
elif args.case == 6:
    def kernel(x,xp):
        return 1.0
elif args.case == 7:
    def kernel(x,xp):
        return 1 + 0.2 * relu(-x) * relu(-xp)
elif args.case == 8:
    def soft_relu(x,rho=3):
        return 1.0/rho * np.log(1 + np.exp(rho*x))
    def kernel(x,xp):
        return np.exp(-0.5 * (x-xp)**2/0.2**2) * soft_relu(1-abs(x)) * soft_relu(1 - abs(xp))  + 1 + 1.0 * soft_relu(-x) * soft_relu(-xp)



# generate random samples of the kernels
n = 300
xvec = np.linspace(-3.0, 3.0, n)
mean = np.zeros((n,))
covar = np.array([[kernel(xvec[i], xvec[j]) for i in range(n)] for j in range(n)])

nsamples = 30
samples = np.random.multivariate_normal(mean, covar, size=nsamples)

colors = plt.cm.jet(np.linspace(0,1,nsamples))

fig, ax = plt.subplots()
plt.style.use(niceplots.get_style())

ymin = np.min(np.array(samples))
ymax = np.max(np.array(samples))
print(f"ymin = {ymin}")
print(f"ymax = {ymax}")
if args.case in [5, 8]:
    for xval in [-1,1]:
        ax.vlines(xval, ymin=ymin, ymax=ymax, linestyles="--", colors="k")

for isample,sample in enumerate(samples):
    ax.plot(xvec, sample,color=colors[isample], linewidth=1)
#ax.set_axis_off(

plt.savefig(f"kernel-demo{args.case}.png", dpi=400)