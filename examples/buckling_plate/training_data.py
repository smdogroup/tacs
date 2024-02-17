from ._generate_plate import generate_plate
from ._static_analysis import run_static_analysis
from ._buckling_analysis import run_buckling_analysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import niceplots

data_dict = {
    "Lx": [],
    "Ly": [],
    "AR": [],
    "thick": [],
    "E": [],
    "nu": [],
    "D": [],
    "Sx": [],
    "Sxy": [],
    "Sy": [],
    "eig1": [],
    "eig2": [],
    "eig3": [],
    "eig4": [],
    "eig5": [],
}

Lx = 1.0
thick = 0.01
E = 70e9
nu = 0.33
D = E * thick**3 / 12.0 / (1 - nu**2)
# for log_Ly in np.linspace(-1.0, 1.0, 10):
#    Ly = 10**log_Ly
for Ly in np.linspace(1.0, 10.0, 1000):
    # collect input data
    data_dict["Lx"] += [Lx]
    data_dict["Ly"] += [Ly]
    data_dict["AR"] += [Ly / Lx]
    data_dict["thick"] += [thick]
    data_dict["E"] += [E]
    data_dict["nu"] += [nu]
    data_dict["D"] += [D]

    generate_plate(Lx=Lx, Ly=Ly, nx=12, ny=12, displacement_control=False)
    avgStresses = run_static_analysis(
        thickness=thick, E=E, nu=nu, displacement_control=False
    )
    funcs = run_buckling_analysis(
        thickness=thick, E=E, nu=nu, sigma=30.0, num_eig=20, displacement_control=False
    )
    print(f"training data funcs = {funcs}")

    # save the average mid-plane stresses in the plate
    for i, key in enumerate(["Sx", "Sy", "Sxy"]):
        data_dict[key] = avgStresses[i]

    # save the first 5 eigenvalues
    for i in range(5):
        data_dict[f"eig{i+1}"] += [funcs[i]]

print(data_dict["AR"])
print(data_dict["eig1"])

# plot some of the training data
plt.style.use(niceplots.get_style())
plt.margins(x=0.04, y=0.04)
plt.plot(data_dict["AR"], data_dict["eig1"])
plt.xlabel("AR")
plt.ylabel(r"$\lambda_1$")
# plt.xscale('log')
plt.savefig("ar-eig.png", dpi=400)


# make a pandas Dataframe and write to csv file
df = pd.DataFrame(data_dict)
df.to_csv("buckle_data.csv")
