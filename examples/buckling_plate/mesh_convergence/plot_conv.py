import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import niceplots

plt.style.use(niceplots.get_style())
df = pd.read_csv("mesh_convergence.csv")
AR = df["AR"].to_numpy()
nelem = df["nelem"].to_numpy()
dt = df["dt[s]"].to_numpy()
lams = [df[f"lam{i}"].to_numpy() for i in range(5)]

for m_AR in [1.0, 3.0, 5.0]:
    fig = plt.figure(f"AR{m_AR}", figsize=(10, 7))
    ax = plt.subplot(111)
    plt.margins(x=0.03, y=0.03)
    mask = AR == m_AR
    for i in range(4, -1, -1):
        ax.plot(nelem[mask], lams[i][mask], "o-", label=r"$k_" + str(i + 1) + r"$")
    # plt.legend()
    plt.xscale("log")
    # plt.yscale("log")
    plt.xlabel("# elements")
    plt.ylabel("k, eigenvalue")
    plt.title(f"Aspect Ratio = {m_AR}")
    plt.ylim(bottom=0.0, top=np.max(lams[4][mask]) + 0.2)
    # if m_AR != 1.0:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(f"AR_{m_AR}_mesh_conv.png", dpi=400)
    plt.close(f"AR{m_AR}")

fig = plt.figure(f"AR{m_AR}-runtime", figsize=(10, 7))
ax = plt.subplot(111)
plt.margins(x=0.03, y=0.03)
for m_AR in [1.0, 3.0, 5.0]:
    mask = AR == m_AR
    ax.plot(nelem[mask], dt[mask], "o-", label=r"$\rho=" + str(round(m_AR)) + r"$")
# plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# elements")
plt.ylabel("runtime (sec)")
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig(f"mesh_conv_runtime.png", dpi=400)
plt.close(f"AR{m_AR}-runtime")
