import numpy as np
import matplotlib.pyplot as plt

# Configure
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"

# Optionally set font to Computer Modern to avoid common missing font errors
params = {
    "axes.labelsize": 18,
    "legend.fontsize": 12,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "text.usetex": True,
}
plt.rcParams.update(params)

plt.rcParams["text.latex.preamble"] = [r"\usepackage{sfmath}"]
plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams["font.size"] = 12
plt.rcParams["font.weight"] = "bold"
plt.rcParams["lines.linewidth"] = 4
plt.rcParams["lines.color"] = "r"

# Make sure everything is within the frame
plt.rcParams.update({"figure.autolayout": True})

# Set marker size
markerSize = 5.0  # 11.0
mew = 1.0

# These are the "Tableau 20" colors as RGB.
tableau20 = [
    (31, 119, 180),
    (174, 199, 232),
    (255, 127, 14),
    (255, 187, 120),
    (44, 160, 44),
    (152, 223, 138),
    (214, 39, 40),
    (255, 152, 150),
    (148, 103, 189),
    (197, 176, 213),
    (140, 86, 75),
    (196, 156, 148),
    (227, 119, 194),
    (247, 182, 210),
    (127, 127, 127),
    (199, 199, 199),
    (188, 189, 34),
    (219, 219, 141),
    (23, 190, 207),
    (158, 218, 229),
]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255.0, g / 255.0, b / 255.0)

# Define colors for flap, torsion and leadlag
flap_color = tableau20[2]
lead_lag_color = tableau20[0]
torsion_color = tableau20[4]

# Data
num_freqs = 12
Omegas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
ydata = np.zeros((num_freqs, len(Omegas)))

# Read the data file
inpFile = open("freqdata.dat", "r")
freqdata = list(inpFile.readlines())
inpFile.close()

# Store the data in a format suitable for plotting
freqindex = 0
for k in xrange(len(Omegas)):
    freq = []
    for i in xrange(num_freqs):
        line = freqdata[freqindex * (num_freqs + 1) + i + 1]
        entry = line.split()
        freq.append(float(entry[0]))

    # Add frequencies to array
    ydata[:, k] = freq
    print(freq)
    freqindex += 1

# identify modes
flap_modes = [0, 2, 5, 6, 7]
torsion_modes = [1, 3, 9]
lead_lag_modes = [4, 8]

# Make plots
fig, ax = plt.subplots()
ax.spines["right"].set_visible(True)
ax.spines["top"].set_visible(True)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
plt.axis([0.085, 1.21, 0, 16])
xlabel = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
ylabel = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
plt.xticks(xlabel)
plt.yticks(ylabel)
plt.xlabel("$\Omega/\Omega_{ref}$")
plt.ylabel("$\omega/\Omega_{ref}$")

for k in flap_modes:
    plt.plot(
        Omegas, ydata[k, :], "-s", mew=mew, ms=markerSize, color=flap_color, mec="black"
    )

for k in lead_lag_modes:
    plt.plot(
        Omegas,
        ydata[k, :],
        "-s",
        mew=mew,
        ms=markerSize,
        color=lead_lag_color,
        mec="black",
    )

for k in torsion_modes:
    plt.plot(
        Omegas,
        ydata[k, :],
        "-s",
        mew=mew,
        ms=markerSize,
        color=torsion_color,
        mec="black",
    )

lines = ax.get_lines()
legend = plt.legend(
    [lines[0], lines[9], lines[5]],
    ["flap", "lag", "torsion"],
    loc="upper left",
    ncol=3,
    framealpha=0.5,
)
ax.add_artist(legend)
plt.savefig("fan.pdf", bbox_inches="tight", pad_inches=0.05)
