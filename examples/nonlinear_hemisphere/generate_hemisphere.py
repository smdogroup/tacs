from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pprint import pprint
from scipy.spatial.transform import Rotation

parser = argparse.ArgumentParser()
parser.add_argument("--ne1", type=int, default=16, help="# elements from top to bottom")
parser.add_argument(
    "--ne2", type=int, default=16, help="# elements around circumference"
)
args = parser.parse_args()

# Overall dimensions
R = 10.0

# Number of elements along each edge of a single panel
ne1 = args.ne1
ne2 = args.ne2

np1 = 1
np2 = 1

# Nodes
n1 = ne1 + 1
n2 = ne2 + 1
theta1 = np.linspace(0, np.deg2rad(72.0), n1)
theta2 = np.linspace(0, np.deg2rad(90.0), n2)
Theta1, Theta2 = np.meshgrid(theta1, theta2)
X = np.zeros_like(Theta1)
Y = np.zeros_like(Theta1)
Z = np.zeros_like(Theta1)

# Node numbering
nid = np.zeros((n2, n1), dtype="intc")
bcnodes = []
count = 1
for i in range(n2):
    Rz = Rotation.from_euler("z", theta2[i])
    for j in range(n1):
        nid[i, j] = count
        # Compute node location on hemisphere
        r = R * np.array([np.cos(theta1[j]), 0.0, np.sin(theta1[j])])
        r = Rz.apply(r)
        X[i, j] = r[0]
        Y[i, j] = r[1]
        Z[i, j] = r[2]
        if j == 0:
            if i == 0:
                bcnodes.append({"nodenum": count, "fixedDOF": "2346"})
            else:
                bcnodes.append({"nodenum": count, "fixedDOF": "246"})
        elif j == n1 - 1:
            bcnodes.append({"nodenum": count, "fixedDOF": "156"})
        count += 1
nodes = np.stack((X, Y, Z), axis=2)
nmat = nodes.reshape((n1 * n2, 3))

# Connectivity
ne1 = n1 - 1
ne2 = n2 - 1
ne = ne1 * ne2
ncomp = 1
conn = {i + 1: [] for i in range(ncomp)}
ie = 1
for i in range(ne2):
    for j in range(ne1):
        compID = i // ne2 * np1 + j // ne1 + 1
        conn[compID].append(
            [ie, nid[i, j], nid[i + 1, j], nid[i + 1, j + 1], nid[i, j + 1]]
        )
        ie += 1


# Write BDF
output_file = "hemisphere.bdf"
fout = open(output_file, "w")


def write_80(line):
    newline = "{:80s}\n".format(line.strip("\n"))
    fout.write(newline)


write_80("SOL 103")
write_80("CEND")
write_80("BEGIN BULK")

# Make component names
compNames = {}
compID = 1
for i in range(np2):
    for j in range(np1):
        compNames[compID] = "PLATE.{:03d}/SEG.{:02d}".format(i, j)
        compID += 1


def write_bulk_line(key, items, format="small"):
    if format == "small":
        width = 8
        writekey = key
    elif format == "large":
        width = 16
        writekey = key + "*"
    line = "{:8s}".format(writekey)
    for item in items:
        if type(item) in [int, np.int64, np.int32]:
            line += "{:{width}d}".format(item, width=width)[:width]
        elif type(item) in [float, np.float64]:
            line += "{: {width}f}".format(item, width=width)[:width]
        elif type(item) is str:
            line += "{:{width}s}".format(item, width=width)[:width]
        else:
            print(type(item), item)
        if len(line) == 72:
            write_80(line)
            line = " " * 8
    if len(line) > 8:
        write_80(line)


# Write nodes
for i in range(n1 * n2):
    write_bulk_line("GRID", [i + 1, 0, nmat[i, 0], nmat[i, 1], nmat[i, 2], 0, 0, 0])

# Write elements
compID = 1
for key in conn:
    famPrefix = "$       Shell element data for family    "
    famString = "{}{:39s}".format(famPrefix, compNames[compID])
    write_80(famString)
    compID += 1
    for element in conn[key]:
        element.insert(1, key)
        write_bulk_line("CQUAD4", element)

# Write boundary conditions
for node in bcnodes:
    write_bulk_line("SPC", [1, node["nodenum"], node["fixedDOF"], 0.0])

write_80("ENDDATA")

fout.close()
